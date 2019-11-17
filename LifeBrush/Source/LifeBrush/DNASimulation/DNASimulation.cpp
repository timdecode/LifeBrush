// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Simulation/Aggregates.h"
#include "MolecularLego/MolecularLego_Relaxation.h"
#include "Simulation/Brownian.h"

#include "DNASimulation/DNASimulation.h"
#include "Simulation/MeshFilamentSimulation.h"

const FName FRNACodon::StartCodon = FName(TEXT("AUG"));
const FName FRNACodon::ACA = FName(TEXT("ACA"));



void UDNASimulation::tick(float deltaT)
{
	// Ribosome stuff
	_initTRNAs();
	_tick_tRNAPayloadSlots(deltaT);

	_initRibosomes();
	_tickRibosomes(deltaT);

	_tickSlots_ribosomeTRNA(deltaT);
	_tickSlots_ribosomeMRNA(deltaT);
	_tickSlots_ribosomeAA(deltaT);

	// RNA Polymerase stuff
	_initRNAPolymerase();

	_tickSlots_RNAPolymerase_DNA(deltaT);
	_tickSlots_RNAPolymerase_RNA(deltaT);

	_tickPolymerase(deltaT);
}

void UDNASimulation::attach()
{
}

void UDNASimulation::detach()
{
}

void UDNASimulation::_initRibosomes()
{
	auto& ribosomes = graph->componentStorage<FRibosomeObject>();

	auto& mRNASlots = graph->componentStorage<FRibosomeMRNASlot>();
	auto& aaSlots = graph->componentStorage<FRibosomeAASlot>();
	auto& tRNASlots = graph->componentStorage<FRibosomeTRNASlot>();


	for (FRibosomeObject& ribosome : ribosomes)
	{
		if( ribosome.slot_mRNA_entrance ) continue;

		FMLAggregateNO * aggregate = FMLAggregateNO::getAggregate(*graph, ribosome.nodeHandle());

		if( !aggregate) continue;

		TArray<FGraphNodeHandle> nodes;
		TArray<FGraphEdgeHandle> edges;
		
		aggregate->edgesAndNodesInAggregate(*graph, nodes, edges);

		for (FGraphNodeHandle other : nodes)
		{
			if (FRibosomeMRNASlot * slot = mRNASlots.componentPtrForNode(other))
			{
				if (slot->slotType == ERibosomeSlotType::entrance)
					ribosome.slot_mRNA_entrance = other;
				else if (slot->slotType == ERibosomeSlotType::current)
					ribosome.slot_mRNA_cur = other;
				else if (slot->slotType == ERibosomeSlotType::exit)
					ribosome.slot_mRNA_exit = other;
			}

			if (FRibosomeAASlot * slot = aaSlots.componentPtrForNode(other))
			{
				ribosome.slot_aa_exit = other;
			}

			if (FRibosomeTRNASlot * slot = tRNASlots.componentPtrForNode(other))
			{
				if (slot->slotType == ERibosomeSlotType::entrance)
					ribosome.slot_tRNA_entrance = other;
				else if (slot->slotType == ERibosomeSlotType::current)
					ribosome.slot_tRNA_cur = other;
				else if (slot->slotType == ERibosomeSlotType::exit)
					ribosome.slot_tRNA_exit = other;
			}
		}
	}
}

void UDNASimulation::_tickRibosomes(float dt)
{
	typedef UMLElementSimulation::BVH_t BVH_T;
	
	auto& ribosomes = graph->componentStorage<FRibosomeObject>();
	auto& elements = graph->componentStorage<FMLElement>();
	auto& codons = graph->componentStorage<FRNACodon>();

	auto& mRNASlots = graph->componentStorage<FRibosomeMRNASlot>();
	auto& tRNASlots = graph->componentStorage<FRibosomeTRNASlot>();
	auto& aaSlots = graph->componentStorage<FRibosomeAASlot>();

	auto& tRNAs = graph->componentStorage<FtRNA>();

	BVH_T& bvh = simulationManager->simulation<UMLElementSimulation>()->elementBVH;

	UTimelineSimulation * timelineSimulation = simulationManager->simulation<UTimelineSimulation>();


	for (FRibosomeObject& ribosome : ribosomes)
	{
		// we didn't initialize for some reason
		if (!ribosome.didInit()) continue;

		ribosome.time -= dt;

		if( ribosome.time > 0.0f ) continue;

		// reset it
		ribosome.time = advancementRate;

		FRibosomeMRNASlot * mRNA_entranceSlot = mRNASlots.componentPtrForNode(ribosome.slot_mRNA_entrance);
		FRibosomeMRNASlot * mRNA_curSlot = mRNASlots.componentPtrForNode(ribosome.slot_mRNA_cur);
		FRibosomeMRNASlot * mRNA_exitSlot = mRNASlots.componentPtrForNode(ribosome.slot_mRNA_exit);

		FRibosomeTRNASlot * tRNA_entranceSlot = tRNASlots.componentPtrForNode(ribosome.slot_tRNA_entrance);
		FRibosomeTRNASlot * tRNA_curSlot = tRNASlots.componentPtrForNode(ribosome.slot_tRNA_cur);
		FRibosomeTRNASlot * tRNA_exitSlott = tRNASlots.componentPtrForNode(ribosome.slot_tRNA_exit);

		FRibosomeAASlot * aaSlot = aaSlots.componentPtrForNode(ribosome.slot_aa_exit);

		auto getTRNAPayloadSlot = [&](FGraphNodeHandle tRNAHandle) -> FtRNAPayloadSlot*
		{
			if (!tRNAHandle) return nullptr;

			FtRNA * tRNA = graph->componentPtr<FtRNA>(tRNAHandle);

			if (!tRNA) return nullptr;

			FtRNAPayloadSlot * tRNAPayloadSlot = graph->componentPtr<FtRNAPayloadSlot>(tRNA->payloadSlot);

			return tRNAPayloadSlot;
		};

		// build the amino acid chain
		if (tRNA_curSlot->occupied())
		{
			FtRNAPayloadSlot * tRNAPayloadSlot = getTRNAPayloadSlot(tRNA_curSlot->pigeon());

			// unbind the payload and grow the chain
			if (tRNAPayloadSlot->occupied())
			{
				_attachAminoAcidToPeptideChain(tRNAPayloadSlot->unbind(*graph), aaSlot);

				// don't seek until we're ejected
				tRNAPayloadSlot->seek = false;
			}
		}
		// we've finished the chain, let it loose
		else
		{
			aaSlot->unbind(*graph);
		}

		// lambda to step advance translation (see below for call)
		auto advance = [&]() {
			if (mRNA_exitSlot->occupied())
			{
				FGraphNodeHandle tRNAHandle = mRNA_exitSlot->pigeon();

				FVector p = graph->node(mRNA_exitSlot->nodeHandle()).position;

				FMLAggregateNO * aggregateNO = FMLAggregateNO::getAggregate(*graph, ribosome.nodeHandle());

				FGraphNodeHandle riboseomeRootHandle = aggregateNO ? aggregateNO->nodeHandle() : ribosome.nodeHandle();


				UEvent_tRNAEjected& pumpEvent = timelineSimulation->recordEvent<UEvent_tRNAEjected>(p);
				pumpEvent.triggeringAgent = riboseomeRootHandle;
				pumpEvent.otherAgents.Add(tRNAHandle);

				// add all the nodes with meshes
				if (aggregateNO)
				{
					// the mesh isn't on the root node, it's on another one
					FGraphNode& aggregateNode = graph->node(aggregateNO->nodeHandle());

					auto& meshes = graph->componentStorage<FGraphMesh>();

					aggregateNode.each<FMLAggregateEO>(*graph, [&](FGraphNodeHandle other, FMLAggregateEO& eo) {
						if (FGraphMesh * mesh = meshes.componentPtrForNode(other))
							pumpEvent.otherAgents.Add(mesh->nodeHandle());
					});
				}


			}


			mRNA_exitSlot->unbind(*graph);
			mRNA_exitSlot->bind(mRNA_curSlot->unbind(*graph), *graph);
			mRNA_curSlot->bind(mRNA_entranceSlot->unbind(*graph), *graph);


			FGraphNodeHandle toEnable = tRNA_exitSlott->unbind(*graph);
			tRNA_exitSlott->bind(tRNA_curSlot->unbind(*graph), *graph);
			tRNA_curSlot->bind(tRNA_entranceSlot->unbind(*graph), *graph);

			if (toEnable)
			{
				FtRNAPayloadSlot * tRNASlot = getTRNAPayloadSlot(toEnable);

				tRNASlot->seek = true;
			}
		};

		// advance when we have a tRNA and an mRNA in the entrance
		if (mRNA_entranceSlot->occupied() && tRNA_entranceSlot->occupied())
		{
			// cache the pigeon (We're going to unbind it)
			FGraphNodeHandle entrance_mRNA = mRNA_entranceSlot->pigeon();

			advance();

			// find the next slot in the filament
			mRNA_entranceSlot->bind(_nextInFilament(entrance_mRNA), *graph);
		}
		// end of sequence, flush
		else if (!mRNA_entranceSlot->occupied())
		{
			advance();
		}

		tRNA_entranceSlot->ignore0 = tRNA_curSlot->pigeon();
		tRNA_entranceSlot->ignore1 = tRNA_exitSlott->pigeon();

		// only mRNA_entrance should seek
		mRNA_entranceSlot->seek = true;
		mRNA_curSlot->seek = false;
		mRNA_exitSlot->seek = false;

		// seek only when we an mRNA in the entrance
		tRNA_entranceSlot->seek = mRNA_entranceSlot->occupied();
		tRNA_curSlot->seek = false;
		tRNA_exitSlott->seek = false;

		aaSlot->seek = false;
	}
}

void UDNASimulation::_initTRNAs()
{
	auto& tRNAs = graph->componentStorage<FtRNA>();
	auto& slots = graph->componentStorage<FtRNAPayloadSlot>();

	// find our payloads
	for (FtRNA& trna : tRNAs)
	{
		if( trna.payloadSlot ) continue;

		FGraphNode& node = graph->node(trna.nodeHandle());

		for (auto ei : node.edges)
		{
			FGraphEdge& edge = graph->edge(FGraphEdgeHandle(ei));

			FGraphNodeHandle other = edge.other(node.handle());

			FtRNAPayloadSlot * payloadSlot = slots.componentPtrForNode(other);

			if (payloadSlot)
			{
				payloadSlot->seek = true;
				trna.payloadSlot = other;

				break;
			}
		}
	}
}

void UDNASimulation::_tick_tRNAPayloadSlots(float dt)
{
	auto& aminoAcids = graph->componentStorage<FAminoAcid>();

	UTimelineSimulation * timelineSimulation = simulationManager->simulation<UTimelineSimulation>();

	FParticleSlot::tickSlots_withCallback<FtRNAPayloadSlot>(dt, [&](FtRNAPayloadSlot& slot, FGraphNodeHandle overlapHandle) {
		FAminoAcid * aminoAcid = aminoAcids.componentPtrForNode(overlapHandle);

		if (aminoAcid && aminoAcid->isValid()) return true;

		return false;
	}, this, ribosomeSeekRadius, [&](FtRNAPayloadSlot& slot, FGraphNodeHandle pigeonHandle) {
		FVector p = graph->node(slot.nodeHandle()).position;

		FMLAggregateNO * aggregateNO = FMLAggregateNO::getAggregate(*graph, slot.nodeHandle());

		FGraphNodeHandle slotRootHandle = aggregateNO ? aggregateNO->nodeHandle() : slot.nodeHandle();


		UEvent_tRNA_pickupAminoAcid& pumpEvent = timelineSimulation->recordEvent<UEvent_tRNA_pickupAminoAcid>(p);
		pumpEvent.triggeringAgent = slotRootHandle;
		pumpEvent.otherAgents.Add(pigeonHandle);



	});
}

void UDNASimulation::_tickSlots_ribosomeMRNA(float dt)
{
	auto& codons = graph->componentStorage<FRNACodon>();

	FParticleSlot::tickSlots<FRibosomeMRNASlot>(dt, [&](FRibosomeMRNASlot& slot, FGraphNodeHandle overlapHandle) {
		FRNACodon * codon = codons.componentPtrForNode(overlapHandle);

		if (!codon) return false;

		bool isStart = codon->codon == FRNACodon::StartCodon;

		if (slot.slotType == ERibosomeSlotType::entrance && isStart)
			return true;
		else
			return false;

		// it's a valid slot
		return true;
	}, this, ribosomeSeekRadius);
}

void UDNASimulation::_tickSlots_ribosomeTRNA(float dt)
{
	auto& tRNAs = graph->componentStorage<FtRNA>();

	// seek
	FParticleSlot::tickSlots<FRibosomeTRNASlot>(dt, [&](FRibosomeTRNASlot& slot, FGraphNodeHandle overlapHandle) {
		FtRNA * tRNA = tRNAs.componentPtrForNode(overlapHandle);

		if (!tRNA) return false;

		if (overlapHandle == slot.ignore0 || overlapHandle == slot.ignore1) return false;

		// for testing, we have a payload
		return true;
		//return tRNA->payload != FGraphNodeHandle::null;
	}, this, ribosomeSeekRadius);
}

void UDNASimulation::_tickSlots_ribosomeAA(float dt)
{
	// don't attach to anything
	FParticleSlot::tickSlots<FRibosomeAASlot>(dt, [&](FRibosomeAASlot& slot, FGraphNodeHandle overlapHandle) {
		return false;
	}, this, ribosomeSeekRadius);
}

void UDNASimulation::_attachAminoAcidToPeptideChain(FGraphNodeHandle newAminoAcidHandle, FRibosomeAASlot * aminoAcidSlot)
{
	checkf(newAminoAcidHandle, TEXT("We should have an amino acid handle"));

	FGraphNodeHandle peptideHandle = aminoAcidSlot->pigeon();

	// connect the chain
	if (peptideHandle)
	{
		FGraphNode & previousAminoAcidNode = graph->node(peptideHandle);

		{
			FGraphEdgeHandle filamentEdge = graph->connectNodes(peptideHandle, newAminoAcidHandle);

			FFilamentConnection& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(filamentEdge);
			filament.group = aminoAcidSlot->group;
			filament.segmentID = aminoAcidSlot->segmentID;
			filament.radius = 0.5f;

			FFlexConnection& flexConnection = graph->addOrReplaceEdgeObject<FFlexConnection>(filamentEdge);
			flexConnection.coefficient = 0.5f;
			flexConnection.length = aminoAcidSegmentLength;
		}

		// doubly linked for stiffness
		if (aminoAcidSlot->lastPigeon)
		{
			FGraphEdgeHandle doubleLink = graph->connectNodes(aminoAcidSlot->lastPigeon, newAminoAcidHandle);

			FFlexConnection& flexConnection = graph->addOrReplaceEdgeObject<FFlexConnection>(doubleLink);
			flexConnection.coefficient = 0.5f;
			flexConnection.length = aminoAcidSegmentLength * 2.0f;
		}
		
		aminoAcidSlot->segmentID++;

		FGraphNode & slotNode = graph->node(aminoAcidSlot->nodeHandle());

		FVector dir = FVector::ForwardVector * aminoAcidSegmentLength;
		dir = slotNode.orientation.RotateVector(dir);
		previousAminoAcidNode.position = slotNode.position + dir;

		if (FVelocityGraphObject * velocity = graph->componentPtr<FVelocityGraphObject>(peptideHandle))
			velocity->linearVelocity = dir * 5.0f;
	}
	else
	{
		aminoAcidSlot->group = simulationManager->simulation<UMeshFilamentSimulation>()->nextGroup(aminoAcidMaterial);
		aminoAcidSlot->segmentID = 0;
		aminoAcidSlot->lastPigeon = FGraphNodeHandle::null;
	}

	// attach it to the slot
	aminoAcidSlot->bind(newAminoAcidHandle, *graph);

	// we're a peptide now, so get rid of our amino acid and mesh
	FGraphNode & newAminoAcidNode = graph->node(newAminoAcidHandle);

	if( newAminoAcidNode.hasComponent<FAminoAcid>() )
		newAminoAcidNode.removeComponent<FAminoAcid>(*graph);

	// and we don't need a mesh
	if (newAminoAcidNode.hasComponent<FGraphMesh>())
		newAminoAcidNode.removeComponent<FGraphMesh>(*graph);
}

FGraphNodeHandle UDNASimulation::_nextInFilament(FGraphNodeHandle handle)
{
	FGraphNode& node = graph->node(handle);

	for (int32 ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		FFilamentConnection * filament = graph->edgeObjectPtr<FFilamentConnection>(edgeHandle);

		if (filament)
		{
			FGraphEdge edge = graph->edge(edgeHandle);

			// the filament is direction, from a to b
			if (FGraphNodeHandle(edge.b) != handle)
				return FGraphNodeHandle(edge.b);
		}
	}
	
	return FGraphNodeHandle::null;
}


bool FRibosomeObject::didInit()
{
	return (
		slot_mRNA_cur && slot_mRNA_entrance && slot_mRNA_exit &&
		slot_tRNA_cur && slot_tRNA_entrance && slot_tRNA_exit && 
		slot_aa_exit);
}

void FParticleSlot::bind(FGraphNodeHandle target, FGraph& graph)
{
	// unbind ourselves first
	unbind(graph);

	if (!target) return;

	// can't bind if the target is bound
	if (graph.componentPtr<FBoundToParticleSlot>(target))
		return;

	// bind it
	graph.node(target).addComponent<FBoundToParticleSlot>(graph);

	_pigeon = target;
}

FGraphNodeHandle FParticleSlot::unbind(FGraph& graph)
{
	FGraphNodeHandle result = _pigeon;

	if (_pigeon && graph.node(_pigeon).hasComponent<FBoundToParticleSlot>())
		graph.node(_pigeon).removeComponent<FBoundToParticleSlot>(graph);

	_pigeon = FGraphNodeHandle::null;

	return result;
}

FGraphNodeHandle FRibosomeAASlot::unbind(FGraph& graph)
{
	FGraphNodeHandle result = Super::unbind(graph);

	lastPigeon = FGraphNodeHandle::null;

	return result;
}

void FRibosomeAASlot::bind(FGraphNodeHandle target, FGraph& graph)
{
	lastPigeon = _pigeon;

	Super::bind(target, graph);
}









bool FRNAPolymeraseObject::didInit()
{
	return slot_mRNA_exit && slot_DNA_exit && slot_DNA_entrance && slot_DNA_cur;
}


void UDNASimulation::_initRNAPolymerase()
{
	auto& polymerases = graph->componentStorage<FRNAPolymeraseObject>();

	auto& dnaSlots = graph->componentStorage<FRNAPolymerase_DNASlot>();
	auto& mrnaSlots = graph->componentStorage<FRNAPolymerase_RNASlot>();

	for (FRNAPolymeraseObject& polymerase : polymerases)
	{
		if (polymerase.didInit()) continue;

		FMLAggregateNO * aggregate = FMLAggregateNO::getAggregate(*graph, polymerase.nodeHandle());

		if (!aggregate) continue;

		TArray<FGraphNodeHandle> nodes;
		TArray<FGraphEdgeHandle> edges;

		aggregate->edgesAndNodesInAggregate(*graph, nodes, edges);

		for (FGraphNodeHandle other : nodes)
		{
			if (FRNAPolymerase_DNASlot * slot = dnaSlots.componentPtrForNode(other))
			{
				if (slot->slotType == ERibosomeSlotType::entrance)
					polymerase.slot_DNA_entrance = other;
				else if (slot->slotType == ERibosomeSlotType::current)
					polymerase.slot_DNA_cur = other;
				else if (slot->slotType == ERibosomeSlotType::exit)
					polymerase.slot_DNA_exit = other;
			}

			if (FRNAPolymerase_RNASlot * slot = mrnaSlots.componentPtrForNode(other))
			{
				polymerase.slot_mRNA_exit = other;
			}
		}
	}
}

void UDNASimulation::_tickPolymerase(float dt)
{
	typedef UMLElementSimulation::BVH_t BVH_T;

	auto& polymerases = graph->componentStorage<FRNAPolymeraseObject>();
	auto& dnaSlots = graph->componentStorage<FRNAPolymerase_DNASlot>();
	auto& rnaSlots = graph->componentStorage<FRNAPolymerase_RNASlot>();
	auto& dnas = graph->componentStorage<FDNAObject>();

	for (FRNAPolymeraseObject& polymerase : polymerases)
	{
		// we didn't initialize for some reason?
		if(!polymerase.didInit()) continue;

		// tick the timer
		polymerase.time -= dt;

		if (polymerase.time > 0.0f) continue;

		// times up, reset it and do a tick on the polymerase
		polymerase.time = advancementRate;

		// RNA Polymerase moves in a 3' to 5' direction

		FRNAPolymerase_DNASlot * dna_fivePrime_start = dnaSlots.componentPtrForNode(polymerase.slot_DNA_entrance);
		FRNAPolymerase_DNASlot * dna_fivePrime_cur = dnaSlots.componentPtrForNode(polymerase.slot_DNA_cur);
		FRNAPolymerase_DNASlot * dna_fivePrime_exit = dnaSlots.componentPtrForNode(polymerase.slot_DNA_exit);

		FRNAPolymerase_RNASlot * rna_slot = rnaSlots.componentPtrForNode(polymerase.slot_mRNA_exit);

		// find a TATA
		if (!dna_fivePrime_start->occupied())
		{
			FGraphNodeHandle tata = _findTATA(polymerase.slot_DNA_entrance);

			if (tata)
			{
				FGraphNodeHandle next = _nextInFilament(tata);

				dna_fivePrime_start->bind(next, *graph);
			}
		}

		// advance
		FGraphNodeHandle next = dna_fivePrime_start->occupied() ?
			_findNextDNA(dna_fivePrime_start->pigeon()) :
			FGraphNodeHandle::null;

		dna_fivePrime_exit->bind(dna_fivePrime_cur->unbind(*graph), *graph);
		dna_fivePrime_cur->bind(dna_fivePrime_start->unbind(*graph), *graph);
		dna_fivePrime_start->bind(next, *graph);

		// we hit the end, start seeking again
		dna_fivePrime_start->seek = next == FGraphNodeHandle::null;

		// attach/detach
		if (dna_fivePrime_cur->occupied())
		{
			FDNAObject * dna = dnas.componentPtrForNode(dna_fivePrime_cur->pigeon());

			if (dna)
				_attachNucleotideToRNAFilament(*dna, *rna_slot);
		}
		else
		{
			rna_slot->unbind(*graph);

			// clear it
			rna_slot->clear();
		}
	}
}

void UDNASimulation::_tickSlots_RNAPolymerase_DNA(float dt)
{
	auto& dnas = graph->componentStorage<FDNAObject>();

	FParticleSlot::tickSlots<FRNAPolymerase_DNASlot>(dt, [&](FRNAPolymerase_DNASlot& slot, FGraphNodeHandle overlapHandle) {
		FDNAObject * dna = dnas.componentPtrForNode(overlapHandle);

		if (!dna) return false;

		// it's a valid slot
		return true;
	}, this, ribosomeSeekRadius);
}

void UDNASimulation::_tickSlots_RNAPolymerase_RNA(float dt)
{
	FParticleSlot::tickSlots<FRNAPolymerase_RNASlot>(dt, [&](FRNAPolymerase_RNASlot& slot, FGraphNodeHandle overlapHandle) {
		return true;
	}, this, ribosomeSeekRadius);
}

FGraphNodeHandle UDNASimulation::_findNextDNA(FGraphNodeHandle start)
{
	return _nextInFilament(start);
}

FGraphNodeHandle UDNASimulation::_findTATA(FGraphNodeHandle handle)
{
	typedef UMLElementSimulation::BVH_t BVH_T;

	auto& tatas = graph->componentStorage<FTATAObject>();
	auto& elements = graph->componentStorage<FMLElement>();

	BVH_T& bvh = simulationManager->simulation<UMLElementSimulation>()->elementBVH;

	FGraphNode node = graph->node(handle);

	FGraphNodeHandle nearestHandle = FGraphNodeHandle::null;
	float minDistSqrd = std::numeric_limits<float>::max();

	unrealAABB::AABB query(node.position, ribosomeSeekRadius);

	bvh.query(query, [&](unsigned int elementIndex) {
		FGraphNodeHandle particleHandle(elementIndex);

		if (!tatas.componentPtrForNode(particleHandle)) return true;

		const FVector position = graph->node(particleHandle).position;

		float distanceSqrd = FVector::DistSquared(position, node.position);

		if (distanceSqrd < minDistSqrd)
		{
			minDistSqrd = distanceSqrd;
			nearestHandle = particleHandle;
		}

		return true; // keep querying
	});

	return nearestHandle;
}

void UDNASimulation::_attachNucleotideToRNAFilament(FDNAObject& dna, FRNAPolymerase_RNASlot& rnaSlot)
{

	FMLElement * element = graph->componentPtr<FMLElement>(dna.nodeHandle());

	if (!element)
		return;

	rnaSlot.unfinishedCodon.Add(element->type);

	FGraphNode slotNode = graph->node(rnaSlot.nodeHandle());

	if (rnaSlot.unfinishedCodon.Num() < 3) return;

	// add a segment to the DNA
	FGraphNodeHandle newHandle(graph->addNode(slotNode.position, slotNode.orientation));
	FGraphNode& node = graph->node(newHandle);

	{
		FFlexParticleObject& particle = node.addComponent<FFlexParticleObject>(*graph);
		FVelocityGraphObject& velocity = node.addComponent<FVelocityGraphObject>(*graph);
		FSingleParticleBrownian& brownian = node.addComponent<FSingleParticleBrownian>(*graph);

		FRNACodon& codon = node.addComponent<FRNACodon>(*graph);

		if (rnaSlot.occupied())
		{
			// this is a hack, we should really translate this the unfinishedCodon array into a codon
			codon.codon = FRNACodon::ACA;
		}
		else
		{
			codon.codon = FRNACodon::StartCodon;

			rnaSlot.segmentGroup = simulationManager->simulation<UMeshFilamentSimulation>()->nextGroup(rnaMaterial);
			rnaSlot.segmentID = 0;
			rnaSlot.lastHandle = FGraphNodeHandle::null;
		}

		FMLElement& element = node.addComponent<FMLElement>(*graph);

		element.type = codon.codon;

		rnaSlot.unfinishedCodon.Empty();
	}

	FGraphNodeHandle last = rnaSlot.pigeon();

	if (last)
	{
		FVector dir = slotNode.orientation.RotateVector(FVector::ForwardVector);
		dir *= rnaSegmentLength;

		// send out the last node
		FGraphNode& lastNode = graph->node(last);
		lastNode.position = slotNode.position + dir;
		lastNode.orientation = slotNode.orientation;

		// create connections
		{
			FGraphEdgeHandle edgeHandle = graph->connectNodes(last, newHandle);

			FFilamentConnection& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(edgeHandle);
			filament.group = rnaSlot.segmentGroup;
			filament.segmentID = rnaSlot.segmentID;
			filament.radius = rnaSegmentRadius;

			FFlexConnection& flexConnection = graph->addOrReplaceEdgeObject<FFlexConnection>(edgeHandle);
			flexConnection.coefficient = 0.5f;
			flexConnection.length = rnaSegmentLength;

			rnaSlot.segmentID++;
		}
		
		// doubly connect the filament
		if (FGraphNodeHandle lastLast = rnaSlot.lastHandle)
		{
			FGraphEdgeHandle edgeHandle = graph->connectNodes(lastLast, newHandle);

			FFlexConnection& flexConnection = graph->addOrReplaceEdgeObject<FFlexConnection>(edgeHandle);
			flexConnection.coefficient = 0.5f;
			flexConnection.length = rnaSegmentLength * 2.0f;
		}
	}

	// bind to the slot
	rnaSlot.lastHandle = rnaSlot.pigeon();
	rnaSlot.bind(newHandle, *graph);
}

FGraphNodeHandle FTATAObject::find5Prime(FGraph& graph)
{
	FGraphNode& node = graph.node(nodeHandle());

	auto flexConnections = graph.edgeStorage<FFlexConnection>();

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle eh(ei);

		FGraphEdge& edge = graph.edge(eh);

		FFlexConnection * connection = flexConnections.objectPtr(eh);
		
		if( !connection ) continue;
		if( !flexConnections.isValid(eh) ) continue;

		if (FGraphNodeHandle(edge.a) == node.handle()) return FGraphNodeHandle(edge.b);
	}

	return FGraphNodeHandle::null;
}

void FRNAPolymerase_RNASlot::clear()
{
	segmentID = 0;
	lastHandle = FGraphNodeHandle::null;
	unfinishedCodon.Empty();
}
