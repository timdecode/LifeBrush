// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"


#include "SwarmGenerator.h"
#include "Simulation/MeshFilamentSimulation.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "Simulation/FlexElements.h"
#include "MolecularLego/MolecularLego_Relaxation.h"
#include "Simulation/Brownian.h"

void USwarmGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::attach(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;
	_graph = &context->graph();

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	USwarmSimulation * swarmSim = _simulationManager->registerSimulation<USwarmSimulation>();
	UMeshFilamentSimulation * filamentSim = _simulationManager->registerSimulation<UMeshFilamentSimulation>();

	_simulationManager->registerSimulation<USingleParticleBrownianSimulation>();
	_simulationManager->attachSimulation<USingleParticleBrownianSimulation>();

	_simulationManager->registerSimulation<UStaticPositionSimulation>();
	_simulationManager->attachSimulation<UStaticPositionSimulation>();

	_simulationManager->attachAllSimulations();

	flexSimulation->play();

}

void USwarmGenerator::detach()
{
	Super::detach();

	flexSimulation->pause();
}


void USwarmGenerator::flexTick(
	float deltaT,
	NvFlexVector<int>& neighbourIndices,
	NvFlexVector<int>& neighbourCounts,
	NvFlexVector<int>& apiToInternal,
	NvFlexVector<int>& internalToAPI,
	int maxParticles)
{
}

FQuat USwarmSimulation::randomQuat()
{
	FVector axis = rand.GetUnitVector();
	float angle = rand.GetFraction() * M_PI * 2.0f;

	return FQuat(axis, angle);
}

void USwarmGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{

}

void USwarmGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	USwarmSimulation * swarmSimulation = _simulationManager->simulation<USwarmSimulation>();

	FVector localPoint = flexSimulationComponent->GetOwner()->GetTransform().InverseTransformPosition(point);

	if( swarmSimulation)
		swarmSimulation->addBrushPoint(localPoint, brushType);
}

void USwarmGenerator::endBrushPath()
{

}

void USwarmGenerator::setPrototype(FGraphNodeHandle exampleHandle)
{
	_brushExampleHandle = exampleHandle;

	USwarmSimulation * swarmSimulation = _simulationManager->simulation<USwarmSimulation>();
	swarmSimulation->setBrushPrototypoe(_brushExampleHandle);
}

void USwarmSimulation::_applyStar(FGraphNodeHandle starHandle)
{
	FGraphNode starNode = starHandle(*graph);

	FBoidStar * boidStar = graph->componentPtr<FBoidStar>(starHandle);

	FSGStarRule * starRule = ruleGraph.componentPtr<FSGStarRule>(boidStar->_exemplarPrototype);

	FBoidRootNode * boidRoot = graph->componentPtr<FBoidRootNode>(starHandle);
	FGraphNodeHandle rootNode = boidRoot ? boidRoot->rootNode : FGraphNodeHandle::null;

	if (!starRule) return;

	FVector forward = starNode.orientation.RotateVector(FVector::ForwardVector);

	TArray<FGraphNodeHandle> ring;

	for (int i = 0; i < starRule->branches.Num(); ++i)
	{
		FSGStarBranchEntry& branch = starRule->branches[i];
		FGraphNodeHandle prototypeHandle = starRule->handlesForBranches[i];

		if (!prototypeHandle) continue;

		FGraphNode prototypeNode = ruleGraph.node(prototypeHandle);

		FVector dir = starNode.orientation.RotateVector(branch.direction);

		FVector position = starNode.position + dir;
		FQuat rotation = FQuat::FindBetween(forward, dir) * prototypeNode.orientation;

		FGraphNodeHandle newHandle = copyElement(prototypeHandle, ruleGraph, *graph, position, rotation, rootNode);

		ring.Add(newHandle);

		FGraphNode& newNode = graph->node(newHandle);

		if (newNode.hasComponent<FBoidGenerator>())
		{
			FBoidGenerator& generator = newNode.component<FBoidGenerator>(*graph);

			generator.last = starHandle;
			generator.starNode = starHandle;

			generator.currentDistance = generator.segmentLength;
		}

		// shoot it in a direction
		FVelocityGraphObject& velocity = newNode.addComponent<FVelocityGraphObject>(*graph);
		velocity.linearVelocity = dir * 8.0f;
	}


}

void USwarmSimulation::_updateSpaceBVH()
{
	auto& storage = graph->componentStorage<FBoid>();

	for (FBoid& boid : storage)
	{
		auto particleIndex = boid.nodeHandle().index;

		FGraphNode& node = graph->node(boid.nodeHandle());

		if (!boid.isValid() && _spaceBVH.containsParticle(particleIndex))
		{
			_spaceBVH.removeParticle(particleIndex);
		}
		else
		{
			if (_spaceBVH.containsParticle(particleIndex))
				_spaceBVH.updateParticle(particleIndex, node.position, boid.radius);
			else
				_spaceBVH.insertParticle(particleIndex, node.position, boid.radius);
		}
	}
}

void USwarmSimulation::_updateRepellerBVH()
{
	auto& storage = graph->componentStorage<FBoidRepeller>();

	for (FBoidRepeller& boid : storage)
	{
		auto particleIndex = boid.nodeHandle().index;

		FGraphNode& node = graph->node(boid.nodeHandle());

		if (!boid.isValid() && _repellerBVH.containsParticle(particleIndex))
		{
			_repellerBVH.removeParticle(particleIndex);
		}
		else
		{
			if (_repellerBVH.containsParticle(particleIndex))
				_repellerBVH.updateParticle(particleIndex, node.position, 1.0f);
			else
				_repellerBVH.insertParticle(particleIndex, node.position, 1.0f);
		}
	}
}



















void USwarmSimulation::setBrushPrototypoe(FGraphNodeHandle handle)
{
	_brushPrototype = handle;
}

void USwarmSimulation::attach()
{
	_initRuleGraph();
	
	_readRuleGraph();
}

void USwarmSimulation::begin()
{
	if (debugDrawBoids)
		_attachSegmentMeshes();
	else
		_detachSegmentMeshes();

	rand.GenerateNewSeed();

	_updateSpaceBVH();
	_updateRepellerBVH();
}

void USwarmSimulation::detach()
{
	_detachSegmentMeshes();
}

void USwarmSimulation::_attachSegmentMeshes()
{
	auto& segmentStorage = graph->componentStorage<FBoidSegment>();
	auto& meshStorage = graph->componentStorage<FGraphMesh>();

	for (FBoidSegment& segment : segmentStorage)
	{
		FGraphMesh * mesh = meshStorage.componentPtrForNode(segment.nodeHandle());

		if( mesh )
			continue;

		mesh = &segment.nodeHandle()(graph).addComponent<FGraphMesh>(*graph);

		mesh->material = material_segment;
		mesh->staticMesh = staticMesh_segment;
	}
}

void USwarmSimulation::_detachSegmentMeshes()
{
	auto& segmentStorage = graph->componentStorage<FBoidSegment>();
	auto& meshStorage = graph->componentStorage<FGraphMesh>();

	for (FBoidSegment& segment : segmentStorage)
	{
		FGraphNode& node = segment.node(*graph);

		if( !node.hasComponent<FGraphMesh>() )
			continue;

		node.removeComponent<FGraphMesh>(*graph);
	}
}

void USwarmSimulation::_initRuleGraph()
{
	ruleGraph.init();
}

void USwarmSimulation::_readRuleGraph()
{
	_nameToRuleHandles.Empty();

	// load the name to rule handle table
	auto& ruleTypes = ruleGraph.componentStorage<FSGRuleTypeName>();

	for (FSGRuleTypeName& ruleType : ruleTypes)
	{
		_nameToRuleHandles.Add(ruleType.typeName, ruleType.nodeHandle());
	}

	// load our rule prototypes
	_readStars();
}

void USwarmSimulation::_readStars()
{
	auto& ruleTypes = ruleGraph.componentStorage<FSGRuleTypeName>();

	auto& starPrototypes = ruleGraph.componentStorage<FSGStarRule>();

	for (FSGStarRule& starPrototype : starPrototypes)
	{
		starPrototype.handlesForBranches.Empty();
		
		for (FSGStarBranchEntry branch : starPrototype.branches)
		{
			FGraphNodeHandle handle = _nameToRuleHandles.FindOrAdd(branch.typeName);

			starPrototype.handlesForBranches.Add(handle);
		}
	}

	// cache the sequences
	auto& filamentRules = ruleGraph.componentStorage<FSGFilamentRule>();

	for(FSGFilamentRule& rule : filamentRules)
	{
		rule.headHandles = _parseSequence(rule.headSequence);
		rule.bodyHandles = _parseSequence(rule.bodySequence);
		rule.tailHandles = _parseSequence(rule.tailSequence);
	}
}

TArray<FGraphNodeHandle> USwarmSimulation::_parseSequence(FString sequence)
{
	TArray<FString> headNames;
	sequence.ParseIntoArray(headNames, TEXT(","));

	TArray<FGraphNodeHandle> handles;

	for (FString& stringName : headNames)
	{
		FName name(*stringName);

		FGraphNodeHandle * headPrototype = _nameToRuleHandles.Find(name);

		if (!headPrototype) continue;

		handles.Add(*headPrototype);
	}

	return handles;
}

void USwarmSimulation::tick(float deltaT)
{
	graph->beginTransaction();

	_updateSpaceBVH();


	for (auto& brushPoint : _brushPoints)
	{
		bool hit = false;

		unrealAABB::AABB query(brushPoint.point, brushSpacing);

		_spaceBVH.query(query, [&hit](unsigned int particleIndex) {
			hit = true;

			return false; // don't continue
		});

		if (hit)
			continue;;


		if (_brushPrototype)
		{
			FVector point = brushPoint.point;
			FQuat rotation = randomQuat();

			FGraphNodeHandle handle = copyElement(_brushPrototype, ruleGraph, *graph, point, rotation, FGraphNodeHandle::null);

			if (FBoid * boid = graph->componentPtr<FBoid>(handle))
			{
				_spaceBVH.insertParticle(handle.index, point, boid->radius);
			}

		}
	}

	_brushPoints.clear();

	_tickFilamentAnchors(deltaT);

	_tickBoidStars(deltaT);

	_tickBoidGenerators();

	_tickBoidSeekers(deltaT);

	_updateRepellerBVH();
	_tickBoidRepellers(deltaT);

	graph->endTransaction();
}

FGraphNodeHandle USwarmSimulation::copyElement(
	FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph,
	FGraph& targetGraph,
	const FVector position,
	const FQuat rotation_in, 
	FGraphNodeHandle boidRootNode)
{
	FGraphNodeHandle result = FGraphCopyContext::copyAggregate(sourceNodeHandle, sourceGraph, targetGraph, position, rotation_in);

	// convert rules into agents
	// all rules also get a boid node so we can track them in the BVH
	FGraphNode& node = targetGraph.node(result);

	float defaultBoidRadius = radius_segment;

	auto addBoid = [&targetGraph, defaultBoidRadius](FGraphNode& node) {
		FBoid * boid = targetGraph.componentPtr<FBoid>(node.handle());

		if (!boid)
		{
			boid = &node.addComponent<FBoid>(targetGraph);
			boid->radius = defaultBoidRadius;
		}

		return boid;
	};

	// add a FBoidStar
	if (node.hasComponent<FSGStarRule>())
	{
		FBoidStar& boidStar = node.addComponent<FBoidStar>(targetGraph);

		boidStar._exemplarPrototype = sourceNodeHandle;

		node.removeComponent<FSGStarRule>(targetGraph);

		addBoid(node);
	}

	// add an FBoidGenerator
	if (FSGFilamentRule * filamentRule = targetGraph.componentPtr<FSGFilamentRule>(node.handle()))
	{
		FBoidGenerator& boidGenerator = node.addComponent<FBoidGenerator>(targetGraph);

		boidGenerator._exemplarPrototype = sourceNodeHandle;
		boidGenerator.filamentRadius = filamentRule->radius;
		boidGenerator.segmentLength = filamentRule->segmentLength;

		node.removeComponent<FSGFilamentRule>(targetGraph);

		addBoid(node);
	}

	if (FBoid * boid = targetGraph.componentPtr<FBoid>(node.handle()))
	{
		_spaceBVH.insertParticle(node.id, node.position, boid->radius);
	}

	// add the root node, and if there isn't one, we are the one.
	if (!boidRootNode)
		boidRootNode = result;

	node.addComponent<FBoidRootNode>(targetGraph, boidRootNode);

	return result;
}

void USwarmSimulation::addBrushPoint(FVector point, ESwarmGenerator_BrushType type)
{
	_brushPoints.emplace_back();

	auto& brushPoint = _brushPoints.back();

	brushPoint.point = point;
	brushPoint.type = type;
}

void USwarmSimulation::_tickFilamentAnchors(float deltaT)
{
	auto& anchors = graph->componentStorage<FSGFilamentAnchor>();
	auto& boidGenerators = graph->componentStorage<FBoidGenerator>();

	TArray<FGraphNodeHandle> toRemove;

	// instantiate SGFilaments from our anchors
	for (FSGFilamentAnchor& anchor : anchors)
	{
		if (!anchor.isValid()) continue;
		
		FGraphNodeHandle * prototype = _nameToRuleHandles.Find(anchor.filamentRuleName);

		if (prototype)
		{
			FGraphNode& anchorNode = graph->node(anchor.nodeHandle());

			FQuat rotation = anchorNode.orientation;

			FVector dir = rotation.RotateVector(FVector::ForwardVector);
			FVector up = rotation.RotateVector(FVector::UpVector);


			FVector position = anchorNode.position + dir;

			FGraphNodeHandle instantiatedHandle = copyElement(*prototype, ruleGraph, *graph, position, rotation, anchor.nodeHandle());

			if (FBoidGenerator * generator = boidGenerators.componentPtrForNode(instantiatedHandle))
			{
				generator->last = anchor.nodeHandle();

				// add a particle to the anchor if it doesn't have one
				FGraphNode& anchorNode = graph->node(anchor.nodeHandle());

				if (!anchorNode.hasComponent<FFlexParticleObject>())
					anchorNode.addComponent<FFlexParticleObject>(*graph);

				FGraphNode& generatorNode = graph->node(instantiatedHandle);

				FVelocityGraphObject& generatorVelocity = generatorNode.addComponent<FVelocityGraphObject>(*graph);
				generatorVelocity.linearVelocity = up * 5.0f;
			}
		}

		// we either instantiated the anchor, of there wasn't one. Time to die.
		toRemove.Add(anchor.nodeHandle());
	}

	// cleanup
	for (FGraphNodeHandle handle : toRemove)
	{
		graph->node(handle).removeComponent<FSGFilamentAnchor>(*graph);
	}
}

void USwarmSimulation::_tickBoidStars(float deltaT)
{
	auto& boidStars = graph->componentStorage< FBoidStar>();

	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	for (FBoidStar& star : boidStars)
	{
		if (!star.isValid())
			continue;

		if (!star._didApply)
		{
			_applyStar(star.nodeHandle());
			star._didApply = true;
		}

		// build the star ring, if it's ready (the generators have added the first elements)
		if (star._inProgressRing.Num() == star.targetRingSize && star.targetRingSize > 0)
		{
			// connect the nodes in the ring, so they point away from each other
			auto& ring = star._inProgressRing;

			for (int i = 0; i < ring.Num(); ++i)
			{
				FGraphNodeHandle cur = ring[i];
				FGraphNodeHandle next = ring[(i + 1) % ring.Num()];

				if (cur == next) continue;

				float distance = FVector::Dist(graph->node(cur).position, graph->node(next).position);

				auto oneRingEdge = graph->connectNodes<FFlexConnection>(cur, next, distance, 0.5f);
			}
		}
	}
}

void USwarmSimulation::_tickBoidRepellers(float deltaT)
{
	auto& boidRepellers = graph->componentStorage<FBoidRepeller>();

	auto& velocities = graph->componentStorage<FVelocityGraphObject>();


	for (FBoidRepeller& repeller : boidRepellers)
	{
		if (!repeller.isValid())
			continue;

		FGraphNode& repellerNode = graph->node(repeller.nodeHandle());

		for (FBoidRepelledPair& pair : repelledPairs)
		{
			if (repeller.species == pair.speciesA || repeller.species == pair.speciesB)
			{
				FName otherSpecies = repeller.species == pair.speciesA ? pair.speciesB : pair.speciesA;

				unrealAABB::AABB query(repellerNode.position, pair.distance);

				_repellerBVH.query(query, [&](unsigned int particleIndex) {
					FGraphNodeHandle otherHandle(particleIndex);

					FGraphNode& otherNode = graph->node(otherHandle);

					FVector dir = otherNode.position - repellerNode.position;

					dir = dir.GetSafeNormal();

					repellerNode.position -= dir * deltaT;
					otherNode.position += dir * deltaT;

					return true;
				});
			}
		}
	}
}

void USwarmSimulation::_tickBoidGenerators()
{
	auto& generatorStorage = graph->componentStorage<FBoidGenerator>();

	std::vector<FGraphNodeHandle> toRemove;

	

	UMeshFilamentSimulation * filamentSimulation = simulationManager->simulation<UMeshFilamentSimulation>();
	auto& filamentRules = ruleGraph.componentStorage<FSGFilamentRule>();

	auto& seekers = graph->componentStorage<FBoidSeeker>();

	// assign filament groups
	for (FBoidGenerator& boid : generatorStorage)
	{
		if (!boid.isValid())
			continue;

		if (boid.filamentGroup < 0)
		{
			FSGFilamentRule * filamentRule = filamentRules.componentPtrForNode(boid._exemplarPrototype);

			if (filamentRule)
				boid.filamentGroup = filamentSimulation->nextGroup(filamentRule->material);
			else
				boid.filamentGroup = filamentSimulation->nextGroup();
		}
	}

	// if the generator doesn't have a last, create one
	for (FBoidGenerator& boid : generatorStorage)
	{
		if (!boid.isValid() || boid.last )
			continue;

		FGraphNode& boidNode = graph->node(boid.nodeHandle());

		FVector p = boidNode.position - boidNode.orientation.RotateVector(FVector::ForwardVector);

		FGraphNodeHandle newHandle = FGraphNodeHandle(graph->addNode(boidNode.position, boidNode.orientation, boidNode.scale));

		FGraphNode& newNode = graph->node(newHandle);

		newNode.addComponent<FFlexParticleObject>(*graph);
		newNode.addComponent<FVelocityGraphObject>(*graph);

		boid.last = newHandle;
	}

	graph->beginTransaction();

	for (FBoidGenerator& boid : generatorStorage)
	{
		if (!boid.isValid() || !boid.last)
			continue;

		FGraphNode& boidNode = graph->node(boid.nodeHandle());

		FGraphNode& lastNode = graph->node(boid.last);

		FVector dir = (boidNode.position - lastNode.position);
		float distance = dir.Size();

		// are we far enough from the last guy to create a new segment?
		// and don't overlap with ourselves
		if (distance < boid.currentDistance + radius_segment * 2.0f)
			continue;

		dir = dir.GetSafeNormal();

		FVector segmentPosition = lastNode.position + dir * boid.currentDistance;

		FGraphNodeHandle theExampleHandle = _handleInFilamentSequence(boid.currentIndex, boid._exemplarPrototype);

		FGraphNodeHandle newHandle;

		// is this the last segment in the chain?
		if (!theExampleHandle)
		{
			// queue the node for removal
			toRemove.push_back(boid.nodeHandle());
		}
		else
		{
			FGraphNode& theExampleNode = ruleGraph.node(theExampleHandle);
			FVector forward = theExampleNode.orientation.RotateVector(FVector::ForwardVector);
			FQuat rotation = FQuat::FindBetween(forward, dir);

			FBoidRootNode * boidRoot = graph->componentPtr<FBoidRootNode>(boid.nodeHandle());
			FGraphNodeHandle rootNode = boidRoot ? boidRoot->rootNode : FGraphNodeHandle::null;

			newHandle = copyElement(theExampleHandle, ruleGraph, *graph, segmentPosition, rotation, rootNode);

			if (FBoidSeeker * seeker = seekers.componentPtrForNode(newHandle))
			{
				seeker->starNode = boid.starNode;
			}
		}

		

		auto oneRingEdge = graph->connectNodes<FFlexConnection>(boid.last, newHandle, boid.currentDistance, 0.5f);

		auto& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(oneRingEdge);
		{
			filament.radius = boid.filamentRadius;
			filament.group = boid.filamentGroup;
			filament.segmentID = boid.currentIndex;
		}

		if (boid.last2ring)
		{
			FVector lastLastPosition = graph->node(boid.last2ring).position;
			float lastLastDist = FVector::Dist(lastLastPosition, segmentPosition);

			graph->connectNodes<FFlexConnection>(boid.last2ring, newHandle, lastLastDist, 0.5f);
		}
		// this is the first chain in the filament, add it to the star ring, if we have one
		else
		{
			if (boid.starNode)
			{
				if (FBoidStar * star = graph->componentPtr<FBoidStar>(boid.starNode))
				{
					star->_inProgressRing.Add(newHandle);
				}
			}
		}

		boid.last2ring = boid.last;
		boid.last = newHandle;
		boid.currentIndex++;

		// update the current distance
		FGraphNodeHandle theNextHandle = _handleInFilamentSequence(boid.currentIndex, boid._exemplarPrototype);
		if (theNextHandle)
		{
			FBoidSegment * nextSegment = ruleGraph.componentPtr<FBoidSegment>(theNextHandle);
			FBoidSegment * segment = ruleGraph.componentPtr<FBoidSegment>(newHandle);

			boid.currentDistance = 0.0f;

			if (nextSegment && segment)
			{
				boid.currentDistance = nextSegment->radius && segment->radius;
			}
			else
				boid.currentDistance = boid.segmentLength;
		}
		// we're done, delete the boid
		else
		{
			// it's the last body in the segment, extend the filament
			filament.bExtension = boid.filamentRadius;

			// queue the node for removal
			toRemove.push_back(boid.nodeHandle());
		}

		// don't touch boidNode anymore, it could be invalid
	}

	graph->endTransaction();

	graph->beginTransaction();
	// kill the boids that are done
	for (auto handle : toRemove)
	{
//		handle(graph).removeComponent<FBoidGenerator>(*graph);
		graph->removeNode(handle);
	}

	graph->endTransaction();
}

FGraphNodeHandle USwarmSimulation::_handleInFilamentSequence(int32 sequenceIndex, FGraphNodeHandle ruleNode)
{
	FSGFilamentRule * filamentRule = ruleGraph.componentPtr<FSGFilamentRule>(ruleNode);

	FGraphNodeHandle result = FGraphNodeHandle::null;

	if (!filamentRule) result;

	if (sequenceIndex < filamentRule->headHandles.Num())
		result = filamentRule->headHandles[sequenceIndex];

	else if (sequenceIndex < filamentRule->headHandles.Num() + filamentRule->numInBody)
	{
		int32 index = (sequenceIndex - filamentRule->headHandles.Num()) % filamentRule->bodyHandles.Num();

		result = filamentRule->bodyHandles[index];
	}
	else if (sequenceIndex < filamentRule->headHandles.Num() + filamentRule->numInBody + filamentRule->tailHandles.Num())
	{
		int32 index = sequenceIndex - (filamentRule->headHandles.Num() + filamentRule->numInBody);

		result = filamentRule->tailHandles[index];
	}
	
	return result;
}

void USwarmSimulation::_tickBoidSeekers(float deltaT)
{
	BVH_t& mlElementsBVH = simulationManager->simulation<UMLElementSimulation>()->elementBVH;

	const float bindingDist = std::pow(radius_segment * 2.0f, 1.0);


	auto& seekers = graph->componentStorage<FBoidSeeker>();
	auto& mlElements = graph->componentStorage<FMLElement>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& blockers = graph->componentStorage<FBoidSeekerBlocker>();
	auto& roots = graph->componentStorage<FBoidRootNode>();

	graph->beginTransaction();

	for (FBoidSeeker& seeker : seekers)
	{
		FGraphNode& seekerNode = graph->node(seeker.nodeHandle());

		auto getRootHandle = [&roots](FGraphNodeHandle handle) {
			FBoidRootNode * rootNode = roots.componentPtrForNode(handle);

			return rootNode ? rootNode->rootNode : FGraphNodeHandle::null;
		};

		FGraphNodeHandle rootNodeHandle = getRootHandle(seeker.nodeHandle());

		unrealAABB::AABB query(seekerNode.position, seekerRadius);

		float minDistanceSqrd = std::numeric_limits<float>::max();
		FGraphNodeHandle minNodeHandle = FGraphNodeHandle::null;

		mlElementsBVH.query(query, [&](unsigned int particleIndex) {
			FGraphNodeHandle handle(particleIndex);

			// abort and keep searching
			if (handle == seeker.nodeHandle()) return true;

			FGraphNode& otherNode = graph->node(handle);

			FMLElement * otherElement = mlElements.componentPtrForNode(handle);
			FBoidSeeker * otherSeeker = seekers.componentPtrForNode(handle);


			FGraphNodeHandle otherNodeHandle = getRootHandle(handle);

			// abort if we are from the same root and keep searching
			if (rootNodeHandle == otherNodeHandle && rootNodeHandle) return true;

			// we might have added a block since we built the BVH
			if (blockers.componentPtrForNode(handle))
				return true;

			// abort if we don't match FMLSpeciesIDs and keep searching
			if ( seeker.bobSpecies.Find(otherElement->type) == INDEX_NONE ) return true;

			// find the nearest node
			float distSqrd = (otherNode.position - seekerNode.position).SizeSquared();

			if (distSqrd < minDistanceSqrd)
			{
				minDistanceSqrd = distSqrd;
				minNodeHandle = handle;
			}

			return true; // continue
		});

		// don't connect to blocked nodes (already bound)
		// don't use our node references, they can be invalidated in this block
		if (minNodeHandle && blockers.componentPtrForNode(minNodeHandle) == nullptr)
		{
			// attract to the neighbor
			FVelocityGraphObject * velocity = velocities.componentPtrForNode(seeker.nodeHandle());

			if (!velocity)
			{
				velocity = &seekerNode.addComponent<FVelocityGraphObject>(*graph);
			}

			FVector dir = (graph->node(minNodeHandle).position - seekerNode.position);
			float dist = dir.Size();
			dir = dir.GetSafeNormal();

			// set velocity (not n-body for now)
			velocity->linearVelocity = dir * 5.0f;

			// create a new node between them
			// and bind the seekers to it
			if (dist < bindingDist)
			{
				FBoidSeeker * otherSeeker = seekers.componentPtrForNode(minNodeHandle);

				if (otherSeeker)
					_connectSeekers(seekerNode.handle(), minNodeHandle, radius_segment * 2.0f);
				else
					_connectSeekerOneWay(seekerNode.handle(), minNodeHandle, radius_segment * 2.0f);
			}
		}
	}

	graph->endTransaction();
}

void USwarmSimulation::_connectSeekers(FGraphNodeHandle aHandle, FGraphNodeHandle bHandle, float desiredDistance)
{
	FVector dir = (graph->node(bHandle).position - graph->node(aHandle).position);
	float dist = dir.Size();
	dir = dir.GetSafeNormal();
	
	FVector p_mid = graph->node(aHandle).position + dir * dist * 0.5f;

	FGraphNodeHandle midHandle = FGraphNodeHandle(graph->addNode(p_mid));

	FGraphNode& aNode = graph->node(aHandle);
	FGraphNode& bNode = graph->node(bHandle);
	FGraphNode& midNode = graph->node(midHandle);


	

	midNode.addComponent<FFlexParticleObject>(*graph);

	if (!aNode.hasComponent<FFlexParticleObject>())
		aNode.addComponent<FFlexParticleObject>(*graph);

	if (!bNode.hasComponent<FFlexParticleObject>())
		bNode.addComponent<FFlexParticleObject>(*graph);

	auto aEdge = graph->connectNodes<FFlexConnection>(aHandle, midHandle, desiredDistance * 0.5f, 0.5f);
	auto bEdge = graph->connectNodes<FFlexConnection>(midHandle, bHandle, desiredDistance * 0.5f, 0.5f);
	
	auto getLastFilament = [&](FGraphNodeHandle handle) {
		FGraphNode& node = graph->node(handle);

		FFilamentConnection * theFilament = nullptr;

		node.each<FFilamentConnection>(*graph, [&](FGraphNodeHandle other, FFilamentConnection& filament)
		{
			theFilament = &filament;
		});

		return theFilament;
	};

	// extend the filaments
	if (FFilamentConnection * lastFilament = getLastFilament(aHandle))
	{
		FFilamentConnection& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(aEdge);

		filament.radius = lastFilament->radius;
		filament.group = lastFilament->group;
		filament.segmentID = lastFilament->segmentID <= 0 ? lastFilament->segmentID - 1 : lastFilament->segmentID + 1;
	}

	if (FFilamentConnection * lastFilament = getLastFilament(bHandle))
	{
		FFilamentConnection& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(bEdge);

		filament.radius = lastFilament->radius;
		filament.group = lastFilament->group;
		filament.segmentID = lastFilament->segmentID <= 0 ? lastFilament->segmentID - 1 : lastFilament->segmentID + 1;
	}

	_cleanupConnection(aNode, bNode);
}

void USwarmSimulation::_connectSeekerOneWay(FGraphNodeHandle seeker, FGraphNodeHandle target, float distance)
{
	auto getLastFilament = [&](FGraphNodeHandle handle) {
		FGraphNode& node = graph->node(handle);

		FFilamentConnection * theFilament = nullptr;

		node.each<FFilamentConnection>(*graph, [&](FGraphNodeHandle other, FFilamentConnection& filament)
		{
			theFilament = &filament;
		});

		return theFilament;
	};

	FGraphNode& seekerNode = graph->node(seeker);
	FGraphNode& targetNode = graph->node(target);

	if (!targetNode.hasComponent<FFlexParticleObject>())
		targetNode.addComponent<FFlexParticleObject>(*graph);

	auto edgeHandlef = graph->connectNodes<FFlexConnection>(seeker, target, distance, 0.5f);

	if (FFilamentConnection * lastFilament = getLastFilament(seeker))
	{
		FFilamentConnection& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(edgeHandlef);

		filament.radius = lastFilament->radius;
		filament.group = lastFilament->group;
		filament.segmentID = lastFilament->segmentID <= 0 ? lastFilament->segmentID - 1 : lastFilament->segmentID + 1;
	}

	_cleanupConnection(seekerNode, targetNode);
}

void USwarmSimulation::_cleanupConnection(FGraphNode& aNode, FGraphNode& bNode)
{
	if (aNode.hasComponent<FBoidSeeker>())
		aNode.removeComponent<FBoidSeeker>(*graph);

	if (bNode.hasComponent<FBoidSeeker>())
		bNode.removeComponent<FBoidSeeker>(*graph);

	//if (aNode.hasComponent<FRandomWalkGraphObject>())
	//	aNode.removeComponent<FRandomWalkGraphObject>(*graph);

	//if (bNode.hasComponent<FRandomWalkGraphObject>())
	//	bNode.removeComponent<FRandomWalkGraphObject>(*graph);

	//if (aNode.hasComponent<FSingleParticleBrownian>())
	//	aNode.removeComponent<FSingleParticleBrownian>(*graph);

	//if (bNode.hasComponent<FSingleParticleBrownian>())
	//	bNode.removeComponent<FSingleParticleBrownian>(*graph);

	aNode.addComponent<FBoidSeekerBlocker>(*graph);
	bNode.addComponent<FBoidSeekerBlocker>(*graph);
}


void ASGPrototypeActor::writeToElement(ElementTuple& element, FGraph& graph)
{
	Super::writeToElement(element, graph);

	handleInRuleGraph = element.handle();

	FSGRuleTypeName& ruleTypeName = element.node.addComponent<FSGRuleTypeName>(graph);

	ruleTypeName.typeName = typeName;
}
