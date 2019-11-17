// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"
#include "Simulation/Aggregates.h"
#include "Simulation/Brownian.h"

#include "Synapse/Synapse.h"

void USynapseSimulation::tick(float deltaT)
{
	
}

void USynapseSimulation::componentAdded(FGraphNodeHandle handle, ComponentType type)
{

}

void USynapseSimulation::componentRemoved(FGraphNodeHandle handle, ComponentType type)
{
	static ComponentType TFAcetylcholineReceptor = componentType<FAcetylcholineReceptor>();
	static ComponentType TFAChSlot = componentType<FAChSlot>();
	static ComponentType TFSodiumPotassiumPump = componentType<FSodiumPotassiumPump>();
	static ComponentType TFVoltageGatedIonChannel = componentType<FVoltageGatedIonChannel>();

	if (type == TFAcetylcholineReceptor)
	{
		if( _receptorBVH.containsParticle(handle.index) )
			_receptorBVH.removeParticle(handle.index);
	}
	else if (type == TFAChSlot)
	{
		if( _receptorSlotsBVH.containsParticle(handle.index))
			_receptorSlotsBVH.removeParticle(handle.index);
	}
	if (type == TFSodiumPotassiumPump)
	{
		if( _sodiumPotassiumPumpBVH.containsParticle(handle.index))
			_sodiumPotassiumPumpBVH.removeParticle(handle.index);
	}
	else if (type == TFVoltageGatedIonChannel)
	{
		if( _voltageGatedChannelsBVH.containsParticle(handle.index))
			_voltageGatedChannelsBVH.removeParticle(handle.index);
	}

}

void USynapseSimulation::flexTick(
	float deltaT,
	NvFlexVector<int>& neighbourIndices, 
	NvFlexVector<int>& neighbourCounts, 
	NvFlexVector<int>& apiToInternal, 
	NvFlexVector<int>& internalToAPI, 
	int maxParticles)
{
	_updateReceptorBVH();
	_updateAChSlotBVH();
	_updateVoltageGatedChannelsBVH();
	_updateSodiumPotassiumPumpBVH();
	_updateAChEBVH();
	_updateLeuTBVH();

	_initReceptors(deltaT);

	_tickAChs(deltaT);
	_expireAChSlots(deltaT);
	_bindToAChSlots(deltaT);

	_tickNaCl_and_potassium(deltaT);
	_openCloseReceptors(deltaT);
	_passNaCl(deltaT);
	_pumpSodium(deltaT);
	_pumpPotassium(deltaT);

	_tickAcetycholinesterase(deltaT);

	_tickLeuT(deltaT);

	_openCloseVoltageGatedChannels(deltaT);

	_tickTimers(deltaT);
}

void USynapseSimulation::attach()
{
	graph->addComponentListener<FAcetylcholineReceptor>(this);
}

void USynapseSimulation::detach()
{

}

void USynapseSimulation::_tickTimers(float deltaT)
{
	auto& timers = graph->componentStorage<FSynapseTimerObject>();

	for (FSynapseTimerObject& timer : timers)
		timer.timeLeft -= deltaT;
}

void USynapseSimulation::_bindToAChSlots(float dt)
{
	auto& achs = graph->componentStorage<FAcetylcholine>();
	auto& achSlots = graph->componentStorage<FAChSlot>();

	UTimelineSimulation * timelineSimulation = simulationManager->simulation<UTimelineSimulation>();

	// bind to AChSlots
	for (FAcetylcholine& ach : achs)
	{
		if (!ach.isValid()) continue;
		if( ach.onCooldown() ) continue;

		FVector p = graph->node(ach.nodeHandle()).position;

		// find the neighborhood 
		unrealAABB::AABB query(p, bindingRadius_ACh_Receptor);

		_receptorSlotsBVH.query(query, [&](unsigned int index_AChReceptor) {
			FGraphNodeHandle receptorHandle(index_AChReceptor);

			FAChSlot * slot = achSlots.componentPtrForNode(receptorHandle);

			if (slot->canBind(ach.nodeHandle(), *graph))
			{
				slot->bind(ach.nodeHandle(), *graph);
				slot->timeLeft = expectedBindingTime_ACh;

				// find the root of the slot
				FMLAggregateNO * slotParent = FMLAggregateNO::getAggregate(*graph, slot->nodeHandle());

				FGraphNodeHandle slotParentHandle = slotParent ? slotParent->nodeHandle() : slot->nodeHandle();
				
			
				UEvent_AChR_bindACh& pumpEvent = timelineSimulation->recordEvent<UEvent_AChR_bindACh>(p);
				pumpEvent.triggeringAgent = ach.nodeHandle();
				pumpEvent.otherAgents.Add(slotParentHandle);


				// we're done, we found something to attach to
				return false;
			}

			return true;
		});
	}

	//	_updateBVH<FAcetylcholine>(_achBVH);

	//for (FAChSlot& slot : achSlots)
	//{
	//	if( !slot.isValid() ) continue;
	//	if (slot.timeLeft >= 0.0f) continue;

	//	FVector p = graph->node(slot.nodeHandle()).position;

	//	// find the neighborhood 
	//	unrealAABB::AABB query(p, bindingRadius_ACh_Receptor);


	//	_achBVH.query(query, [&](unsigned int index_ach) {
	//		FGraphNodeHandle achHandle(index_ach);

	//		FAcetylcholine * ach = achs.componentPtrForNode(achHandle);

	//		if (!ach->onCooldown() && slot.canBind(achHandle, *graph))
	//		{
	//			slot.bind(achHandle, *graph);
	//			slot.timeLeft = expectedBindingTime_ACh;

	//			// we're done, we found something to attach to
	//			return false;
	//		}

	//		return true;
	//	});
	//}

}

void USynapseSimulation::_tickAChs(float dt)
{
	auto& achs = graph->componentStorage<FAcetylcholine>();

	for (FAcetylcholine& ach : achs)
	{
		ach.cooldown -= dt;
	}
}

void USynapseSimulation::_expireAChSlots(float dt)
{
	tickSlots<FAChSlot>(dt, [&](FAChSlot& slot) {});

	auto& achSlots = graph->componentStorage<FAChSlot>();
	auto& achs = graph->componentStorage<FAcetylcholine>();

	for (FAChSlot& slot : achSlots)
	{
		slot.timeLeft -= dt;

		if (slot.timeLeft < 0.0f && slot.occupied())
		{
			FAcetylcholine * ach = achs.componentPtrForNode(slot.pigeon());

			slot.unbind(*graph);

			slot.timeLeft = refractoryTime_AChSlot;
		
			if (ach)
				ach->cooldown = generalRefractoryTime;
		}
	}
}

void USynapseSimulation::_initReceptors(float dt)
{
	auto& receptors = graph->componentStorage<FAcetylcholineReceptor>();
	auto& achSlots = graph->componentStorage<FAChSlot>();

	// init
	for (FAcetylcholineReceptor& receptor : receptors)
	{
		if( receptor.didInit() ) continue;

		// if this guy is in an aggregate, use the root of the aggregate to find our partners
		FMLAggregateNO * aggregate = FMLAggregateNO::getAggregate(*graph, receptor.nodeHandle());

		FGraphNodeHandle rootHandle = aggregate ? aggregate->nodeHandle() : receptor.nodeHandle();

		FGraphNode& node = graph->node(rootHandle);

		for (auto ei : node.edges)
		{
			FGraphEdgeHandle eh(ei);

			FGraphEdge& edge = graph->edge(eh);

			FGraphNodeHandle other = edge.other(rootHandle);

			if (achSlots.componentPtrForNode(other))
			{
				// fill up our references
				if (!receptor.achSlot_a)
					receptor.achSlot_a = other;
				else
					receptor.achSlot_b = other;
			}
		}
	}
}

void USynapseSimulation::_openCloseReceptors(float dt)
{
	auto& receptors = graph->componentStorage<FAcetylcholineReceptor>();
	auto& achSlots = graph->componentStorage<FAChSlot>();

	// open/close them
	for (FAcetylcholineReceptor& receptor : receptors)
	{
		if (!receptor.didInit()) continue;

		FAChSlot * slotA = achSlots.componentPtrForNode(receptor.achSlot_a);
		FAChSlot * slotB = achSlots.componentPtrForNode(receptor.achSlot_b);

		receptor.open = slotA->occupied() && slotB->occupied();
	}
}

void USynapseSimulation::_openCloseVoltageGatedChannels(float dt)
{
	auto& gates = graph->componentStorage<FVoltageGatedIonChannel>();
	auto& nacls = graph->componentStorage<FNaCl>();
	

	// goodsell's book is 15.9cm wide
	// scale is 1,000,000 times
	// in unreal, the book is: 84 units across
	// so 1ur = 0.00189m = 1.89e-3

	// and k/r
	// 9.0e-9 / 0.00189e-6
	// so we have: \bar{k} = 9/1.89 = 4.761, our scaling favtor
	const float k_and_scale = 9.0f / 1.89f;

	// reset voltages, ewww nasty O(N^2)
	for (FVoltageGatedIonChannel& gate : gates )
	{
		if( !gate.isValid() ) continue;

		gate.voltage = FVector::ZeroVector;

		const FVector p_gate = graph->node(gate.nodeHandle()).position;
		FVector v = FVector::ZeroVector;

		FVector up = graph->node(gate.nodeHandle()).orientation.RotateVector(FVector::UpVector);

		// Voltage is point charge related to radius from the charge.
		// $V = kQ / r$
		// We can do our scaling at the end and we can just sum these point charges to get the final Voltage

		float actual = 0.0f;

		gate.upCount = 0;
		gate.downCount = 0;

		for (FNaCl& nacl : nacls)
		{
			if( !nacl.isValid() ) continue;

			FVector p_nacl = graph->node(nacl.nodeHandle()).position;

			FVector dir = p_nacl - p_gate;

			float r = dir.Size();

			if( r > sodiumChannenInactivationSampleRadius )
				continue;

			float actualSign = FVector::DotProduct(up, dir) > 0.0f ? 1.0f : -1.0f;

			actual += actualSign * (1.0 / r);

			gate.upCount += actualSign > 0.0f ? 1 : 0;
			gate.downCount += actualSign <= 0.0f ? 1 : 0;

			if (r < 0.11f)
				r = 0.0f;
			else
				r = 1.0 / r;

			// dir / r is a normalized vector
			// dir / r^2 gives us the field
			dir *= std::pow(r, 2.0f);



			v += dir;
		}

		// scale and apply the k constant
		

		gate.voltage = v;
		
		float sign = FVector::DotProduct(up, v) > 0.0f ? 1.0f : -1.0f;

		float voltage = gate.voltage.Size() * sign;

		gate.averageValues.addValue(voltage);

		float averageVoltage = gate.averageValues.average();

		float deltaVoltage = (averageVoltage - gate.previousAverage) / dt;

		gate.previousAverage = averageVoltage;
		gate.deltaVoltage = deltaVoltage;

		// update state

		// state machine
		//if (gate.state == EVoltageGatedIonChannelState::Closed)
		//{
		//	// open
		//	if (deltaVoltage > sodiumChannelActivation_deltaVoltage)
		//	{
		//		gate.timer = sodiumChannenActivationPeriod;
		//		gate.open = true;
		//		gate.state = EVoltageGatedIonChannelState::Open;
		//	}
		//	else
		//		gate.open = false;
		//}
		//else if (gate.state == EVoltageGatedIonChannelState::Open)
		//{
		//	if (gate.timer < 0.0f)
		//	{
		//		gate.open = false;
		//		gate.state = EVoltageGatedIonChannelState::Inactivated;
		//		gate.timer = sodiumChannenInactivationPeriod;
		//	}
		//	else
		//		gate.open = true;
		//}
		//else if (gate.state == EVoltageGatedIonChannelState::Inactivated)
		//{
		//	if (gate.timer < 0.0f)
		//	{
		//		gate.open = false;
		//		gate.state = EVoltageGatedIonChannelState::Closed;
		//	}
		//}

		//if (gate.state == EVoltageGatedIonChannelState::Closed)
		//{
		//	// open
		//	if (averageVoltage > sodiumChannelActivationVoltage)
		//	{
		//		gate.timer = sodiumChannenActivationPeriod;
		//		gate.open = true;
		//		gate.state = EVoltageGatedIonChannelState::Open;
		//	}
		//	else
		//		gate.open = false;
		//}
		//else if (gate.state == EVoltageGatedIonChannelState::Open)
		//{
		//	// inactivate
		//	if (averageVoltage > sodiumChannelInactivationVoltage)
		//	{
		//		gate.timer = sodiumChannenActivationPeriod;
		//		gate.open = false;
		//		gate.state = EVoltageGatedIonChannelState::Inactivated;
		//	}
		//}
		//else if (gate.state == EVoltageGatedIonChannelState::Inactivated)
		//{
		//	if (averageVoltage < sodiumChannelClosedVoltage)
		//	{
		//		gate.open = false;
		//		gate.state = EVoltageGatedIonChannelState::Closed;
		//	}
		//}

		int total = gate.upCount + gate.downCount;
		float ratio = total > 0 ? float(gate.upCount) / float(total) : 1.0f;

		auto state = gate.state;

		if (state == EVoltageGatedIonChannelState::Closed)
		{
			// open
			if (ratio < sodiumChannelActivationRatio )
			{
				gate.open = true;
				gate.state = EVoltageGatedIonChannelState::Open;
			}
			else
				gate.open = false;
		}
		else if (state == EVoltageGatedIonChannelState::Open)
		{
			// inactivate
			if (ratio < sodiumChannelInactivationRatio)
			{
				gate.open = false;
				gate.state = EVoltageGatedIonChannelState::Inactivated;
			}
		}
		else if (state == EVoltageGatedIonChannelState::Inactivated)
		{
			if (ratio > sodiumChannelClosedRatio)
			{
				gate.open = false;
				gate.state = EVoltageGatedIonChannelState::Closed;
			}
		}

		// set materials
		if (FGraphMesh * mesh = graph->componentPtr<FGraphMesh>(gate.nodeHandle()))
		{
			UMaterialInterface * originalMaterial = mesh->material;

			if (gate.state == EVoltageGatedIonChannelState::Open)
			{
				mesh->material = sodiumChannelMaterial_openState;
			}
			else if (gate.state == EVoltageGatedIonChannelState::Inactivated)
			{
				mesh->material = sodiumChannelMaterial_inactivatedState;
			}
			else if (gate.state == EVoltageGatedIonChannelState::Closed)
			{
				mesh->material = sodiumChannelMaterial_closedState;
			}

			if (mesh->material != originalMaterial)
				mesh->markDirty();
		}




		//FString stateAsString = "";
		//{
		//	const UEnum * AlgorithmEnum = FindObject<UEnum>(ANY_PACKAGE, TEXT("EVoltageGatedIonChannelState"));
		//	int32 int32_algorithm = (int32)gate.state;

		//	stateAsString = AlgorithmEnum->GetDisplayNameTextByValue(int32_algorithm).ToString();
		//}


		//UE_LOG(LogTemp, Warning, TEXT("gate is %s"), *stateAsString);

		//UE_LOG(LogTemp, Warning, TEXT("dV: %f V:%f"), deltaVoltage, averageVoltage);

	}


}

void USynapseSimulation::_tickNaCl_and_potassium(float dt)
{
	auto& nacls = graph->componentStorage<FNaCl>();

	for (FNaCl& nacl : nacls)
	{
		nacl.cooldown -= dt;
	}

	auto& ks = graph->componentStorage<FPotassium>();

	for (FPotassium& k : ks)
	{
		k.cooldown -= dt;
	}
}

void USynapseSimulation::_passNaCl(float dt)
{
	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	_tickChannels_withCallback<FNaCl, FAcetylcholineReceptor>(dt, _receptorBVH, generalRefractoryTime, 0.0f, [&](FNaCl& nacl, FAcetylcholineReceptor& gate) {
		FVector position = graph->node(nacl.nodeHandle()).position;

		UEvent_AChR_passSodium& pumpEvent = timeline->recordEvent<UEvent_AChR_passSodium>(position);
		pumpEvent.triggeringAgent = nacl.nodeHandle();
		pumpEvent.otherAgents.Add(gate.nodeHandle());
	});


	_tickChannels_withCallback<FNaCl, FVoltageGatedIonChannel>(dt, _voltageGatedChannelsBVH, generalRefractoryTime, 0.0f, [&](FNaCl& nacl, FVoltageGatedIonChannel& gate) {
		FVector position = graph->node(nacl.nodeHandle()).position;
		
		UEvent_VolgateGatedSodiumChannel_passSodium& pumpEvent = timeline->recordEvent<UEvent_VolgateGatedSodiumChannel_passSodium>(position);
		pumpEvent.triggeringAgent = nacl.nodeHandle();
		pumpEvent.otherAgents.Add(gate.nodeHandle());
	});
}



void USynapseSimulation::_pumpSodium(float dt)
{

	_tickChannels<FNaCl, FSodiumPotassiumPump>(dt, _sodiumPotassiumPumpBVH, generalRefractoryTime, sodiumPumpRefractoryPeriod);
}

void USynapseSimulation::_pumpPotassium(float dt)
{
	_tickChannels<FPotassium, FSodiumPotassiumPump>(dt, _sodiumPotassiumPumpBVH, generalRefractoryTime, sodiumPumpRefractoryPeriod);
}

void USynapseSimulation::_tickLeuT(float dt)
{
	// when we have 2 NaCls and a Choline transport across the membrane
	// we assume there is always a NaCl gradient

	auto& leuts = graph->componentStorage<FLeuTSodiumCoupledTransporter>();
	auto& cholines = graph->componentStorage<FCholine>();
	auto& nacls = graph->componentStorage<FNaCl>();

	// reset
	for (FLeuTSodiumCoupledTransporter& leut : leuts)
	{
		leut.cholineHandle = FGraphNodeHandle::null;
		leut.naclHandle = FGraphNodeHandle::null;
	}

	// do we have choline?
	for (FCholine& choline : cholines)
	{
		if( !choline.isValid() ) continue;

		FGraphNode& cholineNode = graph->node(choline.nodeHandle());

		// find a receptor 
		unrealAABB::AABB query(cholineNode.position, bindingRadius_ACh_Receptor);

		FGraphNodeHandle channelHandle = FGraphNodeHandle::null;

		_LeuTBVH.query(query, [&](unsigned int index_channel) {
			channelHandle = FGraphNodeHandle(index_channel);

			FLeuTSodiumCoupledTransporter * transporter = leuts.componentPtrForNode(channelHandle);

			// keep searching
			if (transporter->cholineHandle) return true;


			FVector dir = cholineNode.position - graph->node(channelHandle).position;

			FVector up = graph->node(channelHandle).orientation.RotateVector(FVector::UpVector);

			// is it on the bottom?
			if (FVector::DotProduct(up, dir) < 0.0f)
			{
				transporter->cholineHandle = choline.nodeHandle();
				return false; // we found it, so quit the query
			}

			// continue the query
			return true;
		});
	}

	// do we have NaCl?
	for (FNaCl& nacl : nacls)
	{
		if( !nacl.isValid() ) continue;

		FGraphNode& naclNode = graph->node(nacl.nodeHandle());

		// find a receptor 
		unrealAABB::AABB query(naclNode.position, bindingRadius_ACh_Receptor);

		FGraphNodeHandle channelHandle = FGraphNodeHandle::null;

		_LeuTBVH.query(query, [&](unsigned int index_channel) {
			channelHandle = FGraphNodeHandle(index_channel);

			FLeuTSodiumCoupledTransporter * transporter = leuts.componentPtrForNode(channelHandle);

			// keep searching
			if (transporter->naclHandle) return true;


			FVector dir = naclNode.position - graph->node(channelHandle).position;

			FVector up = graph->node(channelHandle).orientation.RotateVector(FVector::UpVector);

			// is it on the bottom?
			if (FVector::DotProduct(up, dir) < 0.0f)
			{
				transporter->naclHandle = nacl.nodeHandle();
				return false; // we found it, so quit the query
			}

			// continue the query
			return true;
		});
	}

	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	// transport
	for (FLeuTSodiumCoupledTransporter& leut : leuts)
	{
		if( !leut.isValid() ) continue;
		if( !(leut.cholineHandle && leut.naclHandle)  ) continue;

		FGraphNode& leutNode = graph->node(leut.nodeHandle());

		FGraphNode& naclNode = graph->node(leut.naclHandle);
		FGraphNode& cholineNode = graph->node(leut.cholineHandle);

		FVector up_a = FVector::UpVector + FVector::RightVector;
		FVector up_b = FVector::UpVector - FVector::RightVector;

		FVector up = leutNode.orientation.RotateVector(FVector::UpVector);

		const float length = 1.0f;

		naclNode.position = leutNode.position + leutNode.orientation.RotateVector(up_a * length);
		cholineNode.position = leutNode.position + leutNode.orientation.RotateVector(up_b * length);

		{

			FMLAggregateNO * leutRoot = FMLAggregateNO::getAggregate(*graph, leut.nodeHandle());

			FGraphNodeHandle rootHandle = leutRoot ? leutRoot->nodeHandle() : leut.nodeHandle();

			UEvent_LeuT_passSodium& pumpEvent = timeline->recordEvent<UEvent_LeuT_passSodium>(cholineNode.position);
			pumpEvent.triggeringAgent = naclNode.handle();
			pumpEvent.otherAgents.Add(rootHandle);
			pumpEvent.otherAgents.Add(cholineNode.handle());
		}

	}
}

void USynapseSimulation::_tickAcetycholinesterase(float dt)
{
	auto& AChEs = graph->componentStorage<FAcetylcholinesterase>();
	auto& AChs = graph->componentStorage<FAcetylcholine>();
	auto& meshes = graph->componentStorage<FGraphMesh>();

	UTimelineSimulation * timelineSimulation = simulationManager->simulation<UTimelineSimulation>();

	graph->beginTransaction();

	for (FAcetylcholinesterase& ache : AChEs)
	{
		ache.cooldown -= dt;
	}

	TArray<FGraphNodeHandle> toKill;

	for (FAcetylcholine& ach : AChs)
	{
		if( !ach.isValid() ) continue;

		FGraphNode& achNode = graph->node(ach.nodeHandle());

		// find an AChE 
		unrealAABB::AABB query(achNode.position, bindingRadius_ACh_Receptor);

		FGraphNodeHandle AChEHandle = FGraphNodeHandle::null;

		_AChEBVH.query(query, [&](unsigned int index_AChE) {
			FGraphNodeHandle handle = FGraphNodeHandle(index_AChE);

			if (handle && AChEs.componentPtrForNode(handle)->onCooldown())
				return true;
				
			AChEHandle = handle;

			// done the query
			return false;
		});

		if (AChEHandle)
		{
			// kill the ACh
			FVector position = achNode.position;
			FGraphNodeHandle achHandle = ach.nodeHandle();

			// but keep it's body in the timeline for path vis
			toKill.Add(achHandle);

			//achNode.addComponent<FCholine>(*graph);


			FGraphNodeHandle newCholine = FGraphNodeHandle(graph->addNode(position));
			{
				FGraphNode& cholineNode = graph->node(newCholine);

				cholineNode.addComponent<FFlexParticleObject>(*graph);

				FGraphMesh& mesh = cholineNode.addComponent<FGraphMesh>(*graph);
				mesh.staticMesh = cholineMesh;
				mesh.material = cholineMaterial;
				mesh.markDirty();
				cholineNode.scale = cholineScale;

				cholineNode.addComponent<FCholine>(*graph);
				cholineNode.addComponent<FVelocityGraphObject>(*graph);
				FSingleParticleBrownian& brownian = cholineNode.addComponent<FSingleParticleBrownian>(*graph);


				UEvent_AChE_breakACh& pumpEvent = timelineSimulation->recordEvent<UEvent_AChE_breakACh>(position);
				pumpEvent.triggeringAgent = achHandle;
				pumpEvent.otherAgents.Add(AChEHandle);
				pumpEvent.otherAgents.Add(newCholine);


				FAcetylcholinesterase * ache = AChEs.componentPtrForNode(AChEHandle);

				ache->cooldown = acetycholinesteraseRefractoryPeriod;
			}


		}

	}

	// nuke the components on the old AChs
	static ComponentType TGraphMesh = componentType<FGraphMesh>();

	for (FGraphNodeHandle heDead : toKill)
	{
		FGraphNode& node = graph->node(heDead);

		if (FGraphMesh * mesh = meshes.componentPtrForNode(heDead))
		{
			mesh->staticMesh = nullptr;
			mesh->markDirty();
		}

		// keep the mesh
		auto components = node.components;
		for (auto tComponent : components)
		{
			if( tComponent != TGraphMesh )
				node.removeComponent(*graph, tComponent);
		}
	}

	graph->endTransaction();
}

void USynapseSimulation::_updateReceptorBVH()
{
	_updateBVH<FAcetylcholineReceptor>(_receptorBVH);
}

void USynapseSimulation::_updateAChSlotBVH()
{
	_updateBVH<FAChSlot>(_receptorSlotsBVH);
}

void USynapseSimulation::_updateVoltageGatedChannelsBVH()
{
	_updateBVH<FVoltageGatedIonChannel>(_voltageGatedChannelsBVH);
}

void USynapseSimulation::_updateSodiumPotassiumPumpBVH()
{
	_updateBVH<FSodiumPotassiumPump>(_sodiumPotassiumPumpBVH);
}

void USynapseSimulation::_updateAChEBVH()
{
	_updateBVH<FAcetylcholinesterase>(_AChEBVH);
}

void USynapseSimulation::_updateLeuTBVH()
{
	_updateBVH<FLeuTSodiumCoupledTransporter>(_LeuTBVH);
}

bool FAChSlot::canBind(FGraphNodeHandle target, FGraph& graph)
{
	bool okay = timeLeft < 0.0f && !occupied();

	return okay && !graph.node(target).hasComponent<FBoundToParticleSlot>();
}

bool FAcetylcholineReceptor::didInit()
{
	return achSlot_a && achSlot_b;
}

float FMovingAverage::average()
{
	float sum = 0.0f;

	for (auto& f : values)
		sum += f;

	if (values.Num() == 0)
		return 0.0f;
	else
		return sum / float(values.Num());
}

void FMovingAverage::addValue(float value)
{
	values.SetNum(maxValues);

	values[nextIndex] = value;

	nextIndex = (nextIndex + 1) % maxValues;
}

void FMovingAverage::clear()
{
	for (auto& f : values)
		f = 0.0f;
}
