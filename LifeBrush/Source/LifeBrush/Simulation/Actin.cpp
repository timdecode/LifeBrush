// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "FlexElements.h"

#include "MolecularLego/MolecularLego.h"
#include "Aggregates.h"

#include "Actin.h"
#include <bitset>

void UMLActinSimulation::attach()
{

}

void UMLActinSimulation::detach()
{

}

void UMLActinSimulation::tick(float deltaT)
{
	if(true) _cacheBonds();

	auto& actins = graph->componentStorage<FMLActin>();
	auto& particles = graph->componentStorage<FMLOrientedParticle>();

	//// increase/decrease strength based on bound angle distances
	//for (FMLActin& actin : actins)
	//{
	//	float cost = 0.0f;

	//	for (FGraphNodeHandle a : actin._bonds)
	//	{
	//		FMLOrientedParticle * particle_a = particles.componentPtrForNode(a);

	//		if (!particle_a) {
	//			continue;
	//		}

	//		FGraphNodeHandle b = particle_a->_bondedPartner;

	//		if (!b) {
	//			// more cost if we are unbound
	//			cost += 0.1f;
	//			continue;
	//		}

	//		const FQuat q_a = graph->node(a).orientation;
	//		const FQuat q_b = UMLParticleSimulation::mirror(graph->node(b).orientation);

	//		cost += q_a.AngularDistance(q_b) * 0.5f;
	//	}

	//	// scale the interaction according to this angle
	//	for (FGraphNodeHandle a : actin._bonds)
	//	{
	//		FMLOrientedParticle * particle_a = particles.componentPtrForNode(a);

	//		particle_a->radiusMultiplier = (1 / (1 + cost));
	//	}
	//}

	_tickRules2();

	//_tickRules();
}

void UMLActinSimulation::_tickRules2()
{
	auto& actins = graph->componentStorage<FMLActin>();
	auto& particles = graph->componentStorage<FMLOrientedParticle>();

	UMLParticleSimulation * mlParticleSimulation = simulationManager->simulation<UMLParticleSimulation>();

	//for (FGraphNodeHandle handle : mlParticleSimulation->recentlyBonded)
	//{
	for( FMLActin& actin : actins )
	{


		std::bitset<FMLActin::_bondsSize> isBonded;
		for (int i = 0; i < FMLActin::_bondsSize; ++i)
		{
			isBonded[i] = particles.componentPtrForNode(actin._bonds[i])->numActiveBonds > 0;
		}

		auto bonded = isBonded.to_ulong();

		auto setActive = [&](int index, bool active) {
			auto particle = particles.componentPtrForNode(actin._bonds[index]);
			particle->active = active;
			particle->isSlave = false;
		};

		auto setActiveSlave = [&](int index, bool active) {
			auto particle = particles.componentPtrForNode(actin._bonds[index]);
			particle->active = active;
			particle->isSlave = true;
		};

		switch (bonded)
		{
		case 0b0000:
			setActiveSlave(0, true);
			setActiveSlave(1, true);
			setActiveSlave(2, false);
			setActiveSlave(3, false);
			break;

		case 0b0001:
			setActive(0, false);
			setActive(1, false);
			setActive(2, false);
			setActive(3, true);
			break;

		case 0b0010:
			setActive(0, false);
			setActive(1, false);
			setActive(2, true);
			setActive(3, false);
			break;

		case 0b0110:
			setActive(0, false);
			setActive(1, false);
			setActive(2, false);
			setActive(3, true);
			break;

		case 0b1001:
			setActive(0, false);
			setActive(1, false);
			setActive(2, true);
			setActive(3, false);
			break;

		case 0b0100:
			setActiveSlave(0, false);
			setActiveSlave(1, false);
			setActiveSlave(2, false);
			setActiveSlave(3, true);
			break;

		case 0b1000:
			setActiveSlave(0, false);
			setActiveSlave(1, false);
			setActiveSlave(2, true);
			setActiveSlave(3, false);
			break;
		
		case 0b1100:
			setActive(0, true);
			setActive(1, true);
			setActive(2, false);
			setActive(3, false);
			break;

		case 0b1110:
			setActive(0, true);
			setActive(1, false);
			setActive(2, false);
			setActive(3, false);
			break;

		case 0b1101:
			setActive(0, false);
			setActive(1, true);
			setActive(2, false);
			setActive(3, false);
			break;

		default:
			break;
		}
	}
}


void UMLActinSimulation::_tickRules()
{
	auto& actins = graph->componentStorage<FMLActin>();
	auto& particles = graph->componentStorage<FMLOrientedParticle>();

	UMLParticleSimulation * mlParticleSimulation = simulationManager->simulation<UMLParticleSimulation>();

	for (FGraphNodeHandle handle : mlParticleSimulation->recentlyBonded)
	{
		auto aggregate = FMLAggregateNO::getAggregate(*graph, handle);
		if (!aggregate) continue;

		FMLActin * actin = actins.componentPtrForNode(aggregate->nodeHandle());

		if (!actin) continue;

		const int type0Hack = 0;
		const int type1Hack = 1;
		const int type2Hack = 2;
		const int type3Hack = 3;

		const float lowMultiplier = 0.4f;
		const float highMultilier = 1.0f;

		std::bitset<FMLActin::_bondsSize> bonded;
		for (int i = 0; i < FMLActin::_bondsSize; ++i)
		{
			bonded[i] = particles.componentPtrForNode(actin->_bonds[i])->numActiveBonds > 0;
		}

		if (!bonded[0] && bonded[1] && bonded[2] && bonded[3])
		{
			if (auto particle = particles.componentPtrForNode(actin->_bonds[0]))
			{
				particle->radiusMultiplier = highMultilier;
			}
		}
		else if (bonded[0] && !bonded[1] && !bonded[2] && bonded[3])
		{
			if (auto particle = particles.componentPtrForNode(actin->_bonds[2]))
			{
				particle->radiusMultiplier = highMultilier;
			}

			if (auto particle = particles.componentPtrForNode(actin->_bonds[1]))
			{
				particle->radiusMultiplier = lowMultiplier;
			}
		}
		else
		{
			for (FGraphNodeHandle h : actin->_bonds)
			{
				particles.componentPtrForNode(h)->interactionMultiplier = lowMultiplier;
				particles.componentPtrForNode(h)->radiusMultiplier = lowMultiplier;
			}
		}
	}
}

void UMLActinSimulation::_cacheBonds()
{
	auto& actinStorage = graph->componentStorage<FMLActin>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();

	// cache and update bonds bonds
	for (FMLActin& actin : actinStorage)
	{
		if (actin._bonds[0] != FGraphNodeHandle::null)
			continue;

		auto rigidBodyHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, actin.nodeHandle());

		int i = 0;
		TArray<FMLOrientedParticle*, TInlineAllocator<16>> bonds;

		rigidBodyHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle handle, FFlexRigidBodyConnection& connection) {
			FMLOrientedParticle * particle = orientedParticles.componentPtrForNode(handle);

			if (!particle) return;

			bonds.Add(particle);
		});

		bonds.Sort([&](FMLOrientedParticle& a, FMLOrientedParticle& b) {
			return a.speciesID < b.speciesID;
		});

		for (int i = 0; i < actin._bonds.size() && i < bonds.Num(); ++i)
		{
			actin._bonds[i] = bonds[i]->nodeHandle();
		}
	}

	_bondsDirty = false;
}
