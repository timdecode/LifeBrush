// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Brownian.h"

void USingleParticleBrownianSimulation::attach()
{
	rand.GenerateNewSeed();
}

void USingleParticleBrownianSimulation::detach()
{

}

void USingleParticleBrownianSimulation::tick(float deltaT)
{
	auto& browns = graph->componentStorage<FSingleParticleBrownian>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	// make sure we have velocities
	for (FSingleParticleBrownian& brown : browns)
	{
		if( !brown.isValid() ) continue;

		if (!velocities.componentPtrForNode(brown.nodeHandle()))
		{
			FGraphNode& node = graph->node(brown.nodeHandle());
			node.addComponent<FVelocityGraphObject>(*graph);
		}
	}

	for (FSingleParticleBrownian& brown : browns)
	{
		if( !brown.isValid() ) continue;
		
		brown.time -= deltaT;

		if (brown.time > 0.0f)
			continue;

		brown.time = rand.FRandRange(minTime, maxTime);

		float scale = FMath::Clamp(1.0f - brown.dampening, 0.0f, 1.0f);

		FVector dv = scale * rand.GetUnitVector() * rand.FRandRange(minSpeed, maxSpeed);

		if (auto velocity = velocities.componentPtrForNode(brown.nodeHandle()))
		{
			velocity->linearVelocity = (velocity->linearVelocity + dv).GetClampedToMaxSize(15.0f);
		}
	}
}

void UGlobalParticleBrownianSimulation::attach()
{
	rand.GenerateNewSeed();
}

void UGlobalParticleBrownianSimulation::detach()
{

}

void UGlobalParticleBrownianSimulation::tick(float deltaT)
{
	if (!enabled) return;

	auto& particles = graph->componentStorage<FFlexParticleObject>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	for (FFlexParticleObject& particle : particles)
	{
		if (!particle.isValid()) continue;

		if (!velocities.componentPtrForNode(particle.nodeHandle()))
		{
			FGraphNode& node = graph->node(particle.nodeHandle());

			node.addComponent<FVelocityGraphObject>(*graph);
		}
	}

	// find _times Num
	int32 timeSize = 0;

	for (FFlexParticleObject& particle : particles)
	{
		if (!particle.isValid()) continue;

		int32 timeIndex = particle.nodeHandle().index;

		if (timeIndex > timeSize)
			timeSize = timeIndex;
	}

	_times.SetNum(timeSize + 1);

	// make sure we have velocities
	for (FFlexParticleObject& particle : particles)
	{
		if (!particle.isValid()) continue;

		int32 timeIndex = particle.nodeHandle().index;

		float& timeLeft = _times[timeIndex];

		timeLeft -= deltaT;

		if (timeLeft > 0.0f)
			continue;

		timeLeft = rand.FRandRange(minTime, maxTime);


		FVelocityGraphObject * velocity = velocities.componentPtrForNode(particle.nodeHandle());

		FVector dv = rand.GetUnitVector() * rand.FRandRange(minSpeed, maxSpeed);

		velocity->linearVelocity += dv;

	}
}
