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
