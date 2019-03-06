// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "Flex/include/NvFlex.h"
#include "Flex/include/NvFlexExt.h"
#include "Flex/include/NvFlexDevice.h"

#include "FlexGraphSimulation_interface.generated.h"

UINTERFACE(Blueprintable)
class UFlexGraphSimulation : public UInterface
{
	GENERATED_BODY()
};

// This is a class for interfacing with a FFlexSimulation.
// \func preTick is called before UObjectSimulation::tick
// \func postTick is called after UObjectSimulation::tick and FlexGraphSimulation::flexTick.
// All ticks are applied in group, successively, to each simulation in a UGraphSimulationManager
// within a FFlexSimulation.
class IFlexGraphSimulation
{
	GENERATED_BODY()

public:
	// No modifications to the graph should be made during a preTick.
	virtual void preTick(float deltaT) {};

	virtual void flexTick(
		float deltaT,
		NvFlexVector<int>& neighbourIndices,
		NvFlexVector<int>& neighbourCounts,
		NvFlexVector<int>& apiToInternal,
		NvFlexVector<int>& internalToAPI,
		int maxParticles
	) {};

	// No modifications to the graph should be made during a postTick
	virtual void postTick(float deltaT) {};
};



