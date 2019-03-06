// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "LifeBrush.h"

#include "NeighbourhoodParameters.generated.h"

USTRUCT( Blueprintable )
struct FNeighbourhoodParameters
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float radius;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 kNearest;
};