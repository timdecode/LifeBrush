//
//  AlgorithmEnums.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2017-11-21.
//  Copyright © 2017 Timothy Davison, Inc. All rights reserved.
//

#pragma once

#include "AlgorithmEnums.generated.h"

UENUM(BlueprintType)
enum class EGenerationMode : uint8
{
    SpaceFilling UMETA(DisplayName="Space Filling"),
    SurfaceProjection UMETA(DisplayName="Surface Projection"),
	SurfaceWalking UMETA(DisplayName="Surface Walking"),
    SurfacePainting UMETA(DisplayName="Surface Painting"),
	SpacePainting UMETA(DisplayName="Space Painting"),
};

UENUM( BlueprintType )
enum class EEnergyCalculationMode : uint8 
{
	Coherence UMETA( DisplayName = "Coherence" ),
	BruteForce UMETA( DisplayName = "Brute Force" )
};
