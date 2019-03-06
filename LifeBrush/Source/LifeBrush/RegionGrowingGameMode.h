// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "GameFramework/GameMode.h"
#include "RegionGrowingGameMode.generated.h"

/**
 * 
 */
UCLASS()
class LIFEBRUSH_API ARegionGrowingGameMode : public AGameMode
{
	GENERATED_BODY()
	
    UPROPERTY(EditAnywhere) int32 elementID = 0;
	
	
};
