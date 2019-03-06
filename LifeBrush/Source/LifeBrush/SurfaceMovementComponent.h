// Copyright 2015, Timothy Davison. All rights reserved.

#pragma once

#include "Components/ActorComponent.h"
#include "RegionGrowingComponent.h"

#include <random>

#include "SurfaceMovementComponent.generated.h"


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API USurfaceMovementComponent : public UActorComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	USurfaceMovementComponent();

	// Called when the game starts
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction ) override;

    // An actor that has a region growing component
    UPROPERTY(EditAnywhere) AActor * regionGrowingActor;
    
    // max acceleration magnitude of the random walk
    UPROPERTY(EditAnywhere) float maxAcceleration = 4.0f;
    
    // max duration of a walk cycle
    UPROPERTY(EditAnywhere) float maxWalkTime = 2.0f;
    
    float t = 0;
    
    float r1 = 0.0f;
    float r2 = 0.0f;
    
    FQuat hackStartRotation;
    
    
private:
    URegionGrowingComponent* _regionGrowingComponent()
    {
        if( regionGrowingActor == nullptr )
            return nullptr;
        
        URegionGrowingComponent * component = regionGrowingActor->FindComponentByClass<URegionGrowingComponent>();
        
        return component;
    }
	
};
