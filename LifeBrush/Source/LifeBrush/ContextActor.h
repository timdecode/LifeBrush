// Copyright 2015, Timothy Davison. All rights reserved.

#pragma once

#include "GameFramework/Actor.h"
#include <memory>
#include "Algorithm/ctpl_stl.h"

#include "ContextActor.generated.h"

UCLASS()
class LIFEBRUSH_API AContextActor : public AActor
{
	GENERATED_BODY()
	
public:
    UPROPERTY(EditAnywhere) float synthesizeDistance;
    
    UPROPERTY(EditAnywhere) float innerQuellingDistance = 0.0f;
	UPROPERTY( EditAnywhere ) float outerQuellingDistance = 0.0f;

	UPROPERTY( EditAnywhere ) float hideSurfaceDistance = 0.0f;
	UPROPERTY( EditAnywhere ) float resetDistance = 0.0f;
    
    UPROPERTY(EditAnywhere) FVector generationLimitsMin = {-1000.0f, -1000.0f, -1000.0f};
    UPROPERTY(EditAnywhere) FVector generationLimitsMax = { 1000.0f,  1000.0f,  1000.0f};


    
public:	
	// Sets default values for this actor's properties
	AContextActor();

	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void Tick( float DeltaSeconds ) override;

    std::shared_ptr<ctpl::thread_pool> threadPool()
    {
        if( _threadPool == nullptr )
        {
            unsigned int n = std::thread::hardware_concurrency();

            _threadPool = std::make_shared<ctpl::thread_pool>();
            _threadPool->resize(n);
        }
        
        return _threadPool;
    }
    
private:
    std::shared_ptr<ctpl::thread_pool> _threadPool;
};
