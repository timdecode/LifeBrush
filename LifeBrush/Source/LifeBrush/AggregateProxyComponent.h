// Copyright 2019, Timothy Davison. All rights reserved.

#pragma once

#include "TimStructBox.h"
#include "ShipEditorSimulation/Graph.h"

#include "AggregateProxyComponent.generated.h"

UENUM(BlueprintType)
enum class EAggregateProxyParticleOrMember : uint8
{
	Particle UMETA(DisplayName = "Particle"),
	Member UMETA(DisplayName = "Member"),
};

UCLASS(ClassGroup = (LifeBrush), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UAggregateProxyComponent : public UStaticMeshComponent
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "Agent Library" )
	TArray<FTimStructBox> graphObjects;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	EAggregateProxyParticleOrMember particleOrMember = EAggregateProxyParticleOrMember::Particle;
};