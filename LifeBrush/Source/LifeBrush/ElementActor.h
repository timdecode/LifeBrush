// Copyright 2016, Timothy Davison. All rights reserved.

#pragma once

#include "GameFramework/Actor.h"

#include "TimStructBox.h"
#include "Algorithm/Element.h"

#include "NeighbourhoodParameters.h"

#include "AggregateProxyComponent.h"

#include "ElementActor.generated.h"

class UGraphSimulationManager;

UENUM(Blueprintable)
enum class ESpaceMode : uint8
{
	Volume,
	Surface,
};

UCLASS()
class LIFEBRUSH_API AElementActor : public AStaticMeshActor
{
	GENERATED_BODY()
	
public:


	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "Agent Library" )
	TArray<FTimStructBox> graphObjects;

	// Optional group name for the element. This lets us quickly select an entire group of elements, if that option is enabled.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "RG")
	FString groupName; 

    // overrides the inferred radius of an element
    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
    float overrideRadius = 1.0f;

    // whether the override radius is affected by the actor's scale
    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
    bool scaleOverrideRadius = true;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	float gradient = 0.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	bool generative = true;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	FNeighbourhoodParameters generationParameters;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	FNeighbourhoodParameters optimizationParameters;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	float minAssignmentDistance = 1.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	float freespaceRadius = 1.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "RG" )
	float generationInnerRadius = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "RG")
	ESpaceMode spaceModeHint = ESpaceMode::Volume;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "RG")
	EAggregateProxyParticleOrMember particleOrMember = EAggregateProxyParticleOrMember::Particle;

public:	
	// Sets default values for this actor's properties
	AElementActor();

	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void Tick( float DeltaSeconds ) override;

	void writeToElement(ElementTuple& element, FGraph& graph);
	void readFromElement( ElementTuple& element, FGraph& graph);
	void readFromActor( AElementActor* elementActor );

	void showSelectionOutline();
	void hideSelectionOutline();

protected:
	void _loadAggregate(FGraphNodeHandle elementNode, FGraph& graph);
};
