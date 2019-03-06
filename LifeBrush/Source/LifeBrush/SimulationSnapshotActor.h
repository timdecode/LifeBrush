// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/Actor.h"
#include "SimulationSnapshotActor.generated.h"

UCLASS()
class LIFEBRUSH_API ASimulationSnapshotActor : public AActor
{
	GENERATED_BODY()
	
protected:
	// Hack to quickly show/hide child ISMCs in the editor
	UPROPERTY( EditAnywhere, BlueprintReadWrite )
	bool visibile_allChildren = true;

	UPROPERTY( VisibleAnywhere, BlueprintReadOnly ) 
	int32 statElementCount = 0;

public:	
	// Sets default values for this actor's properties
	ASimulationSnapshotActor();

protected:
	// Called when the game starts or when spawned
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	virtual bool ShouldTickIfViewportsOnly() const { return true; }

#if WITH_EDITOR
	void PostEditChangeProperty( struct FPropertyChangedEvent& e ) override;
#endif   
	
protected:
	void _toggleVisibility();
};
