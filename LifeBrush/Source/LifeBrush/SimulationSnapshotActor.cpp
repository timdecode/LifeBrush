// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"
#include "SimulationSnapshotActor.h"


// Sets default values
ASimulationSnapshotActor::ASimulationSnapshotActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;


}

// Called when the game starts or when spawned
void ASimulationSnapshotActor::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void ASimulationSnapshotActor::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

	TArray<USceneComponent*> children;
	GetRootComponent()->GetChildrenComponents( false, children );

	statElementCount = 0;

	for(USceneComponent * child : children)
	{
		UInstancedStaticMeshComponent * instance = Cast<UInstancedStaticMeshComponent>( child );

		if(!instance)
			continue;

		statElementCount += instance->GetInstanceCount();
	}
}

#if WITH_EDITOR
void ASimulationSnapshotActor::PostEditChangeProperty( struct FPropertyChangedEvent& e )
{
	FName propertyName = (e.Property != NULL) ? e.Property->GetFName() : NAME_None;

	if(propertyName == GET_MEMBER_NAME_CHECKED( ASimulationSnapshotActor, visibile_allChildren ))
	{
		_toggleVisibility();
	}

	Super::PostEditChangeProperty( e );
}

void ASimulationSnapshotActor::_toggleVisibility()
{
	if(!RootComponent)
		return;

	RootComponent->SetVisibility( visibile_allChildren, true );
}

#endif