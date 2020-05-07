// Copyright 2019, Timothy Davison. All rights reserved.

#include "LifeBrush.h"
#include "MoleculeCreator.h"


// Sets default values for this component's properties
UMoleculeCreator::UMoleculeCreator()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	// ...
}

void UMoleculeCreator::SpawnEActor() {

	

	FActorSpawnParameters params;
	params.Name = "Spawned Molecule";

	FVector location(105.f, 31.f, 97.f);
	FVector scale(0.02, 0.02, 0.02);
	FQuat rotation(0., 0., 0., 0.);
	FTransform trans(rotation, location, scale);
	AElementActor* actor = GetWorld()->SpawnActor<AElementActor>(params);

	actor->SetActorRelativeTransform(trans);
	actor->SetOwner(this->GetOwner());
	SetGraphObjects(actor);

	UStaticMeshComponent* comp = actor->GetStaticMeshComponent();

	comp->SetStaticMesh(mesh);

	
	

}

void UMoleculeCreator::SetGraphObjects(AElementActor* actor) {

	FRandomWalkGraphObject rndWalkGraphObj;
	FFlexParticleObject flxParticleObj;
	

	FTimStructBox rndStructBox;
	rndStructBox.scriptStruct = rndWalkGraphObj.StaticStruct();
	rndStructBox.initStruct();

	FTimStructBox flxPartStrucBox;
	flxPartStrucBox.scriptStruct = flxParticleObj.StaticStruct();
	flxPartStrucBox.initStruct();
	
	

	actor->graphObjects.Push(rndStructBox);
	actor->graphObjects.Push(flxPartStrucBox);


	

}

// Called when the game starts
void UMoleculeCreator::BeginPlay()
{
	Super::BeginPlay();
	
	SpawnEActor();
	// ...
	
}


// Called every frame
void UMoleculeCreator::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	// ...
}

