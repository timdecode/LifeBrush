// Copyright 2019, Timothy Davison. All rights reserved.

#pragma once

#include "ElementActor.h"


#include "CoreMinimal.h"
#include "Simulation/FlexElements.h"
#include "Components/ActorComponent.h"
#include "MoleculeCreator.generated.h"


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UMoleculeCreator : public UActorComponent
{
	GENERATED_BODY()

public:	
	

	//UPROPERTY(EditAnywhere)
		//bool spawnMol;

	UPROPERTY(EditAnywhere)
		UStaticMesh* mesh;

	bool isParticle = true;

	void SpawnEActor();
	void SetGraphObjects(AElementActor*);

	UMoleculeCreator();

protected:
	// Called when the game starts
	virtual void BeginPlay() override;

	

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;

		
	
};
