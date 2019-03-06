// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "VRTool.h"
#include "RegionGrowingToolInit.h"

#include "VRInterface/CollectionSpace.h"
#include "GrabTool.h"

#include "MeshCollectionTool.generated.h"


class UGrabTool;

UCLASS(DefaultToInstanced)
class LIFEBRUSH_API AMeshCollectionSpaceActor : public AActor, public ICollectionSpaceDataSource, public ICollectionSpaceDelegate
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "LifeBrush")
	UCollectionSpace * collectionSpace;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TArray<UStaticMesh*> staticMeshes;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	class UMeshCollectionTool * meshTool = nullptr;

public:
	AMeshCollectionSpaceActor();

	
	virtual bool ShouldTickIfViewportsOnly() const override;
	virtual void Tick(float DeltaSeconds);

public:
	// ICollectionSpaceDataSource
	virtual int32 count_Implementation() override;

	virtual UPrimitiveComponent * primitiveCellAt_Implementation(int32 index) override;

	// ICollectionSpaceDelegate
	virtual void didGrab_Implementation(int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds) override;


	// Others
	virtual void BeginPlay() override;

protected:
	bool _didBegingPlayInEditor = false;
};

// For creating elements from a mesh, or reassigning the mesh of an existing element.
UCLASS(Blueprintable)
class UMeshCollectionTool : public UTool, public IGrabDelegate
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TSubclassOf<class AMeshCollectionSpaceActor> meshCollectionClass;

public:
	void init(FRGC_UToolInitProperties& initProperties, UCameraComponent * camera)
	{
		UTool::init(initProperties);

		_regionGrowingComponent = initProperties.regionGrowingComponent;
		_camera = camera;
	}

	virtual void focused() override;
	virtual void loseFocus() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override;
	virtual void twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override;

	virtual void tick(float dt);

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform lastTransformA,
		FTransform lastTransformB
	) override;

	virtual void faceDown_released();
	virtual void faceUp_released(USceneComponent * interactionPoint /* = nullptr */);


	UFUNCTION()
	void OnRightSelectionPointOverlapBegin(class UPrimitiveComponent* ourComponent, class AActor* OtherActor, class UPrimitiveComponent* OtherComp, int32 OtherBodyIndex, bool bFromSweep, const FHitResult & SweepResult);

	UFUNCTION()
	void OnRightSelectionPointOverlapEnd(class UPrimitiveComponent* ourComponent, class AActor* OtherActor, class UPrimitiveComponent* OtherComp, int32 OtherBodyIndex);

	void didGrab(UCollectionSpace * collectionSpace, int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds);

	// IGrabDelegate
	virtual void didPlace_Implementation(UGrabTool * grabTool, AActor * draggingActor) override;
	virtual void didCancel_Implementation(UGrabTool * grabTool, AActor * draggingActor) override;

protected:
	void _spawnCollectionSpaceActor();

protected:
	class URegionGrowingComponent * _regionGrowingComponent;
	class UCameraComponent * _camera;
	class AMeshCollectionSpaceActor * _collectionSpaceActor;
	class UCollectionSpace * _collectionSpace;

	UPROPERTY()
	AActor * _draggingActor = nullptr;

	UPROPERTY()
	UStaticMeshComponent * _draggingCell = nullptr;

	UPROPERTY()
	UGrabTool * _grabTool = nullptr;
};