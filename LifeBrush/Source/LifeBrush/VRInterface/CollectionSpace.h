// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "CoreMinimal.h"
#include "Components/SceneComponent.h"
#include "CollectionSpace.generated.h"

UCLASS(ClassGroup = (Custom), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UInteractionSpace : public USceneComponent
{
	GENERATED_BODY()

public:
	UInteractionSpace();
 
};








UINTERFACE(Blueprintable)
class UCollectionSpaceDataSource : public UInterface
{
	GENERATED_BODY()
};

class ICollectionSpaceDataSource
{
	GENERATED_BODY()

public:
	UFUNCTION(BlueprintCallable, BlueprintNativeEvent, Category = "VRUI")
	int32 count();

	UFUNCTION(BlueprintCallable, BlueprintNativeEvent, Category = "VRUI")
	UPrimitiveComponent * primitiveCellAt(int32 index);
};

UINTERFACE(Blueprintable)
class UCollectionSpaceDelegate : public UInterface
{
	GENERATED_BODY()
};

class ICollectionSpaceDelegate
{
	GENERATED_BODY()

public:
	UFUNCTION(BlueprintCallable, BlueprintNativeEvent, Category = "VRUI")
	void didGrab(int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds);
};

UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UCollectionSpace : public UInteractionSpace
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UCollectionSpace();

protected:
	// Called when the game starts
	virtual void BeginPlay() override;

public:	
	// Called every frame
	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;

	virtual void insertItemsAt(TArray<uint32> indices);

	virtual void reloadData();

	FVector nearest(UPrimitiveComponent * interactionPoint);

	// Interaction events
	virtual void begin_oneHand(UPrimitiveComponent * interactionPoint);
	virtual void update_oneHand(float dt, UPrimitiveComponent * interactionPoint, FTransform lastTransform);
	virtual void end_oneHand(UPrimitiveComponent * interactionPoint);

	virtual void grab(UPrimitiveComponent * interactionPoint);

	virtual void query(UPrimitiveComponent * interactionPoint);

	virtual TSet<int32> selection();

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	float cellExtents = 10.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	float cellSpacing = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	float damping = 0.1f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	float width = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	int32 rows = 1;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	UStaticMesh * backgroundMesh;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI")
	UMaterialInterface * backgroundMaterial;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI", meta = (MustImplement = "CollectionSpaceDataSource"))
	UObject * dataSource;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI", meta = (MustImplement = "CollectionSpaceDelegate"))
	UObject * delegate;

protected:
	void _layout();
	void _layoutCells();
	void _layoutBackground();

	float _totalWidth();
	float _stepSize();

	float _scaleForCell(UPrimitiveComponent& component);
	FVector _offsetForCell(const FTransform& toCell, UPrimitiveComponent& component);

	void _update_pan(float dt, UPrimitiveComponent * interactionPoint, FTransform lastTransform);

protected:
	UPROPERTY()
	UStaticMeshComponent * _boundsMesh;

protected:
	enum class InteractionMode
	{
		Pan,
		Grab, 
		None
	};

	InteractionMode _interactionMode = InteractionMode::None;

	TArray<UPrimitiveComponent*> _cells;
	TArray<FBox> _localBounds;

	FBox _bounds;

	FVector _velocity = FVector::ZeroVector;

	TSet<int32> _selection;

	static const int Forward = 1;
};
