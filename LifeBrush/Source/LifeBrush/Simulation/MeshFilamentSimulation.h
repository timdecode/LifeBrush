// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"

#include "MeshFilamentSimulation.generated.h"

class UColoredLineFactory;

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFilamentConnection : public FGraphEdgeObject
{
	GENERATED_BODY()

public:
	// extends the filament in aExtension units in the a direction
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float aExtension = 0.0f;

	// extends the filament in bExtension units in the b direction
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float bExtension = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;

	UPROPERTY()
	uint32 group = 0;

	UPROPERTY()
	int32 segmentID = 0;

	UPROPERTY()
	bool visible = true;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFrozenFilament : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FVector position = FVector::ZeroVector;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FQuat orientation = FQuat::Identity;
};

UCLASS( BlueprintType )
class LIFEBRUSH_API UMeshFilamentSimulation : public UObjectSimulation, public EdgeObjectListener
{
	GENERATED_BODY()

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void begin() override;

	virtual void tick(float deltaT) override;
	virtual void tick_paused(float deltaT) override;

	
	uint32 nextGroup(UMaterialInterface * material = nullptr);


	virtual void edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type) override;

	virtual bool canRunInEditor() override { return true; }

	virtual void snapshotToActor(AActor * actor) override;


public:
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void freezeFilaments();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void unfreezeFilaments();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	bool areFilamentsFrozen();

	void setDesaturated(bool desaturated);

protected:
	void _updateFilaments();

	void _updateFrozenComponents();

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * defaultLineMaterial = nullptr;



protected:
	UPROPERTY()
	UColoredLineFactory* _lineFactory;

	UPROPERTY()
	TMap<int32, int32> groupToSection;

	UPROPERTY()
	TMap<UMaterialInterface*, int32> materialToSection;

	UPROPERTY()
	bool isDesaturated = false;

	// It's highly unlikely that we'll ever run out of group-ids with a 32-bit integer.
	UPROPERTY() 
	uint32 _nextGroup = 0;

	size_t lastCount = 0;



	bool _dirty = true;

	bool _areFilamentsFrozen = false;


protected:
	int32 addSection(UMaterialInterface * material);


	void _updateDesaturated();
};