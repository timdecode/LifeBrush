// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"

#include "MeshFilamentSimulation.generated.h"

class UColoredLineFactory;

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFilamentHead : public FGraphEdgeObject
{
	GENERATED_BODY()
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFilamentConnection : public FGraphEdgeObject
{
	GENERATED_BODY()

public:

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	float radius = 0.0f;

	UPROPERTY()
	uint32 group = 0;

	UPROPERTY()
	uint32 segmentID = 0;
};

UCLASS( BlueprintType )
class LIFEBRUSH_API UMeshFilamentSimulation : public UObjectSimulation, public EdgeObjectListener
{
	GENERATED_BODY()

public:
	virtual void attach() override;
	virtual void detach() override;

	virtual void begin() override;

	virtual void tick(float deltaT) override;

	uint64 nextGroup() { auto next = _nextGroup; _nextGroup++; return next; }


	virtual void edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type) override;

protected:
	void _updateFilaments();

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	UMaterialInterface * defaultLineMaterial = nullptr;

protected:
	UPROPERTY()
	UColoredLineFactory * _lineFactory;

	// It's highly unlikely that we'll ever run out of group-ids with a 32-bit integer.
	UPROPERTY() 
	uint32 _nextGroup = 0;

	size_t lastCount = 0;

	bool _needsSort = true;
};