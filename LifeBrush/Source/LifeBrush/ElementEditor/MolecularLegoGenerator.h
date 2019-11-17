// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ElementGenerator.h"

#include "ElementGenerator.h"
#include "Algorithm/Algorithm.h"
#include "Algorithm/PointCloud.hpp"
#include "aabbcc/unrealAABB.h"
#include "LinearPath.h"

#include "MolecularLegoGenerator.generated.h"

UENUM(BlueprintType)
enum class EMolecularLegotGeneratorMode : uint8
{
	AlongPath UMETA(DisplayName = "Along path"),
	Grow UMETA(DisplayName = "Grow path"),
	PlaceOne UMETA(DisplayName = "Place one")
};

UCLASS(ClassGroup = (Custom), DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UMolecularLegoGenerator : public UElementGenerator
{
	GENERATED_BODY()

public:
	SynthesisContext * _context;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

	// The distance between brush point segments.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float segmentLength = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float elementRadius = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	EMolecularLegotGeneratorMode mode;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxBrushRadius = 5.0f;

protected:
	UGraphSimulationManager * _exampleSimulationManager;

	typedef unrealAABB::Tree BVH_t;
	BVH_t elementBVH;

	LinearPath _brushPath;

	TSet<FGraphNodeHandle> _activeElements;

	TMap<FGraphNodeHandle, uint32_t> segmentForHandle;

public:
	virtual ~UMolecularLegoGenerator() {}

	virtual void init(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void tick(float deltaT) override;

	virtual bool wantsFlex()  override { return true; }

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

	virtual std::vector<UClass *> dependencies();

protected:
	void _initPath();
	void _tick(float deltaT);

	bool _hasNearestBrushSegment(FVector point, float radius, size_t& segment_a, size_t& segment_b);

	void _buildElementBVH();

	UGraphSimulationManager * exampleSimulationManager();
	std::vector<FGraphNodeHandle> exampleSelection();

	FGraphNodeHandle _copyElement(FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph, FGraph& targetGraph, const FVector position, const FQuat rotation);

	void _insertElement(FGraphNodeHandle handle, FVector position, float radius, uint32_t segment);

	// Inset the simulation bounds by the particle radius, otherwise we get crazy Flex artifacts.
	FBox _insetBounds();
};