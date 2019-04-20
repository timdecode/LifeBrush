//
//  Created by Timothy Davison on 2018-12-29.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "Algorithm/PointCloud.hpp"

#include "ElementGenerator.h"
#include "Algorithm/Algorithm.h"
#include "aabbcc/AABB.h"

#include "CollagenGenerator.generated.h"

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FCollagen_NC1Head : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	bool satisfiable = true;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FCollagen_7STail : public FGraphObject
{
	GENERATED_BODY()

public:
	auto bound(FGraph& graph) -> bool;
	auto collagenTails(FGraph& graph) ->std::pair<FGraphNodeHandle, FGraphNodeHandle>;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FCollagen_Segment : public FGraphObject
{
	GENERATED_BODY()

public:

};


UCLASS(ClassGroup = (Custom), DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UCollagenGenerator : public UElementGenerator
{
	GENERATED_BODY()

public:
	SynthesisContext * _context;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float chainLength = 8.0f;

	// -------------------------------------------------
	// 7S Head
	// -------------------------------------------------

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius_7S = 1.0f; 

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float scaleFactor_7S = 0.02;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_7S = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_7S = nullptr;

	// -------------------------------------------------
	// NC1 Head
	// -------------------------------------------------

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius_NC1 = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float scaleFactor_NC1 = 0.02;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_NC1 = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_NC1 = nullptr;

	// -------------------------------------------------
	// Segments
	// -------------------------------------------------

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius_segment = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float scaleFactor_segment = 0.02;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_segment = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_segment = nullptr;

protected:
	typedef aabb::Tree<3, float> BVH_t;

	PositionFacePairCloud _brushPoints;
	std::unique_ptr<PositionFacePairCloudAdaptor> _index;

	FVector start;


	BVH_t _brushBVH;

	// Spatial index for NC1 heads and 7S tails, but not collagen segments.
	BVH_t _spaceBVH;

	std::set<FGraphNodeHandle> _frontier;

	FRandomStream rand;

public:
	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void tick(float deltaT) override;

	virtual bool wantsFlex() override { return true; }

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

protected:
	FQuat _randomQuat();

	FGraphNode& _createNC1Head(FVector position, FQuat orientation);
	FGraphNode& _create7STail(FVector position, FQuat orientation);
	FGraphNode& _createSegment(FVector position);

	// Conditionally creates a NC1 Head, if there are none nearby
	FGraphNodeHandle _conditionallySeedNC1Head(FVector position);

	void _loadFrontier();

	void _link(FGraphNodeHandle start, FGraphNodeHandle end);

	FQuat _rotationAroundY(float rad);

	// Does the point overlap the brush points?
	bool _inBounds(FVector p);

	bool _hit(FVector occupied, float radius);

	void _spaceColonization();

	bool _isSatisfied(FGraphNodeHandle handle);
	bool _isSatisfiable(FGraphNodeHandle handle);

	void _process_nc1Head(FGraphNode node);
	void _process_7sTail(FGraphNode node);

	std::vector<FGraphNodeHandle> _chainsOut(FGraphNode& node, FCollagen_NC1Head& nc1);
};