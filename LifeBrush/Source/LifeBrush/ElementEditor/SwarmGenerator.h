// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ElementGenerator.h"

#include "Simulation/FlexGraphSimulation_interface.h"
#include "ShipEditorSimulation/Graph.h"
#include "aabbcc/AABB.h"

#include "SwarmGenerator.generated.h"


USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidGenerator : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	int counter = 0;

	UPROPERTY()
	int maxCount = 0;

	UPROPERTY()
	FGraphNodeHandle last = FGraphNodeHandle::null;
	
	UPROPERTY()
	FGraphNodeHandle last2ring = FGraphNodeHandle::null;

	UPROPERTY()
	FGraphNodeHandle starNode = FGraphNodeHandle::null;

	UPROPERTY()
	uint32 filamentGroup = 0;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidSegment : public FGraphObject
{
	GENERATED_BODY()

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidStar : public FGraphObject
{
	GENERATED_BODY()

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoid : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	float radius = 1.0f;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidSeeker : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	FGraphNodeHandle starNode = FGraphNodeHandle::null;
};




UCLASS(BlueprintType)
class LIFEBRUSH_API USwarmSimulation : public UObjectSimulation, public IFlexGraphSimulation
{
	GENERATED_BODY()

public:
	float minNC1Separation = 6.0f;
	float radius_segment = 1.0f;
	float scaleFactor_segment = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_segment = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_segment = nullptr;


	float radius_7S = 1.0f;
	float scaleFactor_7S = 0.02;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_7S = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_7S = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool debugDrawBoids = true;

public:
	virtual void begin() override;
	virtual void detach() override;

	virtual void tick(float deltaT) override;

	virtual void flexTick(
		float deltaT,
		NvFlexVector<int>& neighbourIndices,
		NvFlexVector<int>& neighbourCounts,
		NvFlexVector<int>& apiToInternal,
		NvFlexVector<int>& internalToAPI,
		int maxParticles
	) override;

protected:
	void _tickBoidStars(float deltaT);
	void _tickBoidGenerators();
	void _tickBoidSeekersFlexNeighbours(int maxParticles, NvFlexVector<int>& apiToInternal, NvFlexVector<int>& neighbourCounts, NvFlexVector<int>& internalToAPI, NvFlexVector<int>& neighbourIndices);
	void _tickBoidSeekers(float deltaT);

	FGraphNodeHandle _createSegment(FVector position, FQuat orientation);
	FGraphNodeHandle _createSeeker(FVector position, FQuat orientation, FGraphNodeHandle starNode);

	void _loadStarBVH();

	void _attachSegmentMeshes();
	void _detachSegmentMeshes();

protected:
	typedef aabb::Tree<3, float> BVH_t;

	// Spatial index for NC1 heads and 7S tails, but not collagen segments.
	BVH_t _starBVH;
};

UENUM(BlueprintType)
enum class ESwarmGenerator_BrushType : uint8
{
	Star UMETA(DisplayName = "Star"),
	Anchor UMETA(DisplayName = "Anchor"),
};

UCLASS( BlueprintType )
class LIFEBRUSH_API USwarmGenerator : public UElementGenerator, public IFlexGraphSimulation
{
	GENERATED_BODY()

public:
	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void flexTick(
		float deltaT,
		NvFlexVector<int>& neighbourIndices,
		NvFlexVector<int>& neighbourCounts,
		NvFlexVector<int>& apiToInternal,
		NvFlexVector<int>& internalToAPI,
		int maxParticles
	);

	virtual bool wantsFlex() override { return true; }

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ESwarmGenerator_BrushType brushType;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float minNC1Separation = 6.0f;

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
	FGraphNodeHandle _createAnchor(FVector position, FQuat orientation);

	FGraphNodeHandle _createStar(FVector position, FQuat orientation);
	FGraphNodeHandle _createWanderer(FVector position, FQuat orientation);

	FQuat _randomQuat();

	void _reloadBVH();

protected:
	SynthesisContext * _context;

	FGraph* _graph;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

	FRandomStream rand;




	typedef aabb::Tree<3, float> BVH_t;

	// Spatial index for NC1 heads and 7S tails, but not collagen segments.
	BVH_t _spaceBVH;
};