// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ElementGenerator.h"

#include "Simulation/FlexGraphSimulation_interface.h"
#include "ShipEditorSimulation/Graph.h"
#include "aabbcc/unrealAABB.h"
#include "ElementActor.h"
#include "MolecularLego/MolecularLego.h"

#include "SwarmGenerator.generated.h"

// -----------------------------------------------------------
// Swarm Grammar Prototypes
// -----------------------------------------------------------

UCLASS()
class LIFEBRUSH_API ASGPrototypeActor : public AElementActor
{
	GENERATED_BODY()

public:

	virtual void writeToElement(ElementTuple& element, FGraph& graph) override;

public:
	FGraphNodeHandle handleInRuleGraph;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName typeName;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSGRuleTypeName : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	FName typeName;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSGFilamentRule : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float segmentLength = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString headSequence;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString bodySequence;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	int32 numInBody = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString tailSequence;

	UPROPERTY()
	TArray<FGraphNodeHandle> headHandles;

	UPROPERTY()
	TArray<FGraphNodeHandle> bodyHandles;

	UPROPERTY()
	TArray<FGraphNodeHandle> tailHandles;
};

// Attach this to a node that you want to be the start of a filament. The rule will be instantiated
// and the FBoidGenerator will set it's FBoidGenerator::last node handle to this node. 
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSGFilamentAnchor : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName filamentRuleName;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSGStarBranchEntry
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FVector direction;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName typeName;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSGStarRule : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TArray<FSGStarBranchEntry> branches;

	// transient cache of names converted to handles
	// this is a parallel array for FSGStarProtoype::prototypes
	TArray<FGraphNodeHandle> handlesForBranches;
};

// -----------------------------------------------------------
// Swarm Grammar Executors
// -----------------------------------------------------------

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidGenerator : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	int currentIndex = 0;

	UPROPERTY()
	int maxCount = 0;

	UPROPERTY()
	FGraphNodeHandle last = FGraphNodeHandle::null;
	
	UPROPERTY()
	FGraphNodeHandle last2ring = FGraphNodeHandle::null;

	UPROPERTY()
	FGraphNodeHandle starNode = FGraphNodeHandle::null;

	UPROPERTY()
	int32 filamentGroup = -1;

	UPROPERTY()
	float filamentRadius = 0.5f;

	UPROPERTY()
	float segmentLength = 0.5f;

	UPROPERTY()
	FGraphNodeHandle currentExemplar = FGraphNodeHandle::null;

	UPROPERTY()
	float currentDistance = 0.0f;

public:
	// A handle to the prototype in the exemplar
	UPROPERTY()
	FGraphNodeHandle _exemplarPrototype;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidSegment : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 0.5f;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidStar : public FGraphObject
{
	GENERATED_BODY()

public:
	// A handle to the prototype in the exemplar
	UPROPERTY()
	FGraphNodeHandle _exemplarPrototype;

	UPROPERTY()
	bool _didApply = false;

	UPROPERTY()
	int32 targetRingSize = 0;

	UPROPERTY()
	TArray<FGraphNodeHandle> _inProgressRing;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoid : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;
};

// All boids instantiated by a grammar will share the same root. 
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidRootNode : public FGraphObject
{
	GENERATED_BODY()

public:
	FBoidRootNode() {}
	FBoidRootNode(FGraphNodeHandle rootNode) : rootNode(rootNode) {}

	// The node that created the boid structure
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FGraphNodeHandle rootNode = FGraphNodeHandle::null;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidSeeker : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	FGraphNodeHandle starNode = FGraphNodeHandle::null;

	// Species that we can bind with
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TArray<FName> bobSpecies;
};

// Add this to a node to prevent further bindings from an FBoidSeeker
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoidSeekerBlocker : public FGraphObject
{
	GENERATED_BODY()

public:
};


UENUM(BlueprintType)
enum class ESwarmGenerator_BrushType : uint8
{
	Star UMETA(DisplayName = "Star"),
	Anchor UMETA(DisplayName = "Anchor"),
};

UCLASS(BlueprintType)
class LIFEBRUSH_API USwarmSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FGraph ruleGraph;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float brushSpacing = 7.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float seekerRadius = 3.0f;
	
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius_segment = 1.8f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")		
	float scaleFactor_segment = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material_segment = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh_segment = nullptr;

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
	// Misc
	// -------------------------------------------------
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool debugDrawBoids = true;



public:
	virtual void attach() override;
	virtual void begin() override;
	virtual void detach() override;

	virtual void tick(float deltaT) override;


public:
	void addBrushPoint(FVector point, ESwarmGenerator_BrushType type);

	FGraphNodeHandle _createAnchor(FVector position, FQuat orientation);

	void _applyStar(FGraphNodeHandle starHandle);

	FQuat randomQuat();

	void _updateBVH();

	void setBrushPrototypoe(FGraphNodeHandle handle);

	FGraphNodeHandle copyElement(FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph, FGraph& targetGraph, const FVector position, const FQuat rotation_in, FGraphNodeHandle boidRootNode);

protected:
	void _tickFilamentAnchors(float deltaT);

	void _tickBoidStars(float deltaT);

	void _tickBoidGenerators();
	FGraphNodeHandle _handleInFilamentSequence(int32 sequenceIndex, FGraphNodeHandle ruleNode);
	
	void _tickBoidSeekers(float deltaT);

	void _connectSeekers(FGraphNodeHandle a, FGraphNodeHandle b, float distance);
	void _connectSeekerOneWay(FGraphNodeHandle seeker, FGraphNodeHandle target, float distance);
	void _cleanupConnection(FGraphNode& aNode, FGraphNode& bNode);

	void _attachSegmentMeshes();
	void _detachSegmentMeshes();

	void _initRuleGraph();

	void _readRuleGraph();

	void _readStars();

	TArray<FGraphNodeHandle> _parseSequence(FString sequence);



public:
	typedef unrealAABB::Tree BVH_t;

protected:

	FRandomStream rand;

	struct BrushPoint
	{
		FVector point;
		ESwarmGenerator_BrushType type;
	};
	std::vector<BrushPoint> _brushPoints;

	FGraphNodeHandle _brushPrototype;

	TMap<FName, FGraphNodeHandle> _nameToRuleHandles;

public:
	// Spatial index for NC1 heads and 7S tails, but not collagen segments.
	BVH_t _spaceBVH;
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

	void setPrototype(FGraphNodeHandle exampleHandle);

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ESwarmGenerator_BrushType brushType;



protected:
	SynthesisContext * _context;

	FGraphNodeHandle _brushExampleHandle = FGraphNodeHandle::null;

	FGraph* _graph;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

public:
};