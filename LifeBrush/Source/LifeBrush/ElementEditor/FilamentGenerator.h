//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#pragma once

#include "ElementGenerator.h"

#include "ShipEditorSimulation/ObjectSimulation.h"
#include "LinearPath.h"
#include "ElementEditor/SwarmGenerator.h"

#include "FilamentGenerator.generated.h"

// Add this graph object to a node to create a new filament segment species.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFilamentPrototypeSegment : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString species;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FFilamentPrototype
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float segmentLength = 1.0f;

	// The head sequence. The sequence is comma separated, the words are FFilamentPrototypeSegment::species names.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString headSequence;

	// The body sequence.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString bodySequence;

	// The end sequence.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FString terminusSequence;

	// The number of body segments to generate from the body sequence when bodyLengthType == FromNumSegments.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	int numSegments = 0;

	// The colour of the segments
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * bodyMaterial = nullptr;
};

UCLASS()
class LIFEBRUSH_API AFilamentPrototypeActor : public AStaticMeshActor
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FFilamentPrototype filamentPrototype;
};

UCLASS(ClassGroup = (Custom), DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UFilamentGenerator : public UElementGenerator
{
	GENERATED_BODY()

public:
	SynthesisContext * _context;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float brushSegmentLength = 0.5f;

public:
	virtual ~UFilamentGenerator() {}

	virtual void init(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void tick(float deltaT) override;
	virtual void tickPaused(float deltaT) override;

	float _correctedSegmentLength(float radius);
	virtual bool wantsFlex()  override { return true; }

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;

	void _beginBrushPath(FVector localPoint, float radius, FSurfaceIndex surfaceIndex);

	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

public:
	void setFilamentPrototype(FFilamentPrototype filamentPrototype);
	void clearFilamentPrototype();

	void togglePlay();

protected:
	void _readSegments();

	FGraph * _exampleGraph();

	FGraphNodeHandle _copyElement(
		FGraphNodeHandle sourceNodeHandle, 
		FGraph& sourceGraph, 
		FGraph& targetGraph, 
		const FVector position,
		const FQuat rotation_in,
		FGraphNodeHandle rootNode);

	// parse the sequence into example graph node handles
	TArray<FGraphNodeHandle> _parseSequence(FString sequence);

	void _initPath();

	float _sequenceLength(int32 numInSequence);

	void _updateFilamentPositions();
	void _updateConnections();

	FGraphNodeHandle _nextPrototypeHandle();

	// Inset the simulation bounds by the particle radius, otherwise we get crazy Flex artifacts.
	FBox _insetBounds();

	float _filamentLength();

protected:
	FRandomStream rand;

	bool _hasPrototype = false;
	FFilamentPrototype _prototype;

	TArray<FGraphNodeHandle> _filament;
	TArray<FQuat> _filamentOriginalRotations;

	bool _forceEndPath = false;
	bool _firstOutside = false;

	int32 _filamentGroup = 0;

	TMap<FName, FGraphNodeHandle> _segmentTypeToNode;

	UGraphSimulationManager * _exampleSimulationManager = nullptr;

	// the path
	LinearPath _brushPath;

	TArray<FGraphNodeHandle> _headSequence;
	TArray<FGraphNodeHandle> _tailSequence;
	TArray<FGraphNodeHandle> _bodySequence;
};