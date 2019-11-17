//
//  LifeBrush
//
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "ElementGenerator.h"

#include "Algorithm/ctpl_stl.h"
#include "Algorithm/Algorithm_RegionGrowing.h"
#include "RegionGrowingComponent.h"

#include "RegionGrowingGenerator.generated.h"

UCLASS(ClassGroup = (Custom), DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API URegionGrowingGenerator : public UElementGenerator
{
	GENERATED_BODY()

protected: // Threading
	ctpl::thread_pool _generationWorker;

	std::atomic<bool> _generationWorkFinished = { true };
	std::function<void()> _generationWorkerMainThreadWork;

protected: // Core data structures
	std::unique_ptr<Algorithm_RegionGrowing> _algorithm;

	SynthesisContext * _outputContext;

	// The simulation manager that we use to draw the scene
	UPROPERTY()
	UGraphSimulationManager * _outputSimulationManager;

	UPROPERTY()
	UGraphSimulationManager * _exemplarSimulationManager;

protected:
	// Cache of the points to remove
	std::vector< PositionRadiusFace > _toRemove;

	std::unordered_map<AElementActor*, FGraphNodeHandle> _actorToElement;


public:
	virtual ~URegionGrowingGenerator() {}

	virtual void init(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void tick(float deltaT) override;

	virtual void start() override;
	virtual void stop() override;

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

	virtual void beginEraserPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void eraseInRadiusAt(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endEraserPath() override;

	void setExampleSelection(std::vector<AElementActor*> elementActors);
	std::vector<FGraphNodeHandle> exampleSelection();

	FGraphNodeHandle handleForExampleActor(AElementActor * elementActor);

	UGraphSimulationManager * exampleSimulationManager();


	void setGenerationMode(EGenerationMode generationMode);
	EGenerationMode generationMode();

public:
	virtual FString parametersString();

protected:
	void loadParameters();
	void syncExemplarFromElementActors();

	FString _optimizationParametersString(FNeighbourhoodParameters& params);

public:
	UPROPERTY(EditAnywhere) bool perElementParameters = false;

	UPROPERTY(EditAnywhere) bool pauseSynthesis = false;

	UPROPERTY(EditAnywhere) bool enableRoundSummaries = false;
	UPROPERTY(EditAnywhere) bool calculate_kCoherenceEnergy = false;
	UPROPERTY(EditAnywhere) bool calculate_bruteForceEnergy = false;
	UPROPERTY(EditAnywhere) bool useCinpactEnergy = false;
	UPROPERTY(EditAnywhere) float cinpactCellSize = 1.0f;

	// A value of 1.0 will not relax overlaps
	// A value < 1.0 will increase the amount of overlap allowed
	// The comparison in Algorithm::_overlaps is something like
	//     bool overlaps = e.position - ex.position < (e.radius + ex.radius) * relaxation;
	UPROPERTY(EditAnywhere) float relaxation = 0.8f;

	UPROPERTY(EditAnywhere) FNeighbourhoodParameters generationParameters; // search parameters
	UPROPERTY(EditAnywhere) FNeighbourhoodParameters optimizationParameters;

	UPROPERTY(EditAnywhere) bool useTypeVoting = true;
	UPROPERTY(EditAnywhere) bool removeOverlaps = false;

	// An expensive option to reassign the source example for each element along a brush path.
	UPROPERTY(EditAnywhere) bool enableReassignment = false;

	UPROPERTY(EditAnywhere) float generationRadius; // the radius out to which elements are generated around the seeds
	UPROPERTY(EditAnywhere) float frozenElementRadius = 0.0f;

	UPROPERTY(EditAnywhere) int32 kCoherence = 5;

	UPROPERTY(EditAnywhere) float brushSize = 30.0f;

	UPROPERTY(EditAnywhere) float selectionRadius = 30.0f;

	UPROPERTY(EditAnywhere) bool enableVolumeSurfaceInteraction = true;


	UPROPERTY(EditAnywhere) float voidSize = 10.0f; // the size of the void brush (the void brush prevents elements from being generated in the painted regions)

	// Used in the closest point assignment, when points are too far from each other
	// then they are not paired
	UPROPERTY(EditAnywhere) float minAssignmentDistance = 50.0f;
	UPROPERTY(EditAnywhere) float freespaceRadius = 15.0f;
	UPROPERTY(EditAnywhere) float typeCost = 30.0f;
	UPROPERTY(EditAnywhere) float gradientTerm = 1.0f;

	UPROPERTY(EditAnywhere) float kCoherenceClusteringPenalty = 0.0f;
	UPROPERTY(EditAnywhere) float kCoherenceClusteringPenaltyRadius = 5.0f;

	UPROPERTY(EditAnywhere) bool ignoreBadSuggestions = false;
	UPROPERTY(EditAnywhere) float ignoreBadSuggestionsDistanceFactor = 10.0f;

	UPROPERTY(EditAnywhere) float sourceHistogramRadius = 0.0f;
	UPROPERTY(EditAnywhere) float sourceHistogramWeight = 0.0f;
	UPROPERTY(EditAnywhere) float samplingDistanceWeight = 0.0f;

	UPROPERTY(EditAnywhere) float stretchRadius = 100.0f;

	UPROPERTY(EditAnywhere) bool sketchBasedForceTangentPlane = false;

	UPROPERTY(EditAnywhere) bool flipSurfaceNormals = false;

	UPROPERTY(EditAnywhere) bool disableOptimization = false;
	UPROPERTY(EditAnywhere) int optimizationRounds = 1;
	UPROPERTY(EditAnywhere) bool useGlobalOptimization = false;

	UPROPERTY(VisibleAnywhere, BlueprintReadOnly) int32 statElementCount = 0;

protected:
	UPROPERTY(EditAnywhere) EGenerationMode _generationMode;
};