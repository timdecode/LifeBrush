// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include <vector>
#include <set>
#include <memory>

#include "Simulation/FlexGraphSimulation_interface.h"

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/ObjectSimulation.h"
#include "ShipEditorSimulation/GraphSnapshot.h"
#include "ShipEditorSimulation/MeshSimulation.h"

#include "Visualization/EdgeFactory.h"

#include "Algorithm/FastBVH/BVH.h"

#include "Timeline.generated.h"

class UEdgeFactory;

UCLASS(BlueprintType, DefaultToInstanced)
class LIFEBRUSH_API USEGraphEvent : public UObject
{
	GENERATED_BODY()

public:
	virtual ~USEGraphEvent() {}

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	int32 frame = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	FGraphNodeHandle triggeringAgent;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	TArray<FGraphNodeHandle> otherAgents;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	FVector position;

	UInstancedStaticMeshComponent * _transientISMC = nullptr;
	float _transientScale = 1.0f;
};




UCLASS(BlueprintType, DefaultToInstanced)
class LIFEBRUSH_API USEGraphFrame : public UObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	int32 number = 0;

	UPROPERTY(Transient)
	TArray<USEGraphEvent*> events;
};





UCLASS(BlueprintType, DefaultToInstanced)
class LIFEBRUSH_API USEGraphTimeline : public UObject
{
	GENERATED_BODY()

public:
	UPROPERTY(Transient)
	TMap<int32, USEGraphFrame*> sparseFrames;
};






USTRUCT(BlueprintType)
struct LIFEBRUSH_API FEventGlyph 
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	UStaticMesh * glyphMesh = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	UMaterialInterface * glyphMaterial = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	float glyphScale = 0.05f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	TSubclassOf<USEGraphEvent> eventClass;
};






UCLASS(BlueprintType, DefaultToInstanced)
class LIFEBRUSH_API UTimelineSimulation : public UObjectSimulation, public NodeListener, public IFlexGraphSimulation
{
	GENERATED_BODY()

	UTimelineSimulation()
	{
		timeline = CreateDefaultSubobject<USEGraphTimeline>(TEXT("timeline"));
	}

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	FEventGlyph defaultGlyph;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	TArray<FEventGlyph> glyphPrototypes;

	UPROPERTY(Instanced)
	USEGraphTimeline * timeline;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	int32 maxEventsHistory = 900; 

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Timeline")
	int32 maxPathHistory = 900;

protected:
	virtual void attach() override;

public:
	virtual void clear() override;

	virtual void tick(float deltaT) override;

	virtual void snapshotToActor(AActor * actor) override;

	virtual void preTick(float deltaT) override;
	virtual void postTick(float deltaT) override;

	void beginFrame();

	void endFrame();

	template<typename EventType>
	EventType& recordEvent(const FVector& position)
	{
		_didGenerateEventsInFrame = true;

		USEGraphEvent * graphEvent = NewObject<EventType>(_currentFrame);

		_currentFrame->events.Emplace(graphEvent);

		graphEvent->frame = _currentFrame->number;
		graphEvent->position = position;

		return static_cast<EventType&>(*graphEvent);
	}

	std::vector<USEGraphEvent*> eventsOverlappingPosition(const FVector position, float radius);

	void traceEvents(std::vector<USEGraphEvent*> events, int maxHops = 1);

	// -1, use maxFramesBack to find minFrame
	void showAllEvents();

	void setGlyphVisibility(bool visibility);

protected:
	void _clearVisualizations();

	void _loadEvent(USEGraphEvent& graphEvent);

	void _loadFrameEvents(USEGraphFrame& frame);

	void _loadAllFrameEvents(int32 minFrame = 0);

	void _initGlyphIndicesForEventClass();

	// _initGlyphIndicesForEventClass must be called first
	void _initISMCs();

	void _expandFrame(USEGraphFrame * frame, TSet<FGraphNodeHandle> &agents, TSet<USEGraphEvent *> &includedEvents);

	bool _parseFrame(USEGraphFrame * frame, TSet<FGraphNodeHandle>& agents_out, TSet<USEGraphEvent*>& events_out);

	void _expandFrames(
		USEGraphEvent * graphEvent,
		const std::vector<int32>& sortedFrames,
		TSet<FGraphNodeHandle> &agents_out,
		TSet<USEGraphEvent *> &includedEvents_out, 
		int32 maxFrameDistance
	);
	void _expandFrames2(
		USEGraphEvent * graphEvent,
		const std::vector<int32>& sortedFrames,
		TSet<FGraphNodeHandle> &agents_out,
		TSet<USEGraphEvent *> &includedEvents_out,
		int32 maxFrameDistance
	);

	void _showAgentPathLines(TSet<FGraphNodeHandle>& agentSet, TSet<USEGraphEvent*>& eventSet);


	FEventGlyph& glyphForEventClass(TSubclassOf<USEGraphEvent> eventClass);

	UInstancedStaticMeshComponent& meshForEvent(USEGraphEvent& graphEvent);


protected:
	UPROPERTY(Transient)
	USEGraphFrame * _currentFrame;

	// will contain [class,-1] if there is no glyph for an event class
	std::unordered_map<UClass*, int32> _glyphIndicesForEventClass;

	std::unordered_map<UClass*, UInstancedStaticMeshComponent*> _instancedStaticeMeshes;

	std::set<USEGraphEvent*> _loadedEvents;

	bool _cached_showAgentPathLines = false;

	bool _didGenerateEventsInFrame = false;
};



USTRUCT(BlueprintType)
struct LIFEBRUSH_API FPositionSnapshot
{
	GENERATED_BODY()
public:
	TArray<FVector> agentPositions;

	TBitArray<> validPositions;
	
	int32 frameNumber = 0;
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UVisualization_AgentPathLines : public UObjectSimulation, public NodeListener
{
	GENERATED_BODY()

protected:
	virtual void attach();

public:
	virtual void clear();

	virtual void tick(float deltaT) override;

	virtual void tick_paused(float deltaT) override;

	virtual void snapshotToActor(AActor * actor) override;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	UMaterialInterface * defaultLineMaterial = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	bool showAgentPathLines = false;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	int32 historySize = 30;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	float pathRadius = 0.2f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	float secondsBetweenSegments = 0.1f;

	UPROPERTY()
	UEdgeFactory * _dynamicEdgeFactory = nullptr;

	UPROPERTY()
	UColoredLineFactory * _agentPathFactory;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	bool hidePaths = true;

	USceneComponent * camera;

	FTransform toWorld;

	void buildBVH(ColoredQuadFactory& quadFactory);

	std::unordered_map<UMaterialInterface*, FColoredLineBuilder> coloredLinesForAgents(const std::vector<FGraphNodeHandle>& agents, float radius, int32 minFrame, int32 maxFrame);

	void showTotalHistoryForAgents(const std::vector<FGraphNodeHandle>& agents, int32 minFrame, int32 maxFrame);
	void hideTotalHistory();

	void captureFrame();

public:
	virtual void nodeAdded(FGraphNode& node) override;

	virtual void nodeRemoved(FGraphNode& oldNode) override;

protected:
	void _initPipeFactory();

	void _initCache();

	void _tickVisibility(float deltaT);

	void _clearDynamicEdgeFactory();

	void _clearAgentPathFactory();

	void _updateDynamicPaths();

	void _cachePositions();

	void _fadeBackground(bool fade);

	void _setOtherActorsDesaturation(bool desaturated);

	int32 _sectionForMaterial(UMaterialInterface* material);

	FLinearColor _colorForNode(FGraphNodeHandle node);


protected:
	struct CachcedPosition
	{
		FVector position;
		FVector velocity;
	};

	std::unordered_map<UMaterialInterface*, int32> _materialToSection;

	std::unordered_map<FGraphNodeHandle, bool> _cachedVisibility;

	TArray<FVector> lastPositionCache;

	bvh::BVH pathBVH;

	// No, let's not persist this.
//	UPROPERTY()
	TArray<FPositionSnapshot> totalHistory;

	std::unordered_set<FGraphNodeHandle> _totalHistoryAgents;

	float _timer = 0.0f;

	int32 historyIndex = 0;

	UPROPERTY()
	UTexture2D * _valuesTexture = nullptr;
};