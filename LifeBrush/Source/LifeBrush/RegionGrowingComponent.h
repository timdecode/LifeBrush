// Copyright (c) 2017 Timothy Davison. All rights reserved.

#pragma once

#include "Components/ActorComponent.h"
#include "Algorithm//SynthesisContext.h"
#include "Algorithm/Algorithm.h"
#include "Algorithm/Algorithm_RegionGrowing.h"
#include "Algorithm/Algorithm_PatchCopy.h"
#include "Algorithm/ElementKDTree.h"

#include "Algorithm/ctpl_stl.h"

#include <map>
#include <memory>
#include <ctime>

#include "tcodsMeshInterface.h"
#include "DepthReader.h"
#include "AlgorithmEnums.h"
#include "InstanceManager.h"
#include "NeighbourhoodParameters.h"
#include "ElementActor.h"
#include "MeshInterfaceMode.h"

#include "RegionGrowingComponent.generated.h"

class URuntimeMeshComponent;
class AElementActor;

namespace tcods { class Mesh; }

static const FName TCODSMarker = FName(TEXT("TCODSMarker"));

struct OutputDomainExport
{
public:
	OutputDomainExport()
	{
		graph.init();
	}

	void read(FGraph& graph_in, FTransform transform = FTransform::Identity)
	{
		graph.clear();

		// copy all  the nodes
		auto sourceNodesVec = graph_in.nodeHandles();
		TArray<FGraphNodeHandle> sourceNodes(&sourceNodesVec[0], sourceNodesVec.size());

		// copy all the edges
		auto sourceEdgesVec = graph_in.edgeHandles();
		TArray<FGraphEdgeHandle> sourceEdges(&sourceEdgesVec[0], sourceEdgesVec.size());

		FGraphCopyContext::copySubgraph(graph_in, graph, sourceNodes, sourceEdges, transform);
	}

	void readSubgraph(FGraph& sourceGraph, TArray<FGraphNodeHandle>& nodes, TArray<FGraphEdgeHandle>& edges, FTransform transform = FTransform::Identity)
	{
		graph.clear();

		FGraphCopyContext::copySubgraph(sourceGraph, graph, nodes, edges, transform);
	}

	void write(FGraph& graph_in, FTransform transform = FTransform::Identity)
	{
		// copy all  the nodes
		auto sourceNodesVec = graph.nodeHandles();
		TArray<FGraphNodeHandle> sourceNodes(&sourceNodesVec[0], sourceNodesVec.size());

		// copy all the edges
		auto sourceEdgesVec = graph.edgeHandles();
		TArray<FGraphEdgeHandle> sourceEdges(&sourceEdgesVec[0], sourceEdgesVec.size());

		FGraphCopyContext::copySubgraph(graph, graph_in, sourceNodes, sourceEdges, transform);
	}

	FGraph graph;

	std::shared_ptr<tcodsMeshInterfaceBase> meshInterface;
};



USTRUCT()
struct FSingularity
{
    GENERATED_USTRUCT_BODY()
    
    FSingularity() {}
    FSingularity(uint32 vertexIndex_, float k_) { vertexIndex = vertexIndex_; k = k_; }
    FSingularity(const FSingularity& other) { vertexIndex = other.vertexIndex; k = other.k; }
    
    UPROPERTY(EditAnywhere) uint32 vertexIndex = 0;

	UPROPERTY(EditAnywhere) uint32 sectionIndex = 0;

    UPROPERTY(EditAnywhere) float k = 0.0f;
};

UENUM( BlueprintType )
enum class EAlgorithm : uint8
{
	RegionGrowing UMETA(DisplayName="Region Growing"),
	PatchCopy UMETA(DisplayName="Patch Copy"),
};

class RGOcclusionTester : public OcclusionBase
{
public:
	virtual ~RGOcclusionTester() {}

	virtual bool isVisible( PositionFace point, float radius ) override;


public:
    UDepthReader * reader = nullptr;
    FVector cameraLocation;
	std::shared_ptr<tcodsMeshInterfaceBase> _meshInterface = nullptr;
};



USTRUCT()
struct FElementType
{
	GENERATED_BODY()

public:
	FElementType() {}

	FElementType( AActor * actorTemplate )
	{

	}

	FElementType( UStaticMesh * mesh, UMaterialInterface * material )
	{

	}

	UPROPERTY() UStaticMesh * mesh = nullptr;
	UPROPERTY() UMaterialInterface * material = nullptr;

	UPROPERTY() int16 integerType = 0;



	struct InferredElementComparator
	{
		bool operator()( const FElementType& a, const FElementType& b ) const
		{
			if(a.mesh != b.mesh)
				return a.mesh < b.mesh;

			if(a.material != b.material)
				return a.material < b.material;

			return false;
		}
	};
};



USTRUCT()
struct FElementTypeLookup
{
	GENERATED_BODY()

public:
	int16 getOrGenerateType( AActor * actor )
	{
		FElementType query( actor );
		
		auto found = typeSet.find( query );

		if(found == typeSet.end())
		{
			auto type = elementTypes.Num();

			query.integerType = type;

			typeSet.insert( query );

			return type;
		}
		else
			return found->integerType;
	} 

	void PostSerialize( const FArchive& Ar )
	{
		for(auto type : elementTypes)
			typeSet.insert(type);
	}

protected: 
	std::set<FElementType, FElementType::InferredElementComparator> typeSet;

	UPROPERTY() TArray<FElementType> elementTypes;
};

template<>
struct TStructOpsTypeTraits<FElementTypeLookup> : public TStructOpsTypeTraitsBase2<FElementTypeLookup>
{
	enum
	{
		WithPostSerialize = true,
	};
};

UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API URegionGrowingComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere ) bool perElementParameters = false;


    UPROPERTY( EditAnywhere ) bool pauseSynthesis = false;

	UPROPERTY( EditAnywhere ) bool enableRoundSummaries = false;
	UPROPERTY( EditAnywhere ) bool calculate_kCoherenceEnergy = false;
	UPROPERTY( EditAnywhere ) bool calculate_bruteForceEnergy = false;
	UPROPERTY( EditAnywhere ) bool useCinpactEnergy = false;
	UPROPERTY( EditAnywhere ) float cinpactCellSize = 1.0f;


	// The duration at which to pause synthesis. If this is negative, it is ignored.
	UPROPERTY( EditAnywhere ) float stopDuration = -1.0f;

	// A value of 1.0 will not relax overlaps
	// A value < 1.0 will increase the amount of overlap allowed
	// The comparison in Algorithm::_overlaps is something like
	//     bool overlaps = e.position - ex.position < (e.radius + ex.radius) * relaxation;
	UPROPERTY( EditAnywhere ) float relaxation = 0.8f;




    UPROPERTY(EditAnywhere) FNeighbourhoodParameters generationParameters; // search parameters
    UPROPERTY(EditAnywhere) FNeighbourhoodParameters optimizationParameters;

	UPROPERTY( EditAnywhere ) bool useTypeVoting = true;
	UPROPERTY( EditAnywhere ) bool removeOverlaps = false;

	// An expensive option to reassign the source example for each element along a brush path.
	UPROPERTY(EditAnywhere) bool enableReassignment = false;

    UPROPERTY(EditAnywhere) float generationRadius; // the radius out to which elements are generated around the seeds
	UPROPERTY(EditAnywhere) float frozenElementRadius = 0.0f;

    UPROPERTY(EditAnywhere) int32 kCoherence = 5;

	UPROPERTY(EditAnywhere) int32 nThreads = 2;
    
    UPROPERTY(EditAnywhere) EGenerationMode generationMode;	

	UPROPERTY(EditAnywhere) EAlgorithm algorithmType = EAlgorithm::RegionGrowing;
    
	UPROPERTY(EditAnywhere) EMeshInterfaceMode meshInterfaceMode = EMeshInterfaceMode::StaticMesh;

	UPROPERTY(EditAnywhere) AActor * chunkedStaticMeshActor = nullptr;

    UPROPERTY(EditAnywhere, BlueprintReadWrite) AActor * exemplar = nullptr;
	UPROPERTY(EditAnywhere) AActor * targetExemplar = nullptr;
   
    UPROPERTY(EditAnywhere) float brushSize = 30.0f;

	UPROPERTY(EditAnywhere) float selectionRadius = 30.0f;

	UPROPERTY(EditAnywhere) bool enableVolumeSurfaceInteraction = true;


	UPROPERTY(EditAnywhere) float voidSize = 10.0f; // the size of the void brush (the void brush prevents elements from being generated in the painted regions)
    
    UPROPERTY(EditAnywhere) float hideSurfaceDistance = 0.0f;
	UPROPERTY(EditAnywhere) float resetDistance = 0.0f;

    // Synthesize elements when they are synthesizeDistance from the camera.
    // Use a negative number to synthesize to infinity on manual stepping.
    UPROPERTY(EditAnywhere) float synthesizeDistance;
    
    UPROPERTY(EditAnywhere) float innerQuellingDistance = 0.0f;
	UPROPERTY(EditAnywhere) float outerQuellingDistance = 0.0f;
    
    UPROPERTY(EditAnywhere) FVector generationLimitsMin = {-1000.0f, -1000.0f, -1000.0f};
    UPROPERTY(EditAnywhere) FVector generationLimitsMax = { 1000.0f,  1000.0f,  1000.0f};
    
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

    UPROPERTY(EditAnywhere) int32 batchSize = 1;
    
    UPROPERTY(EditAnywhere) int32 sampleSubdivisions = 64;
    UPROPERTY(EditAnywhere) bool flipSurfaceNormals = false;

	UPROPERTY( EditAnywhere ) bool forceRotation = false;
	UPROPERTY( EditAnywhere ) FQuat forcedRotation = FQuat::Identity;
    
    UPROPERTY(EditAnywhere) float sourceHistogramRadius = 0.0f;
    UPROPERTY(EditAnywhere) float sourceHistogramWeight = 0.0f;
    UPROPERTY(EditAnywhere) float samplingDistanceWeight = 0.0f;
    
    UPROPERTY( EditAnywhere ) float stretchRadius = 100.0f;
    
    UPROPERTY( EditAnywhere ) bool rayCastSeeds = true;
    UPROPERTY( EditAnywhere ) float rayCastDivisions = 9.0f;

	UPROPERTY( EditAnywhere ) bool sketchBasedForceTangentPlane = false;

	UPROPERTY( EditAnywhere ) FIntVector patchInitializationExemplarDivisions = FIntVector( 1, 1, 1 ); 
	// The percentage of the calculate patch extents to over-copy elements (leads to overlaps in the margins around patches, which is necessary for some exemplars).
	UPROPERTY( EditAnywhere ) float patchMarginPercentage = 0.0f;

    UPROPERTY(EditAnywhere) AActor * depthReaderActor = nullptr;
    UDepthReader * depthReader()
    {
        if( !depthReaderActor )
            return nullptr;
        
        return (UDepthReader*)depthReaderActor->GetComponentByClass(UDepthReader::StaticClass());
    }
    
    // singularities and seeds are the in the coordinate space of the untransformed mesh
    UPROPERTY(EditAnywhere) TArray< FSingularity > singularities;

	UPROPERTY( EditAnywhere ) bool useSeeds = false;
	UPROPERTY( EditAnywhere ) bool seedsIgnoreMeshInBoundary = false;
    UPROPERTY(EditAnywhere) TArray< FVector > seeds;
    
    UPROPERTY(EditAnywhere) AActor * _output = nullptr;
    AActor * output()
    {
        if( _output == nullptr )
        {
            _output = GetWorld()->SpawnActor<class AActor>(this->GetOwner()->GetClass(), FVector::ZeroVector, FRotator::ZeroRotator, FActorSpawnParameters());
            _output->SetRootComponent(NewObject<USceneComponent>(_output));
            _output->GetRootComponent()->SetMobility(EComponentMobility::Static);
			_output->AttachToActor( this->GetOwner(), FAttachmentTransformRules::KeepRelativeTransform );
#if WITH_EDITOR
            _output->SetActorLabel(FString(TEXT("Output")));
#endif
        }
        
        return _output;
    }
	UPROPERTY(EditAnywhere) bool buildTrivialConnections = true;

    
    UPROPERTY(EditAnywhere) bool drawDistanceField = false;
    UPROPERTY(EditAnywhere) bool drawHitSamples = false;
    UPROPERTY(EditAnywhere) bool drawNearest = false;
	UPROPERTY(EditAnywhere) bool drawSurface = false;
	UPROPERTY(EditAnywhere) bool drawHalfEdge = false;
	UPROPERTY(EditAnywhere) bool drawBounds = false;

	UPROPERTY(EditAnywhere) bool showMatchingExamples = false;
	UPROPERTY(EditAnywhere) FBox showMatchingExamplesBounds;

	UPROPERTY(EditAnywhere) UStaticMesh * debug_matchingPointMesh = nullptr;
	UPROPERTY(EditAnywhere) UMaterialInterface * debug_matchingPointMaterial = nullptr;
	UPROPERTY(EditAnywhere) float debug_matchingPointSize = 0.1f;

	UPROPERTY( EditAnywhere ) bool disableOptimization = false;
	UPROPERTY( EditAnywhere ) int optimizationRounds = 1;
	UPROPERTY( EditAnywhere ) bool useGlobalOptimization = false;

    UPROPERTY() TArray<UStaticMeshComponent*> singularityDebugSpheres;
    UPROPERTY() TArray<UStaticMeshComponent*> seedDebugSpheres;
    
    UPROPERTY(EditAnywhere) bool debugHorizon = false;
    UPROPERTY(EditAnywhere) UMaterialInstance * debugHorizonMaterial = nullptr;

	UPROPERTY( EditAnywhere ) UMaterialInterface * directionFieldMaterial = nullptr;
	UPROPERTY( EditAnywhere ) UMaterialInterface * debugSurfaceMaterial = nullptr;

	UPROPERTY( EditAnywhere ) FBox trimBounds;

    UPROPERTY(VisibleAnywhere, BlueprintReadOnly) int32 statElementCount = 0;
    
    UStaticMeshComponent * staticMeshComponent();
	URuntimeMeshComponent * debugMeshComponent(); // for debug drawing
    
    URuntimeMeshComponent * meshInterfaceRuntimeMesh = nullptr;


    
	// Sets default values for this component's properties
	URegionGrowingComponent();

	// Called when the game starts
	virtual void InitializeComponent() override;
    
	virtual void BeginDestroy() override;
	
	// Called every frame
	virtual void TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction ) override;

	void _tickVisibility();

	virtual void SetActive( bool newActive, bool reset = false)
	{
		Super::SetActive( newActive, reset );

		if(newActive == false && reset == true)
		{
			LoadExemplar();
		}
	}



    UFUNCTION(BlueprintCallable, Category=Generation) void LoadExemplar();
    UFUNCTION(BlueprintCallable, Category=Generation) void Generate();
    UFUNCTION(BlueprintCallable, Category=Generation) void GlobalOptimization();
	UFUNCTION(BlueprintCallable, Category=Generation) void ClearOutput();
	UFUNCTION(BlueprintCallable, Category=Generation) void ClearElementsKeepPaths();

	UFUNCTION( BlueprintCallable, Category = Generation ) void OutputStats();
	UFUNCTION( BlueprintCallable, Category = Generation ) void OutputSpaceSeparatedStats();
	UFUNCTION( BlueprintCallable, Category = Generation ) void OutputParameters();

	UFUNCTION( BlueprintCallable, Category = Generation ) void Trim();




#if WITH_EDITOR
    void PostEditChangeProperty(struct FPropertyChangedEvent& e) override;
#endif    

    void addSingularity(const FVector& vertex, float k);
    void removeSingularityAt(unsigned int index);
    
    void addSeed(const FVector& localPoint);
    void removeSeedAt(unsigned int index);
    
    bool isMeshInterfaceReady() { return _didInit; };
    std::shared_ptr<tcodsMeshInterfaceBase> meshInterface() { return _meshInterface; };
	void updateMeshInterface();


    class ASimulationSnapshotActor * createGraphicalSnapshotActor(UWorld * targetWorld);
    
	OutputDomainExport exportOutputDomain();
	void loadElementDomain( FGraphSnapshot& snapshot );

	void setElementVisibility(bool visible);

    // ----------------------------------------------------------
    // Stretching
    // ----------------------------------------------------------
    void startStretch( FVector point );
    void updateStretch( FVector point );
    void endStretch( FVector point );
    
    // ----------------------------------------------------------
    // Exemplar Painting
    // ----------------------------------------------------------
    void startPaint();
    void addBrushPoint( FVector point, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface );
	void addBrushPointWithRadius( FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface );
    void endPaint();

	void startErase();
	void eraseAt( FVector point, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface );
	void eraseInRadiusAt( FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface);
	void endErase();

	Algorithm::ExampleSelectionPtr getExampleSelection();
	void updateExampleSelection( std::vector<AElementActor*> newSelection, float weight );

	void generateExemplar();

	FBox generationBoundsWorld();


	FVector toExemplar( FVector worldSpace );
	FVector toWorld( FVector exemplarSpace );

	InstanceManager * getInstanceManager() { return &_instanceManager; }
    
	// Subclasser's API
public:
	virtual Algorithm* getOrCreateAlgorithm();

	virtual void loadParameters();
	virtual FString parametersString();

	virtual void didClear() {}

protected:
	virtual void initAlgorithms();

public:
	struct ElementTypeDescription
	{
		ElementTypeDescription( AElementActor * fromActor )
		{
			UStaticMeshComponent * staticMeshComponent = fromActor->GetStaticMeshComponent();

			if(!staticMeshComponent)
				return;

			mesh = staticMeshComponent->GetStaticMesh();

			for(int i = 0; i < staticMeshComponent->GetNumMaterials(); ++i)
				materials.Add( staticMeshComponent->GetMaterial( i ) );

			if(staticMeshComponent->GetNumMaterials())
				material = staticMeshComponent->GetMaterial( 0 );
		}

		ElementTypeDescription( AActor * fromActor )
		{
			UStaticMeshComponent * staticMeshComponent = fromActor->FindComponentByClass<UStaticMeshComponent>();

			if(!staticMeshComponent)
				return;

			mesh = staticMeshComponent->GetStaticMesh();

			for(int i = 0; i < staticMeshComponent->GetNumMaterials(); ++i)
				materials.Add( staticMeshComponent->GetMaterial( i ) );

			if( staticMeshComponent->GetNumMaterials() )
				material = staticMeshComponent->GetMaterial( 0 );
		}

		ElementTypeDescription( const ElementTypeDescription& other )
		{
			mesh = other.mesh;
			material = other.material;
			materials = other.materials;
		}

		ElementTypeDescription()
		{
			mesh = nullptr;
			material = nullptr;
		}

		UStaticMesh * mesh = nullptr;
		UMaterialInterface * material = nullptr;

		TArray<UMaterialInterface*> materials;
	};

	struct InferredElementTypeDescriptionComparator
	{
		bool operator()(const URegionGrowingComponent::ElementTypeDescription& a, const URegionGrowingComponent::ElementTypeDescription& b) const
		{
			if (a.mesh != b.mesh)
				return a.mesh < b.mesh;

			if (a.material != b.material)
				return a.material < b.material;

			return false;
		}
	};

protected:
    void _init();

	void _initContextAndGraph();
    void _initRuntimeMesh();

    virtual void _loadResult(const AlgorithmResult& result);

	virtual void _mainThreadGenerationWorkDone() {};

	void _privateClear();

	// clears the output entity/element to instanced static mesh mappings
	void _clearInstancesData();

    void _updateInstance(FGraphNodeHandle element);

	auto _computeInstanceTransform( FGraphNode& element )->FTransform;

    void _destroyChildren(AActor* actor);
    
    void _updateSeeds();
    void _drawSeeds();
    void _destroySeeds();
    
	void _updateGenerationLimits();
	void _updateDrawBounds();
	void _drawBounds();
	void _undrawBounds();

	void _updateDrawHalfEdge();
	void _drawHalfEdge();
	void _undrawHalfEdge();

	void _updateShowMatching();

	void _updateDrawDirectionField();
    void _drawDirectionField();
    void _undrawDirectionField();
    
    bool _initMeshInterface();
    void _updateDirectionField();
    
    auto _getGenerationLimits() -> Algorithm::AABB;
    
    // Called during LoadExemplar
    void _loadInstancedStaticMeshes();
    
    void _loadSingleSampleExemplar(Domain& exemplarDomain);

	void _loadTypeDescriptionsFromNewElementsAndSetTypes( const std::vector<FGraphNodeHandle>& elements, FGraph& graph );
	void _reassignFaceIndices(std::vector<FGraphNodeHandle>& elements, FGraph& graph);

	// Generation
    auto _spaceFillingStartingPositions(Algorithm::AABB& limits, FVector& cameraLocation) -> std::vector<PositionFace>;
    auto _surfaceStartingPositions(Algorithm::AABB& limits, FVector& cameraLocation, FVector& cameraDirection) -> std::vector<PositionFace>;
    
   
	void _removeElements( std::vector<FGraphNodeHandle> toRemove );

	FString _optimazationParametersString( FNeighbourhoodParameters& params );

	std::vector<FGraphNodeHandle> _actorsToElements( std::vector<AElementActor*>& actors );

    
protected:
	double _duration = 0.0f;

    int tickCount;

	FGraph _synthesisGraph;
	std::unique_ptr<SynthesisContext> _synthesisContext;

	std::unique_ptr<Algorithm> algorithm = nullptr;
    


	// Example selection
	std::unordered_map<AElementActor*, FGraphNodeHandle> _exampleActorToElement;
    std::map<int16, ElementTypeDescription> _elementTypesToDescriptions;
   
    std::map<int16, AActor*> _entityIDsToActors; // if an entity has a mesh, this table looks up the actor with that mesh
    std::map<int16, UInstancedStaticMeshComponent*> _entityIDsToInstancedMeshes;
    std::map<int16, int32> _entityIndexForID; // for each actor assigned to an entity id, we also assign the element samples take from that actor
    std::map<Entity*, int32> _entitiesToInstanceIndices;
    std::map<int16, FVector> _entityIDsToMeshOffsets;
    
    

    
    struct CachedElement
    {
        FTransform transform;
		int16 type;
        int32 instance;
    };

	

	InstanceManager _instanceManager;
    
    std::unordered_map<FGraphNodeHandle, CachedElement> _cache;
    ElementKDTree _elementKDTree = ElementKDTree(_synthesisGraph);   // because algorithm is running on a background thread, we use this as a cache
                                    // of the last state of the element positions
    
    std::vector<std::vector<FGraphNodeHandle>> _previouslyCulled;
    
    UStaticMeshComponent * _nearestPointDebugMesh = nullptr;
    UStaticMeshComponent* nearestPointDebugMesh();
    
    UInstancedStaticMeshComponent * _seedPointsInstancedMeshComponent = nullptr;
    UInstancedStaticMeshComponent * seedPointsInstancedMeshComponent();
    
	UInstancedStaticMeshComponent * _debug_matchingISMC = nullptr;
	UInstancedStaticMeshComponent& debug_matchingISMC();

    std::shared_ptr<tcodsMeshInterfaceBase> _meshInterface = nullptr;
    bool _didInit = false;
    
    std::atomic<bool> _readyToSubmitNextGenerateCall = {true};
    ctpl::thread_pool _generationWorker;
	ctpl::detail::Queue<std::function<void()>> _mainThreadWork;

    std::function<void()> _generationWorkDone;
    
    double generationTime;
    
 
    FVector _exemplarPosition = FVector(0.0f, 0.0f, 0.0f);
    
    std::shared_ptr<RGOcclusionTester> _occlusionTester = nullptr;
    
	int _resetGuard = 0;

    // -----------------
    // Stretching
    // -----------------
    FVector _stretchStart;
    FVector _stretchLast;

	FVector * _exemplarPoint;

    struct StretchElement
    {
        size_t index;
		size_t section;
        FVector originalPosition;
        float scale;
    };

    std::vector<StretchElement> _stretchVertices; // tcods::Vertex index and start position for a given vertex
    URuntimeMeshComponent * _debugMesh = nullptr;

	// Painting
	std::vector<PositionFace> _accumulatedStartingPositions;

	std::vector< PositionRadiusFace > _toRemove;



	UPROPERTY(EditAnywhere, BlueprintReadOnly)
	UBoxComponent * _boundsBox = nullptr;
};
