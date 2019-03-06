// Fill out your copyright notice in the Description page of Project Settings.

#include "LifeBrush.h"
#include "RegionGrowingComponent.h"
#include "ElementActor.h"
#include "ContextActor.h"
#include "Utility.h"
#include "tcodsMeshInterface.h"
#include "VolumeComponent.h"

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "ShipEditorSimulation/GraphSnapshot.h"

#include <vector>
#include <map>
#include <cmath>
#include <unordered_map>
#include <ctime>

#include "Algorithm/Domain.h"
#include "Algorithm/Entity.hpp"

// unreal includes
#include "StaticMeshResources.h"
#include "RuntimeMeshComponent.h"
#include "RuntimeMeshLibrary.h"
#include "SimulationSnapshotActor.h"
//#include "InstancedStaticMeshComponent.h"


#include "Algorithm/tcods/MeshIO.h"
#include "Algorithm/tcods/Problem.h"

#if WITH_EDITOR
#include "LevelEditor.h"
#endif

// fuck windows http://stackoverflow.com/questions/118774/is-there-a-clean-way-to-prevent-windows-h-from-creating-a-near-far-macro
#undef near

static double _pi() { return std::acos(-1.0); }

// we don't do `enum class MeshSections : int32` because we want this to be an int, so we embed it in a struct for the namespace
struct MeshSections  
{
	enum Type : int32
	{
		DirectionField = 0,
		Samples,
		HalfEdge,
		Bounds
	};
};


// Notes on Unreal's Signed Distance Fields
// ----------------------------------------
//
// - FDistanceFieldVolumeData
// - void FDistanceFieldAsyncQueue::Build(FAsyncDistanceFieldTask* Task, FQueuedThreadPool& ThreadPool)
// - MeshUtilities::GenerateSignedDistanceFieldVolumeData
// - FMeshDistanceFieldAsyncTask::DoWork
//
// All the magic happens in MeshUtilities.cpp
//
// Essentially, Unreal uses a Monte Carlo sampling method to determine the distance
// to the nearest surface inside of the mesh. They accelerate the Monte Carlo sampling
// with a tree data structure.
//
// The volume data is floating point. So this is the floating point transform that we were looking for.


bool RGOcclusionTester::isVisible( PositionFace facePoint, float radius )
{
	const FVector& point = unreal(facePoint.position);

	float distanceSqrd = FVector::DistSquared( cameraLocation, point );

	if(!reader->inFrustum( point ))
		return false;

	FVector dir = (point - cameraLocation).GetSafeNormal();

	if(facePoint.surfaceIndex.isOnSurface())
	{
		auto face = _meshInterface->mesh(facePoint.surfaceIndex.sectionIndex).face( facePoint.surfaceIndex.faceIndex );

		FVector normal = unreal(face->normal());

		if(FVector::DotProduct( normal, dir ) <= 0)
			return false;
	}


	auto hit = _meshInterface->getIntersectionAndFace( cameraLocation, dir );

	if(!hit.first)
		return false;

	float distanceToHitSqrd = FVector::DistSquared( point, hit.second.point );

	if(distanceToHitSqrd < radius * radius)
		return true;
	else
		return false;

	//float pointDepth = reader->depthOf( point );
	//float depth = reader->depthAt( point );
 //   
	//bool visible = pointDepth < depth;   

	//// dumb, but using for fast breakpoints
	//if(visible)
	//	return true;
	//else
	//	return false;
}






// Sets default values for this component's properties
URegionGrowingComponent::URegionGrowingComponent()
{
    // Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
    // off to improve performance if you don't need them.
    bWantsInitializeComponent = true;
    PrimaryComponentTick.bCanEverTick = true;
    bIsActive = true;
    bTickInEditor = true;
}

void URegionGrowingComponent::BeginDestroy()
{
    Super::BeginDestroy();
}

// Called when the game starts
void URegionGrowingComponent::InitializeComponent()
{
    Super::InitializeComponent();
    

}

void URegionGrowingComponent::setElementVisibility( bool visible )
{
	_instanceManager.setVisibility( visible );
}

// -----------------------------------------------------------------------------
// Sketch
// -----------------------------------------------------------------------------

void URegionGrowingComponent::_removeElements(std::vector<FGraphNodeHandle> toRemove)
{
	if(toRemove.size() == 0)
		return;

	algorithm->removeOutputElements( toRemove );

	// update instances
	AlgorithmResult result;

	auto& output = algorithm->output();

	std::vector<FGraphNodeHandle>& remaining = result.generated;
	for(auto& e : output)
		remaining.push_back( e.nodeHandle() );

	_entityIDsToActors.clear();
	_elementKDTree.clear();

	_instanceManager.clear();

	_loadResult( result );
}

void URegionGrowingComponent::startStretch( FVector point )
{
    _stretchStart = point;
    _stretchLast = point;

    _stretchVertices.clear();


    const float stretchRadiusSqrd = stretchRadius * stretchRadius;

	for (auto& pair : _meshInterface->_sections)
	{
		auto& mesh = *pair.second.get();
		size_t section = pair.first;

		for (auto& vertex : mesh.vertices)
		{
			FVector v = unreal(vertex.position);

			float ds = FVector::DistSquared(v, point);

			if (ds < stretchRadiusSqrd)
			{
				float scale = 1.0f - sqrt(ds) / stretchRadius;

				_stretchVertices.emplace_back(StretchElement{ size_t(vertex.index), section, v, scale });
			}
		}
	}



    pauseSynthesis = true;

	// wait for the background thread to finish
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	// find any nearby elements
	auto& output = algorithm->output();

	auto toRemove = output.nearestInRadius( eigen( point ), stretchRadius );

	_removeElements( toRemove );
}

void URegionGrowingComponent::updateStretch( FVector point )
{
    UE_LOG( LogTemp, Warning, TEXT( "updateStretch %s"), *point.ToString() ); 

    auto& lookup = _meshInterface->_halfEdgeVertex_to_runtimeMeshVertices;

    float stretchRadiusSqrd = stretchRadius * stretchRadius;

    FVector worldDelta = point - _stretchLast;

    FTransform transform = meshInterfaceRuntimeMesh->GetComponentTransform();
    FVector delta = transform.InverseTransformVector( worldDelta );

	for (auto& pair : _meshInterface->_sections)
	{
		auto& mesh = *pair.second.get();
		size_t section = pair.first;

		auto meshSection = meshInterfaceRuntimeMesh->BeginSectionUpdate(section);

		for (auto& stretchElement : _stretchVertices)
		{
			// mapping from half-edge vertices to unreal mesh vertices (unreal duplicates)
			for (auto i : lookup[stretchElement.index])
			{
				meshSection->SetPosition(i, meshSection->GetPosition(i) + delta * stretchElement.scale);
			}

			// update the tcods::mesh
			mesh.vertices[stretchElement.index].position += to_tcods(worldDelta * stretchElement.scale);
		}

		meshSection->Commit(true, false, false, false, false);
	}

    _stretchLast = point;
}

void URegionGrowingComponent::endStretch( FVector point )
{
    // wait for the background thread to finish
    _generationWorker.push( []( int threadID ) mutable {} ).get();

	_meshInterface->rebuildBvh();

    // find any nearby elements
    auto& output = algorithm->output();

    auto nearby = output.nearestInRadius( eigen( point ), stretchRadius * 1.2 );

    algorithm->addToHorizon( nearby );

    pauseSynthesis = false;
}

// -----------------------------------------------------------------------------
// Painting
// -----------------------------------------------------------------------------

void URegionGrowingComponent::startPaint()
{
	// hack for now
	// we should do something more intelligent, like checking if the brush stroke has completed growing
	_generationWorker.push( [this]( int threadID ) mutable
	{
		algorithm->clearBrushPoints();
	} );

	_toRemove.clear();
}

void URegionGrowingComponent::addBrushPoint( FVector point, FSurfaceIndex surfaceIndex /* = FSurfaceIndex::OffSurface */ )
{
	addBrushPointWithRadius( point, brushSize, surfaceIndex );
}

void URegionGrowingComponent::addBrushPointWithRadius( FVector point, float radius, FSurfaceIndex surfaceIndex /* = FSurfaceIndex::OffSurface */ )
{
	auto eigen_point = eigen( point );

	_generationWorker.push( [this, eigen_point, radius, surfaceIndex]( int threadID ) mutable
	{
		PositionRadiusFace brushPoint( eigen_point, radius, surfaceIndex);
		algorithm->addBrushPoint( brushPoint );
		algorithm->setCurrentBrushPoint( brushPoint );
	} );

	_accumulatedStartingPositions.emplace_back( eigen( point ), surfaceIndex);
}

void URegionGrowingComponent::endPaint()
{
	_generationWorker.push( [this]( int threadID ) mutable
	{
		algorithm->clearCurrentBrushPoint();
	} );
}

void URegionGrowingComponent::startErase()
{
	_generationWorker.push([this](int threadID) mutable
	{
		algorithm->clearBrushPoints();
	});

	_toRemove.clear();
}

void URegionGrowingComponent::eraseAt( FVector point, FSurfaceIndex surfaceIndex )
{
	if (_toRemove.size())
	{
		auto& last = _toRemove.back();
		
		const float distSqrd = (eigen(point) - last.position).squaredNorm();
		if (distSqrd < last.radius * 0.2)
			return;
	}


	_toRemove.emplace_back( eigen(point), brushSize, surfaceIndex);
}

void URegionGrowingComponent::eraseInRadiusAt( FVector point, float radius, FSurfaceIndex surfaceIndex /* = FSurfaceIndex::OffSurface */ )
{
	_toRemove.emplace_back( eigen( point ), radius, surfaceIndex);
}


void URegionGrowingComponent::endErase()
{

}

std::vector<FGraphNodeHandle> URegionGrowingComponent::_actorsToElements( std::vector<AElementActor*>& actors )
{
	std::vector<FGraphNodeHandle> elements;

	for(AElementActor * actor : actors)
	{
		if(_exampleActorToElement.find( actor ) == _exampleActorToElement.end())
			continue;

		FGraphNodeHandle element = _exampleActorToElement[actor];

		if(!element)
			continue;

		elements.push_back( element );
	}

	return elements;
}

Algorithm::ExampleSelectionPtr URegionGrowingComponent::getExampleSelection( )
{
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	return algorithm->getExampleSelection();
}

void URegionGrowingComponent::updateExampleSelection(std::vector<AElementActor*> newSelection, float weight )
{
	std::vector<FGraphNodeHandle> elements = _actorsToElements( newSelection );

	_generationWorker.push( [this, elements, weight]( int threadID ) mutable {
		algorithm->updateExampleSelection( elements, weight );
	} );
}



void URegionGrowingComponent::generateExemplar()
{
	if(!targetExemplar)
		return;

	auto& output = algorithm->output();

	output.graph.view<FElementObject>().each([&](FGraphNode& elementNode, FElementObject& element)
	{
		UInstancedStaticMeshComponent * insancedMesh = _instanceManager.meshForType(element.type);


		AElementActor * newActor = GetWorld()->SpawnActor<AElementActor>();

		UStaticMeshComponent * newActorMesh = newActor->GetStaticMeshComponent();

		newActorMesh->SetStaticMesh(insancedMesh->GetStaticMesh());
		newActorMesh->SetMaterial(0, insancedMesh->GetMaterial(0));
		newActorMesh->SetMobility(EComponentMobility::Movable);
		newActorMesh->SetCollisionProfileName(TEXT("NoCollision"));

		newActor->RegisterAllComponents();
		newActor->SetMobility(EComponentMobility::Movable);

		newActor->AttachToActor(targetExemplar, FAttachmentTransformRules::KeepRelativeTransform);

		newActor->SetActorLocationAndRotation(elementNode.position, elementNode.orientation);
		newActor->SetActorScale3D(FVector(1.0f));

		ElementTuple elementTuple(elementNode, element);

		newActor->readFromElement(elementTuple, output.graph);
	});

	for(auto& element : output)
	{

	}
}

FVector URegionGrowingComponent::toExemplar( FVector worldSpace )
{
	FVector exemplarSpace = worldSpace;
	exemplarSpace.Z = worldSpace.Z - exemplar->GetActorLocation().Z;

	return exemplarSpace;
}

FVector URegionGrowingComponent::toWorld( FVector exemplarSpace )
{
	FVector worldSpace = exemplarSpace;
	worldSpace.Z = exemplarSpace.Z + exemplar->GetActorLocation().Z;

	return worldSpace;
}

void URegionGrowingComponent::initAlgorithms()
{
	if(_synthesisContext->threadPool == nullptr)
	{
		//unsigned int n = std::thread::hardware_concurrency();

		_synthesisContext->threadPool = std::make_shared<ctpl::thread_pool>();
		_synthesisContext->threadPool->resize( nThreads );
	}

	getOrCreateAlgorithm();
}

void URegionGrowingComponent::_init()
{
	_initContextAndGraph();

	initAlgorithms();

    _initRuntimeMesh();
    _initMeshInterface();
    _updateDirectionField();
    
    _updateDrawDirectionField();
    _updateSeeds();
	_updateDrawHalfEdge();
	_updateDrawBounds();
    
    _generationWorker.resize(1);

    
#ifdef EIGEN_VECTORIZE
    UE_LOG(LogTemp, Warning, TEXT("Eigen is vectorized"));
#endif

	_didInit = true;
}


void URegionGrowingComponent::_initContextAndGraph()
{
	if (_synthesisContext)
		return;

	_synthesisGraph.init();

	_synthesisContext = std::make_unique<SynthesisContext>(_synthesisGraph);
}

// Called every frame
void URegionGrowingComponent::TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction )
{
    Super::TickComponent( DeltaTime, TickType, ThisTickFunction );

	_instanceManager.rotateBuffers();

#if WITH_EDITOR
	_updateGenerationLimits();
#endif
    
    if( !_didInit )
    {
        _init();
        LoadExemplar();
    }

	if(stopDuration > 0.0f && _duration > stopDuration)
	{
		pauseSynthesis = true;
	}

    if( synthesizeDistance > 0 && !pauseSynthesis )
        Generate();



	_tickVisibility();
}

void URegionGrowingComponent::_tickVisibility()
{
	if(meshInterfaceRuntimeMesh)
		meshInterfaceRuntimeMesh->SetVisibility( drawSurface );

	if(drawNearest)
	{
		auto world = GetWorld();
		if(world == nullptr)
			return;

		auto viewLocations = world->ViewLocationsRenderedLastFrame;
		if(viewLocations.Num() == 0)
			return;

		FVector location = viewLocations[0];

		UStaticMeshComponent * mesh = staticMeshComponent();
		if(mesh == nullptr)
			return;

		UStaticMeshComponent * point = nearestPointDebugMesh();

		if(point == nullptr)
			return;

		auto nearest = _meshInterface->nearestPointOnMesh( location );

		point->SetVisibility( true );
		point->SetWorldLocation( nearest.point );
	}

	if(!drawNearest)
	{
		UStaticMeshComponent * point = nearestPointDebugMesh();

		if(point != nullptr)
			point->SetVisibility( false );
	}

	auto world = GetWorld();
	if(world == nullptr)
		return;

	auto viewLocations = world->ViewLocationsRenderedLastFrame;
	if(viewLocations.Num() == 0)
		return;

	FVector location = viewLocations[0];

	auto nearestOnMesh = _meshInterface->nearestPointOnMesh(location);

	bool hasNearest = nearestOnMesh.surfaceIndex.isOnSurface();

	float nearestElementDistance = (nearestOnMesh.point - location).Size();


	// hide elements that are too near to the camera for cut-away views
	if(innerQuellingDistance > 0.0f)
	{
		std::set<UInstancedStaticMeshComponent*> dirty;

		// restore the previously culled
		for(int i = 0; i < _previouslyCulled.size(); i++)
		{
			auto& elements = _previouslyCulled[i];

			if(elements.size() == 0)
				continue;

			auto mesh = _instanceManager.meshForType( i );

			if (!mesh)
				continue;

			for(FGraphNodeHandle e : elements)
			{
				auto& cached = _cache[e];

				mesh->UpdateInstanceTransform( cached.instance, cached.transform, true );
			}

			dirty.insert( mesh );
		}

		// cull some new guys
		auto near = _elementKDTree.nearestInRadius( eigen( location ), innerQuellingDistance );

		//std::vector<Element*> near;
		//auto& allElements = algorithm->output();

		//float innerSqrd = innerQuellingDistance * innerQuellingDistance;

		//float outerSqrd = outerQuellingDistance == 0.0f ? std::numeric_limits<float>::max() : outerQuellingDistance * outerQuellingDistance;

		//auto eigenLocation = eigen( location );
		//for(auto& e : allElements)
		//{
		//	float d = (e->position - eigenLocation).squaredNorm();

		//	if(d < innerSqrd || d > outerSqrd)
		//		near.push_back( e.get() );
		//}

		int maxType = 0;
		Element::eachElement(near, algorithm->output().graph, [&](FElementObject& e)
		{
			if (e.type > maxType)
				maxType = e.type;
		});

		std::vector<std::vector<FGraphNodeHandle>> categorized( maxType + 1 );

		Element::eachElement(near, algorithm->output().graph, [&](FElementObject& e)
		{
			categorized[e.type].push_back(e.nodeHandle());
		});

		_previouslyCulled = categorized;

		for(int i = 0; i < categorized.size(); i++)
		{
			auto& elements = categorized[i];

			if(elements.size() == 0)
				continue;

			auto mesh = _instanceManager.meshForType( i );

			if (mesh)
			{
				for (FGraphNodeHandle e : elements)
				{
					auto& cached = _cache[e];

					FTransform t;
					mesh->GetInstanceTransform(cached.instance, t, true);
					t.SetTranslation(FVector(100000.0f, 0.0f, 0.0f));

					mesh->UpdateInstanceTransform(cached.instance, t, true);
				}

				dirty.insert(mesh);
			}
		}

		for(auto mesh : dirty)
			mesh->MarkRenderStateDirty();

		// hide the mesh?
		if(hideSurfaceDistance > 0.0f)
		{
			auto mesh = staticMeshComponent();
			if(!mesh)
				return;

			if(hasNearest && nearestElementDistance < hideSurfaceDistance * hideSurfaceDistance)
			{

				if(mesh->IsVisible())
					mesh->SetVisibility( false );

				if(meshInterfaceRuntimeMesh && meshInterfaceRuntimeMesh->IsVisible())
					meshInterfaceRuntimeMesh->SetVisibility( false );
			}
			else
			{
				if(!mesh->IsVisible())
					mesh->SetVisibility( true );

				if(meshInterfaceRuntimeMesh && !meshInterfaceRuntimeMesh->IsVisible())
					meshInterfaceRuntimeMesh->SetVisibility( true );
			}
		}
	}

	// reset if we're faraway
	if(_elementKDTree.size() && resetDistance >= 0.0001f && hasNearest && nearestElementDistance > resetDistance * resetDistance)
	{
		ClearOutput();
	}
}





auto URegionGrowingComponent::_getGenerationLimits() -> Algorithm::AABB
{
    Algorithm::AABB limits;
    limits.aabb = Eigen::AlignedBox3f(eigen(generationLimitsMin), eigen(generationLimitsMax));

    return limits;
}

void URegionGrowingComponent::_privateClear()
{
	_accumulatedStartingPositions.clear();

	// rebuild the nearest in limits index in case the limits have changed
	Algorithm::AABB limits = _getGenerationLimits();

	_elementKDTree.clear();

	generationTime = 0.0f;


	_instanceManager.clear();

	for(auto instancedComponent : _entityIDsToInstancedMeshes)
		instancedComponent.second->ClearInstances();

	_clearInstancesData();

	// clear the previous results
	_destroyChildren( output() );

	loadParameters();

	_generationWorkDone = nullptr;
	_readyToSubmitNextGenerateCall = true;

	_duration = 0.0f;

	didClear();
}

void URegionGrowingComponent::ClearOutput()
{
	if( !_didInit )
		LoadExemplar();

	// wait for work to be finished
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	_privateClear();

	algorithm->clearPainting();
	algorithm->clearOutput();
}

void URegionGrowingComponent::ClearElementsKeepPaths()
{
	if(!_didInit)
		LoadExemplar();

	// wait for work to be finished
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	_privateClear();

	algorithm->clearOutput();
}

void URegionGrowingComponent::OutputStats()
{
	FString stats = algorithm->roundSummary().statsString();

	UE_LOG( LogTemp, Warning, TEXT("%s"), *stats );
}

void URegionGrowingComponent::OutputSpaceSeparatedStats()
{
	FString stats = algorithm->roundSummary().separatedValueStatsString(" ");

	UE_LOG( LogTemp, Warning, TEXT( "%s" ), *stats );
}

FString URegionGrowingComponent::_optimazationParametersString( FNeighbourhoodParameters& params )
{
	FString result = "";

	result += FString::Printf( TEXT( "   radius %f\n" ), params.radius );
	result += FString::Printf( TEXT( "   kNearest %f\n" ), params.kNearest );

	return result;
}


FString URegionGrowingComponent::parametersString()
{
	FString params = "";

	{
		const UEnum * AlgorithmEnum = FindObject<UEnum>( ANY_PACKAGE, TEXT( "EAlgorithm" ) );
		int32 int32_algorithm = (int32)algorithmType;

		params += FString::Printf( TEXT( "algorithmType: %s\n" ), *(AlgorithmEnum->GetDisplayNameTextByValue( int32_algorithm ).ToString()) );
	}

	params += FString::Printf( TEXT( "stopDuration %f\n" ), stopDuration );

	params += FString::Printf( TEXT( "perElementParameters %d\n" ), perElementParameters );

	params += FString::Printf( TEXT( "disableOptimization %d\n" ), disableOptimization );
	params += FString::Printf( TEXT( "optimizationRounds %d\n" ), optimizationRounds );
	params += FString::Printf( TEXT( "useGlobalOptimization %d\n" ), useGlobalOptimization );

	params += FString::Printf( TEXT( "enableRoundSummaries %d\n" ), enableRoundSummaries );
	params += FString::Printf( TEXT( "calculate_kCoherenceEnergy %d\n" ), calculate_kCoherenceEnergy );
	params += FString::Printf( TEXT( "calculate_bruteForceEnergy %d\n" ), calculate_bruteForceEnergy );

	{
		const UEnum * GenerationMode = FindObject<UEnum>( ANY_PACKAGE, TEXT( "EGenerationMode" ) );
		int32 int32_generationMode = (int32)generationMode;

		params += FString::Printf( TEXT( "generationMode: %s\n" ), *(GenerationMode->GetDisplayNameTextByValue( int32_generationMode ).ToString()) );
	}

	params += FString::Printf( TEXT( "sketchBasedForceTangentPlane %d\n" ), sketchBasedForceTangentPlane );

	params += FString::Printf( TEXT( "generationParameters\n" ) );
	params += _optimazationParametersString( generationParameters );

	params += FString::Printf( TEXT( "optimizationParameters\n" ) );
	params += _optimazationParametersString( optimizationParameters );

	params += FString::Printf( TEXT( "useTypeVoting %d\n" ), useTypeVoting );
	params += FString::Printf( TEXT( "removeOverlaps %d\n" ), removeOverlaps );
	params += FString::Printf( TEXT( "enableReassignment %d\n"), enableReassignment);
	params += FString::Printf( TEXT( "generationInnerRadius %f\n" ), generationRadius );
	params += FString::Printf( TEXT( "kCoherence %d\n" ), kCoherence );
	params += FString::Printf( TEXT( "minAssignmentDistance %f\n" ), minAssignmentDistance );
	params += FString::Printf( TEXT( "freespaceRadius %f\n" ), freespaceRadius );
	params += FString::Printf( TEXT( "typeCost %f\n" ), typeCost );
	params += FString::Printf( TEXT( "gradientTerm %f\n" ), gradientTerm );
	params += FString::Printf( TEXT( "relaxation %f\n" ), relaxation );
	params += FString::Printf( TEXT( "flipSurfaceNormals %d\n" ), flipSurfaceNormals );
	params += FString::Printf( TEXT( "ignoreBadSuggestions %d\n" ), ignoreBadSuggestions );
	params += FString::Printf( TEXT( "ignoreBadSuggestionsDistanceFactor %f\n" ), ignoreBadSuggestionsDistanceFactor );

	params += FString::Printf( TEXT( "kCoherenceClusteringPenalty %f\n" ), kCoherenceClusteringPenalty );
	params += FString::Printf( TEXT( "kCoherenceClusteringPenaltyRadius %f\n" ), kCoherenceClusteringPenaltyRadius );

	params += FString::Printf( TEXT( "nThreads %d\n" ), nThreads );

	params += FString::Printf( TEXT( "forceRotation %d\n" ), forceRotation );
	params += FString::Printf( TEXT( "forcedRotation %s\n" ), *forcedRotation.ToString() );

	params += FString::Printf( TEXT( "seedsIgnoreMeshInBoundary %d\n" ), seedsIgnoreMeshInBoundary );

	params += FString::Printf( TEXT( "sourceHistogramRadius %f\n" ), sourceHistogramRadius );
	params += FString::Printf( TEXT( "sourceHistogramWeight %f\n" ), sourceHistogramWeight );
	params += FString::Printf( TEXT( "samplingDistanceWeight %f\n" ), samplingDistanceWeight );

	params += FString::Printf( TEXT( "patchInitializationExemplarDivisions %s\n" ), *patchInitializationExemplarDivisions.ToString() );
	params += FString::Printf( TEXT( "patchMarginPercentage %f\n" ), patchMarginPercentage );

	params += FString::Printf( TEXT( "brushSize %f\n" ), brushSize );

	return params;
}

void URegionGrowingComponent::OutputParameters()
{
	FString params = parametersString();

	UE_LOG( LogTemp, Warning, TEXT( "Parameters\n------------\n%s" ), *params );
}

void URegionGrowingComponent::Trim()
{
	// wait for work to be finished
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	std::vector<FGraphNodeHandle> toRemove;

	for(auto& e : algorithm->output())
	{
		FGraphNode& node = algorithm->output().graph.node(e.nodeHandle());
		FVector p = node.position;

		if( !trimBounds.IsInside( p ) )
			toRemove.push_back( e.nodeHandle() );
	}

	this->_removeElements( toRemove );
}

Algorithm* URegionGrowingComponent::getOrCreateAlgorithm()
{
	Algorithm * result = nullptr;

	if(algorithmType == EAlgorithm::RegionGrowing)
		algorithm = std::make_unique<Algorithm_RegionGrowing>(*_synthesisContext);
	else if(algorithmType == EAlgorithm::PatchCopy)
		algorithm = std::make_unique<Algorithm_PatchCopy>(*_synthesisContext);

	return result;
}

void URegionGrowingComponent::LoadExemplar()
{
    if( !_didInit )
        _init();

	// select algorithm
	getOrCreateAlgorithm();

    // wait for work to be finished
    _generationWorker.push( [](int threadID) mutable {} ).get();

	_privateClear();
    
    algorithm->clear();

    _elementTypesToDescriptions.clear();
    _entityIDsToActors.clear();
    _elementKDTree.clear();
    
    if( !exemplar )
        return;
    
    if( Cast<AStaticMeshActor>(exemplar) != nullptr )
    {
        UE_LOG(LogTemp, Warning, TEXT("URegionGrowingComponent::LoadExemplar - Unreal's SpawnActor does not work with AStaticMeshActor, sorry."));
        return;
    }
    
    
	Domain& domain = algorithm->exemplar();

    _loadSingleSampleExemplar(domain);

	algorithm->loadExemplar();

	// make the default example selection the entire exemplar
	// This is where we should create our default example selection.

    _loadInstancedStaticMeshes();
    
    _generationWorkDone = nullptr;
    _readyToSubmitNextGenerateCall = true;
}

void URegionGrowingComponent::_loadSingleSampleExemplar(Domain& exemplarDomain)
{
	_exampleActorToElement.clear();

    TArray<USceneComponent*> exemplarSceneComponents;
    exemplar->GetRootComponent()->GetChildrenComponents(false, exemplarSceneComponents);
    
  
	std::map<ElementTypeDescription, int, InferredElementTypeDescriptionComparator> actorToTypes;

    // infer our element types
    int nextIndex = 0;
    for( auto sceneComponent : exemplarSceneComponents )
    {
		AElementActor * elementActor = Cast<AElementActor>( sceneComponent->GetOwner() );
        
        if( !elementActor  )
            continue;

		ElementTypeDescription typeDescription( elementActor );
        
        if( actorToTypes.find( typeDescription ) == actorToTypes.end() )
        {
            actorToTypes[typeDescription] = nextIndex;
            _elementTypesToDescriptions[nextIndex] = typeDescription;
            nextIndex++;
        }
    }
   
    // build the elements
    for( auto sceneComponent : exemplarSceneComponents )
    {
		AElementActor * elementActor = Cast<AElementActor>( sceneComponent->GetOwner() );
       
        if( !elementActor )
            continue;
        
        int elementID = actorToTypes[elementActor];
              
		FVector position = toExemplar(sceneComponent->GetComponentLocation());
		FQuat orientation = sceneComponent->GetComponentRotation().Quaternion();
		float scale = sceneComponent->GetComponentScale().X;

		auto h = exemplarDomain.insert(position, orientation, scale);
		ElementTuple element(h, exemplarDomain.graph);

		element.element.type = elementID;
        
        UStaticMeshComponent * mesh = (UStaticMeshComponent*)elementActor->GetComponentByClass(UStaticMeshComponent::StaticClass());
        
        if( mesh == nullptr || mesh->GetStaticMesh() == nullptr )
            continue;
         
        if(elementActor)
			elementActor->writeToElement( element, exemplarDomain.graph );

		_exampleActorToElement[elementActor] = h;
    }
}

void URegionGrowingComponent::_clearInstancesData()
{
	_cache.clear();
	_previouslyCulled.clear();
}

void URegionGrowingComponent::_loadInstancedStaticMeshes()
{
	_instanceManager.clear();
	_instanceManager.clearMeshTypes();

	// deallocate the old instanced meshes
	for(auto pair : _instanceManager.elementTypesToInstancedMeshes)
	{
		UInstancedStaticMeshComponent * instance = pair.second;
		instance->UnregisterComponent();
		instance->DestroyComponent(); // or DestroyComponent?
	}

    // deallocate the old instanced meshes
    for( auto pair : _entityIDsToInstancedMeshes )
    {
        UInstancedStaticMeshComponent * instance = pair.second;
        instance->UnregisterComponent();
        instance->DestroyComponent(); // or DestroyComponent?
    }

    // we need full floating point precision for our instances. disable this garbage
//    UE_LOG(LogTemp, Warning, TEXT("Disabling VET_Half2 support so that instanced static meshes have the required precision"));    
//    GVertexElementTypeSupport.SetSupported(VET_Half2,false); 

    
    _entityIDsToInstancedMeshes.clear();

	_clearInstancesData();

    // Build a table of actors that are the same in terms of mesh and materials.
    // We'll use this table to assign one instanced static mesh component to multiple
    // similar actors.
    std::map<ElementTypeDescription, UInstancedStaticMeshComponent*, InferredElementTypeDescriptionComparator> typeDescriptionsToISMCs;
    
    // allocate the instanced meshes for each element type
    for( const auto& pair : _elementTypesToDescriptions )
    {
        auto elementType = pair.first;

		ElementTypeDescription typeDescription = pair.second;

        auto foundInstance = typeDescriptionsToISMCs.find( typeDescription );
        if( foundInstance == typeDescriptionsToISMCs.end() )
        {         
            UStaticMesh * mesh = typeDescription.mesh;
            if( mesh == nullptr )
                return;
  
            UMaterialInterface * material0 = typeDescription.material;
            if( material0 == nullptr )
                return;
            
            UInstancedStaticMeshComponent * instance = NewObject<UInstancedStaticMeshComponent>(GetOwner());
            
            if( instance == nullptr )
                return;

            // To get the instanced static meshes to save add this
            //https://answers.unrealengine.com/questions/330978/actorcomponent-constructed-from-c-not-saved.html
            // GetOwner()->AddInstanceComponent(instance);
            
            instance->SetStaticMesh(mesh);
            instance->SetMobility(EComponentMobility::Static);
			instance->UseDynamicInstanceBuffer = true;
			instance->KeepInstanceBufferCPUAccess = true;
            instance->SetCollisionProfileName(TEXT("NoCollision"));
			instance->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
			if(material0->GetName() == TEXT( "Invisible" ))
				instance->SetVisibility( false );

			instance->RegisterComponent();

			int i = 0;
			for(UMaterialInterface * material : typeDescription.materials)
			{
				instance->SetMaterial( i, material );
				++i;
			}
            
			_instanceManager.setMeshForType( instance, elementType );

			typeDescriptionsToISMCs[typeDescription] = instance;
        }
        else
        {
            UInstancedStaticMeshComponent * instance = foundInstance->second;
			_instanceManager.setMeshForType( instance, elementType );
		}
    }
    
    typeDescriptionsToISMCs.clear();
    
    // allocate the instanced meshes for each entity type
    for( const auto& pair : _entityIDsToActors )
    {
        auto elementType = pair.first;
        auto actor = pair.second;

		ElementTypeDescription typeDescription( actor );
        
        auto foundInstance = typeDescriptionsToISMCs.find( typeDescription );
        if( foundInstance == typeDescriptionsToISMCs.end() )
        {
            UStaticMeshComponent * meshComponent = (UStaticMeshComponent*)actor->GetComponentByClass(UStaticMeshComponent::StaticClass());
            if( meshComponent == nullptr )
                return;
            
            UStaticMesh * mesh = meshComponent->GetStaticMesh();
            if( mesh == nullptr )
                return;
            
            
            UMaterialInterface * material = meshComponent->GetMaterial(0);
            if( material == nullptr )
                return;
            
            UInstancedStaticMeshComponent * instance = NewObject<UInstancedStaticMeshComponent>(GetOwner());
            
            if( instance == nullptr )
                return;
            
            instance->SetStaticMesh(mesh);
            instance->SetMobility(EComponentMobility::Static);
            instance->SetCollisionProfileName(TEXT("NoCollision"));
			instance->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
            instance->SetMaterial(0, material);
			instance->RegisterComponent();

            _entityIDsToInstancedMeshes[elementType] = instance;
            typeDescriptionsToISMCs[actor] = instance;
        }
        else
        {
            UInstancedStaticMeshComponent * instance = foundInstance->second;
            _entityIDsToInstancedMeshes[elementType] = instance;
        }
    }
}

void URegionGrowingComponent::Generate()
{
      
	while(!_mainThreadWork.empty())
	{
		std::function<void()> f;

		_mainThreadWork.pop( f );

		f();
	}
        
    if( _generationWorkDone )
    {
        _generationWorkDone();
        _generationWorkDone = nullptr;
    }

	if(!_readyToSubmitNextGenerateCall)
		return;

	std::clock_t start = std::clock();
        
    auto world = GetWorld();
    if(world == nullptr)
        return;
        
        
    FVector cameraLocation;
    FRotator cameraRotation;
    FVector cameraDirection;
        
        
    SceneViewAndFamily::viewLocationRotation(GetWorld(), cameraLocation, cameraRotation);
        
    cameraDirection = cameraRotation.RotateVector(FVector::ForwardVector);
        
        
    UDepthReader * reader = depthReader();
    if( reader )
    {
        _occlusionTester->reader = reader;
        _occlusionTester->cameraLocation = cameraLocation;
		_occlusionTester->_meshInterface = _meshInterface;

        //reader->updateBuffer();
    }
        
    _readyToSubmitNextGenerateCall = false;
        
    float radius = synthesizeDistance;

    Algorithm::AABB limits = _getGenerationLimits();

    // determine our starting points for generation from a grid projected out
    // into the scene
    std::vector<PositionFace> startingPositions;
        
    if( generationMode == EGenerationMode::SpaceFilling )
        startingPositions = _spaceFillingStartingPositions(limits, cameraLocation);
    else if( generationMode == EGenerationMode::SurfaceProjection || generationMode == EGenerationMode::SurfaceWalking)
        startingPositions = _surfaceStartingPositions(limits, cameraLocation, cameraDirection);
    else if( generationMode == EGenerationMode::SurfacePainting )
    {
		if(_accumulatedStartingPositions.size())
			startingPositions.push_back( _accumulatedStartingPositions[0] );
    }

	// copy the toRemove
	auto toRemove = _toRemove;
	_toRemove.clear();
            
    _generationWorker.push( [this,startingPositions,limits,radius,toRemove](int threadID) mutable
    {            
		algorithm->disableOptimization = disableOptimization;

		algorithm->generationMode = generationMode;
		algorithm->useTypeVoting = useTypeVoting;

		// remove points
		AlgorithmResult removalResult;
		{
			std::set<FGraphNodeHandle> removedElements;
			for(auto& p : toRemove)
			{
				auto erasedThisTime = algorithm->eraseAt( p );

				removedElements.insert( erasedThisTime.begin(), erasedThisTime.end() );
			}

			removalResult.removed.insert(removalResult.removed.end(), removedElements.begin(), removedElements.end() );
		}

		// then generate
		AlgorithmResult generationResult;

        generationResult = algorithm->generate(startingPositions, radius, limits);

		// append the generation and removal results
		generationResult.append( removalResult );
            
		RoundRecord round = algorithm->roundSummary().back();

		// do the calculations here, as the roundSummary duration could change in the meantime
		double roundDuration = round.duration();
		double instantaneousGenerationRate = double( round.elementsGenerated ) / round.duration();
		double totalDuration = algorithm->roundSummary().duration();

		double totalGenerationRate = double( round.elementsTotal ) / totalDuration;

        _generationWorkDone = [this, 
			generationResult, 
			round, 
			instantaneousGenerationRate, 
			totalGenerationRate, 
			totalDuration, 
			roundDuration] ()
        {
            _loadResult(generationResult);

			_mainThreadGenerationWorkDone();

			_updateShowMatching();

			if(round.elementsGenerated)
				UE_LOG( LogTemp, Warning, TEXT( "instantaneous generation rate %f" ), instantaneousGenerationRate );

			if(round.elementsTotal)
				UE_LOG( LogTemp, Warning, TEXT( "total generation rate %f total time %f" ), totalGenerationRate, totalDuration );

			_duration += roundDuration;

			_readyToSubmitNextGenerateCall = true;
        };
    });
    
    statElementCount = algorithm->output().size();
}

void URegionGrowingComponent::GlobalOptimization()
{
    AlgorithmResult result = algorithm->globalOptimization();
    
    _loadResult(result);
}

auto URegionGrowingComponent::_spaceFillingStartingPositions(Algorithm::AABB& limits, FVector& cameraLocation) -> std::vector<PositionFace>
{
	std::vector<PositionFace> startingPositions;

	if(useSeeds)
	{
		for(FVector p : seeds)
		{
			startingPositions.emplace_back( eigen( p ), FSurfaceIndex::OffSurface );
		}
	}
	else
	{
		SceneViewAndFamily sceneViewAndFamily = SceneViewAndFamily( GetWorld() );
		FSceneView * sceneView = sceneViewAndFamily.sceneView;

		float radiusSquared = synthesizeDistance * synthesizeDistance;
		float resetSquared = resetDistance * resetDistance;

		if(sceneView != nullptr)
		{
			FIntRect rect = sceneView->ViewRect;
			FIntPoint size = rect.Size();

			float dx = float( size.X ) / 3.0f;
			float dy = float( size.Y ) / 3.0f;

			if(dx < 10.0f)
				dx = 10.0f;

			if(dy < 10.0f)
				dy = 10.0f;

			for(float x = rect.Min.X + dx; x < float( size.X ); x += dx)
			{
				for(float y = rect.Min.Y + dy; y < float( size.Y ); y += dy)
				{
					FVector2D screenPosition = FVector2D( x, y );

					FVector worldOrigin;
					FVector worldDirection;

					sceneView->DeprojectFVector2D( screenPosition, worldOrigin, worldDirection );

					auto hit = limits.intersection( eigen( worldOrigin ), eigen( worldDirection ) );

					Eigen::Vector3f point = eigen( worldOrigin ) + eigen( worldDirection ) * (hit.tMin * 1.0001);
					float distance = (point - eigen( cameraLocation )).squaredNorm();

					if(hit.hit && distance < radiusSquared && distance < resetSquared)
						startingPositions.emplace_back( point, FSurfaceIndex::OffSurface ); // space filling naturally doesn't have faces
				}
			}
		}

		if(limits.contains( eigen( cameraLocation ) ))
			startingPositions.emplace_back( eigen( cameraLocation ), FSurfaceIndex::OffSurface);
	}
    
    return startingPositions;
}

auto URegionGrowingComponent::_surfaceStartingPositions(Algorithm::AABB& limits, FVector& cameraLocation, FVector& cameraDirection) 
-> std::vector<PositionFace>
{
    std::vector<PositionFace> startingPositions;
    
    SceneViewAndFamily sceneViewAndFamily = SceneViewAndFamily(GetWorld());
    FSceneView * sceneView = sceneViewAndFamily.sceneView;
    
    float radiusSquared = synthesizeDistance * synthesizeDistance;
    float resetSquared = resetDistance * resetDistance;
    
    if( sceneView != nullptr && rayCastSeeds )
    {
        FIntRect rect = sceneView->ViewRect;
        FIntPoint size = rect.Size();
        
		if(rayCastDivisions == 1)
		{
			float x = rect.Min.X + size.X / 2;
			float y = rect.Max.Y + size.Y / 2;

			FVector2D screenPosition = FVector2D( x, y );

			FVector worldOrigin;
			FVector worldDirection;

			sceneView->DeprojectFVector2D( screenPosition, worldOrigin, worldDirection );


			auto intersection = _meshInterface->getIntersectionAndFace( worldOrigin, worldDirection );

			tcodsMeshInterfaceBase::SurfacePoint& surfacePoint = intersection.second;

			if( intersection.first )
				startingPositions.emplace_back( eigen(surfacePoint.point), surfacePoint.surfaceIndex );
		}
		else
		{
			float dx = float( size.X ) / rayCastDivisions;
			float dy = float( size.Y ) / rayCastDivisions;

			if(dx < 10.0f)
				dx = 10.0f;

			if(dy < 10.0f)
				dy = 10.0f;

			for(float x = rect.Min.X + dx; x < rect.Max.X; x += dx)
			{
				for(float y = rect.Min.Y + dy; y < rect.Max.Y; y += dy)
				{
					FVector2D screenPosition = FVector2D( x, y );

					FVector worldOrigin;
					FVector worldDirection;

					sceneView->DeprojectFVector2D( screenPosition, worldOrigin, worldDirection );


					auto intersection = _meshInterface->getIntersectionAndFace( worldOrigin, worldDirection );

					if(!intersection.first)
						continue;

					tcodsMeshInterfaceBase::SurfacePoint& surfacePoint = intersection.second;

					startingPositions.emplace_back(eigen(surfacePoint.point), surfacePoint.surfaceIndex);
				}
			}
		}
    }
    
    
    


	auto nearest = _meshInterface->nearestPointOnMesh( cameraLocation );
	if(limits.contains( eigen( nearest.point ) ))
	{
		Eigen::Vector3f nearestPoint = eigen( nearest.point );
		Eigen::Vector3f camera = eigen( cameraLocation );

		float distance = (camera - nearestPoint).squaredNorm();

		if(distance < synthesizeDistance * synthesizeDistance && distance < resetDistance * resetDistance)
			startingPositions.emplace_back( nearestPoint, nearest.surfaceIndex );
	}
   
    return startingPositions;
}



UStaticMeshComponent* URegionGrowingComponent::staticMeshComponent()
{
    return GetOwner()->FindComponentByClass<UStaticMeshComponent>();
}

URuntimeMeshComponent * URegionGrowingComponent::debugMeshComponent() // for debug drawing
{
    if(_debugMesh == nullptr )
    {
        _debugMesh = NewObject<URuntimeMeshComponent>(GetOwner());
        _debugMesh->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
        _debugMesh->SetWorldScale3D(this->GetOwner()->GetActorScale3D());
        _debugMesh->SetMaterial(0, directionFieldMaterial);
        _debugMesh->RegisterComponent();
    }
    
    UStaticMeshComponent * mesh = staticMeshComponent();
    if( mesh )
        _debugMesh->SetWorldTransform(mesh->GetComponentToWorld());
    
    return _debugMesh; 
}

void URegionGrowingComponent::_initRuntimeMesh() // for debug drawing
{
    meshInterfaceRuntimeMesh = NewObject<URuntimeMeshComponent>( GetOwner() );
    meshInterfaceRuntimeMesh->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
    meshInterfaceRuntimeMesh->SetWorldScale3D( this->GetOwner()->GetActorScale3D() );
    // DISABLED this, we don't use that depth-actor visibility stuff anymore
	//runtimeMesh->bRenderCustomDepth = true;
	meshInterfaceRuntimeMesh->SetShouldSerializeMeshData(false);
    meshInterfaceRuntimeMesh->RegisterComponent();

    UStaticMeshComponent * mesh = staticMeshComponent();
    if(mesh)
        meshInterfaceRuntimeMesh->SetWorldTransform( mesh->GetComponentToWorld() );

	meshInterfaceRuntimeMesh->SetMaterial( 0, debugSurfaceMaterial );
}

UStaticMeshComponent* URegionGrowingComponent::nearestPointDebugMesh()
{
    if( _nearestPointDebugMesh == nullptr )
    {
        
        UStaticMeshComponent * meshComponent = staticMeshComponent();
        
        if( meshComponent == nullptr )
            return nullptr;
        
        const auto transform = meshComponent->GetComponentTransform();
        const auto bounds = meshComponent->Bounds;
        
        auto spherePath = FName(TEXT("/Engine/BasicShapes/Sphere"));
        UStaticMesh * sphere = LoadObjectFromPath<UStaticMesh>(spherePath);
        
        auto materialPath = FName(TEXT("/Engine/BasicShapes/BasicShapeMaterial"));
        UMaterial * material = LoadObjectFromPath<UMaterial>(materialPath);
        UMaterialInstanceDynamic * materialInstance = UMaterialInstanceDynamic::Create(material, this);
        
        materialInstance->SetVectorParameterValue(FName(TEXT("Color")), FLinearColor(0.0f, 1.0f, 0.0f, 1.0f));
        
        const float pointSize = (bounds.SphereRadius / 50.0f) / 1000.0f; // sphere radius is 100, and we'd like about 100 spheres to fit across the mesh
        
        _nearestPointDebugMesh = NewObject<UStaticMeshComponent>(GetOwner());
        
        _nearestPointDebugMesh->SetStaticMesh(sphere);
        _nearestPointDebugMesh->SetMobility(EComponentMobility::Movable);
        _nearestPointDebugMesh->RegisterComponent();
		_nearestPointDebugMesh->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
        
        _nearestPointDebugMesh->SetWorldScale3D(FVector(pointSize));
        _nearestPointDebugMesh->SetMaterial(0, materialInstance);
    }
    
    return _nearestPointDebugMesh;
}


void URegionGrowingComponent::loadParameters()
{
    if( !_didInit && !_initMeshInterface() )
        return;

	algorithm->perElementParameters = perElementParameters;
    
	algorithm->enableRoundSummaries = enableRoundSummaries;
	algorithm->calculate_kCoherenceEnergy = calculate_kCoherenceEnergy;
	algorithm->calculate_bruteForceEnergy = calculate_bruteForceEnergy;
	algorithm->cellExtents = cinpactCellSize;
	algorithm->useCinpactEnergy = useCinpactEnergy;

    algorithm->generationMode = generationMode;
	algorithm->sketchBasedForceTangentPlane = sketchBasedForceTangentPlane;
    algorithm->generationParameters = generationParameters;
    algorithm->optimizationParameters = optimizationParameters;
	algorithm->useTypeVoting = useTypeVoting;
	algorithm->removeOverlaps = removeOverlaps;
	algorithm->enableReassignment = enableReassignment;
    algorithm->generationInnerRadius = generationRadius;
	algorithm->frozenElementRadius = frozenElementRadius;
    algorithm->kCoherence = kCoherence;
    algorithm->minAssignmentDistance = minAssignmentDistance;
    algorithm->freespaceRadius = freespaceRadius;
    algorithm->typeCost = typeCost;
	algorithm->gradientTerm = gradientTerm;
    algorithm->relaxation = relaxation;
    algorithm->flipSurfaceNormals = flipSurfaceNormals;
    algorithm->ignoreBadSuggestions = ignoreBadSuggestions;
	algorithm->ignoreBadSuggestionsDistanceFactor = ignoreBadSuggestionsDistanceFactor;

	algorithm->kCoherenceClusteringPenalty = kCoherenceClusteringPenalty;
	algorithm->kCoherenceClusteringPenaltyRadius = kCoherenceClusteringPenaltyRadius;

	algorithm->nThreads = nThreads;

	algorithm->enableVolumeSurfaceInteraction = enableVolumeSurfaceInteraction;

	algorithm->disableOptimization = disableOptimization;
	algorithm->optimizationRounds = optimizationRounds;
	algorithm->useGlobalOptimization = useGlobalOptimization;

	algorithm->forceRotation = forceRotation;
	algorithm->forcedRotation = forcedRotation;

	algorithm->seedsIgnoreMeshInBoundary = seedsIgnoreMeshInBoundary;
    
    algorithm->sourceHistogramRadius = sourceHistogramRadius;
    algorithm->sourceHistogramWeight = sourceHistogramWeight;
    algorithm->samplingDistanceWeight = samplingDistanceWeight;

	algorithm->patchInitializationExemplarDivisions = patchInitializationExemplarDivisions;
	algorithm->patchMarginPercentage = patchMarginPercentage;
   
    // We need to convert our seeds into Algorithm space. Right now
    // algorithm space is just Unreal's world space.
    UStaticMeshComponent * staticMesh = staticMeshComponent();
    FTransform transform = staticMesh->GetComponentToWorld();
    
    algorithm->setMeshInterface(_meshInterface);
    
    if( depthReader() )
        _occlusionTester = std::shared_ptr<RGOcclusionTester>(new RGOcclusionTester());
    else
        _occlusionTester = nullptr;
    
    algorithm->setOcclusionTester(_occlusionTester);
}

auto URegionGrowingComponent::_computeInstanceTransform( FGraphNode& element ) -> FTransform
{
	FVector newLocation = element.position;

	FQuat rotation = element.orientation;
	FVector scale = FVector(element.scale);

	return FTransform( rotation, newLocation, scale );
}

void URegionGrowingComponent::_updateInstance(FGraphNodeHandle handle)
{
	if(_cache.find(handle) == _cache.end())
		return;

	ElementTuple element(handle, _synthesisContext->graph());

	FTransform transform = _computeInstanceTransform(element.node);

	auto& cached = _cache[handle];
     
    cached.transform = transform;

	if (cached.instance < 0)
		return;
    
	if(cached.type == element.element.type)
	{
		_instanceManager.updateInstance( cached.instance, cached.type, transform );
	}
	else
	{
		_instanceManager.removeInstance( cached.instance, cached.type );

		auto newInstance = _instanceManager.addInstance(element.element.type, transform );

		cached.instance = newInstance;
		cached.type = element.element.type;
	}
}

void URegionGrowingComponent::_loadResult(const AlgorithmResult& result)
{
    auto sceneComponent = output()->GetRootComponent();
    
    FActorSpawnParameters spawningParameters;
    spawningParameters.Owner = GetOwner();
    
    auto& outputDomain = algorithm->output();
    std::vector<FGraphNodeHandle> elements;
    elements.reserve(outputDomain.size());
    for( auto& e : outputDomain )
        elements.push_back(e.nodeHandle());
    
    _elementKDTree.setElements(elements);
    
    bool dirty = false;
        
    std::vector<FGraphNodeHandle> generatedElements;
    std::vector<FGraphNodeHandle> modifiedElements;
	std::vector<FGraphNodeHandle> removedElements;
        

    generatedElements = std::move(result.generated);
    modifiedElements = std::move(result.modified);
	removedElements = std::move( result.removed );
    
	FGraph& graph = _synthesisContext->graph();
        
    for( FGraphNodeHandle handle : generatedElements )
    {
		ElementTuple generated(handle, graph);

		FTransform transform = _computeInstanceTransform(generated.node);
            		
		int32 instance = _instanceManager.addInstance(generated.element.type, transform );
            
        auto& cached = _cache[handle];
        cached.type = generated.element.type;
        cached.instance = instance;
        cached.transform = transform;
            
        dirty = true;
    }
        
    for( FGraphNodeHandle modified : modifiedElements )
    {
        _updateInstance(modified);
            
        dirty = true;
    }

	for(FGraphNodeHandle removed : removedElements)
	{
		auto found = _cache.find( removed );

		if(found == _cache.end())
			continue;

		auto& cached = found->second;

		_instanceManager.removeInstance( cached.instance, cached.type );

		_cache.erase( found );

		dirty = true;
	}
             
    // we need to let unreal know that we've made changes to each of the
    // ``UInstanceStaticMeshComponent''s.
        
    if( dirty )
    {
		_instanceManager.markDirty();
            
        for( auto & pair : _entityIDsToInstancedMeshes )
            pair.second->MarkRenderStateDirty();
    }
}

void URegionGrowingComponent::_destroyChildren(AActor * actor)
{
    if( !actor->GetRootComponent() )
        return;
    
    TArray<USceneComponent*> children;
    actor->GetRootComponent()->GetChildrenComponents(false, children);
    
    auto world = GetWorld();
    bool editor = GetWorld()->WorldType == EWorldType::Editor;
    
    std::vector<ULevel*> levels;
    
    for( auto child : children )
    {
        AActor * childActor = child->GetOwner();
        
        levels.push_back(actor->GetLevel());
        
        if( editor )
            world->EditorDestroyActor(childActor, false);
        else
            world->DestroyActor(childActor, false, false /* we'll call modify level */ );
    }
    
    for( auto level : levels )
        world->ModifyLevel(level);
}

#if WITH_EDITOR
void URegionGrowingComponent::PostEditChangeProperty(struct FPropertyChangedEvent& e)
{
    FName propertyName = (e.Property != NULL) ? e.Property->GetFName() : NAME_None;
    
    if( propertyName == GET_MEMBER_NAME_CHECKED(URegionGrowingComponent, drawDistanceField) )
    {
        _updateDrawDirectionField();
        _updateSeeds();
    }

    
	if(propertyName == GET_MEMBER_NAME_CHECKED( URegionGrowingComponent, drawHalfEdge ))
	{
		_updateDrawHalfEdge();
	}

	if( propertyName == GET_MEMBER_NAME_CHECKED( URegionGrowingComponent, drawBounds ) ||
		propertyName == GET_MEMBER_NAME_CHECKED( URegionGrowingComponent, generationLimitsMin ) ||
		propertyName == GET_MEMBER_NAME_CHECKED( URegionGrowingComponent, generationLimitsMax )		)
	{
		_updateDrawBounds();
	}
    
    Super::PostEditChangeProperty(e);
}
#endif

// -----------------------------------------------------------------------------
/// Draw Distance Field
// -----------------------------------------------------------------------------

void URegionGrowingComponent::_updateDirectionField()
{
	using namespace tcods;
	using namespace std;

	if (!_meshInterface)
		return;

	for (auto& pair : _meshInterface->_sections)
	{
		Mesh& mesh = *(pair.second.get());

		mesh.fieldAngle = 0.0;

		// we need to reset all of our k's
		for (auto& v : mesh.vertices)
			v.k = 0.0;
	}



	for (const auto singularity : singularities)
	{
		if ( !_meshInterface->hasSection(singularity.sectionIndex) )
		{
			UE_LOG(LogTemp, Warning, TEXT("Warning: singularity references a mesh section that does not exist, %d"), singularity.sectionIndex );
			continue;
		}
		
		Mesh& mesh = _meshInterface->mesh(singularity.sectionIndex);



		int V = mesh.vertices.size();

		int i = singularity.vertexIndex;
		double k = singularity.k;

		if (i < 0 || i >= V)
		{
			UE_LOG(LogTemp, Warning, TEXT("Warning: singularity requested at vertex %d (mesh has only %d vertices!)"), i + 1, V);
		}
		else
		{
			mesh.vertex(i)->k = k;
		}
	}

	//vector<Generator> generators;

	//int G = mesh.nGenerators();
	//for (const auto g : generators)
	//{
	//	int i = g.first;
	//	double k = g.second;

	//	if (i < 0 || i >= G)
	//	{
	//		UE_LOG(LogTemp, Warning, TEXT("Warning: additional holonomy requested around generator %d (mesh has only %d generators!)"), i + 1, G);
	//	}
	//	else
	//	{
	//		mesh.generatorIndices[i] = k;
	//	}
	//}



	if (buildTrivialConnections)
	{
		for (auto& pair : _meshInterface->_sections)
		{
			pair.second->computeTrivialConnection();
		}
	}
}



bool URegionGrowingComponent::_initMeshInterface()
{
	if (_meshInterface)
		return true;
    
	if (meshInterfaceMode == EMeshInterfaceMode::StaticMesh)
	{
		auto staticMeshInterface = std::make_shared<tcodsUStaticMeshInterface>();

		_meshInterface = staticMeshInterface;
		_meshInterface->sampleSubdivisions = sampleSubdivisions;

		UStaticMeshComponent * meshComponent = staticMeshComponent();

		if (meshComponent == nullptr)
			return false;

		UStaticMesh * uMesh = meshComponent->GetStaticMesh();

		if (uMesh == nullptr)
			return false;

		FTransform toWorld = meshComponent->GetComponentTransform();

		staticMeshInterface->buildMesh(*uMesh, toWorld, generationBoundsWorld(), meshInterfaceRuntimeMesh);

		meshInterfaceRuntimeMesh->SetMaterial(0, meshComponent->GetMaterial(0));

	}
	else if (meshInterfaceMode == EMeshInterfaceMode::ChunkedMesh)
	{
		auto chunkMeshInterface = std::make_shared<tcodsChunkedGridMeshInterface>();

		_meshInterface = chunkMeshInterface;

		if (!chunkedStaticMeshActor)
			return false;

		//updateMeshInterface();
	}

	return true;
}

void URegionGrowingComponent::_updateDrawDirectionField()
{
    if( drawDistanceField )
        _drawDirectionField();
    else
        _undrawDirectionField();
}

void URegionGrowingComponent::_drawDirectionField()
{
    using namespace tcods;
    using namespace std;
	using namespace DDG;
    
    if( !_didInit )
        return;
    
    UStaticMeshComponent * meshComponent = staticMeshComponent();
    
    if( meshComponent == nullptr )
        return;

	const auto transform = meshComponent->GetComponentTransform();
	const auto invTransform = transform.Inverse();
	const auto bounds = meshComponent->Bounds;
	float scale = (bounds.SphereRadius / 100.0f) * 1.2f;

	TArray<FVector> outVertices;
	TArray<int32> outIndices;
	int32 currentIndex = 0;

	// draw arrows for the direction fields
	for (auto& pair : _meshInterface->_sections)
	{
		uint32_t section = pair.first;
		Mesh& mesh = *pair.second.get();

		const vector<Face>& faces(mesh.faces);

		for (FaceCIter i = faces.begin(); i != faces.end(); i++)
		{
			HalfEdgeIter he = i->he;

			// don't write vectors on boundary faces
			if (he->onBoundary)
			{
				continue;
			}

			Vector p = i->barycenter();

			double alpha = i->alpha;

			Vector e1, e2; i->frame(e1, e2);
			Vector n = cross(e2, e1);
			Vector u = e1 * cos(alpha) + e2 * sin(alpha);

			// the vector orthogonal to u in the plane of e1, e2
			const double pi_2 = _pi() * .5;
			Vector uo = e1 * cos(alpha + pi_2) + e2 * sin(alpha + pi_2);

			Vector o = p + n * scale * .1; // lift the arrow off of the surface a bit
			Vector p0 = o + u * scale;
			Vector p1 = o + uo * scale * .1;
			Vector p2 = o - uo * scale * .1;

			// The verices that we get from the above operations are already in world
			// space, but because our procedural mesh component is a child of our rootComponent
			// we need to apply our inverse transform to them.
			const FVector p0_ = invTransform.TransformPosition(unreal(p0));
			const FVector p1_ = invTransform.TransformPosition(unreal(p1));
			const FVector p2_ = invTransform.TransformPosition(unreal(p2));

			outVertices.Add(p0_);
			outVertices.Add(p1_);
			outVertices.Add(p2_);

			outIndices.Add(currentIndex++);
			outIndices.Add(currentIndex++);
			outIndices.Add(currentIndex++);
		}
	}

	URuntimeMeshComponent * procedural = debugMeshComponent();
	procedural->CreateMeshSection(MeshSections::DirectionField, outVertices, outIndices, TArray<FVector>() /* normals */, TArray<FVector2D>() /* uvs*/, TArray<FColor>() /* colors */, TArray<FRuntimeMeshTangent>(), false);
	procedural->UpdateBounds();
    
   
    
    // remove the unused singularity meshes
    for( auto component : singularityDebugSpheres )
    {
        component->RemoveFromRoot();
        component->UnregisterComponent();
        component->DestroyComponent();
    }
    
    singularityDebugSpheres.Empty();
    
    // add new singularity meshes
    auto spherePath = FName(TEXT("/Engine/BasicShapes/Sphere"));
    UStaticMesh * sphere = LoadObjectFromPath<UStaticMesh>(spherePath);
    
    auto materialPath = FName(TEXT("/Engine/BasicShapes/BasicShapeMaterial"));
    UMaterial * material = LoadObjectFromPath<UMaterial>(materialPath);
    UMaterialInstanceDynamic * materialInstance = UMaterialInstanceDynamic::Create(material, this);
    
    materialInstance->SetVectorParameterValue(FName(TEXT("Color")), FLinearColor(1.0f, 0.0f, 0.0f, 1.0f));
    
    const float pointSize = (bounds.SphereRadius / 50.0f) / 100.0f; // sphere radius is 100, and we'd like about 100 spheres to fit across the mesh
    
    for( auto& singularity : singularities )
    {
		Vertex& v = *_meshInterface->mesh(singularity.sectionIndex).vertex(singularity.vertexIndex);
        
        UStaticMeshComponent * sphereComponent = NewObject<UStaticMeshComponent>(GetOwner());
        
        sphereComponent->SetStaticMesh(sphere);
        sphereComponent->SetMobility(EComponentMobility::Movable);
        sphereComponent->RegisterComponent();
		sphereComponent->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
        
        sphereComponent->SetWorldLocation(unreal(v.position));
        sphereComponent->SetWorldScale3D(FVector(pointSize));
        sphereComponent->SetMaterial(0, materialInstance);
        
        singularityDebugSpheres.Add(sphereComponent);
    }
    
}

void URegionGrowingComponent::_undrawDirectionField()
{
	URuntimeMeshComponent * procedural = debugMeshComponent();
    procedural->ClearMeshSection( MeshSections::DirectionField );
    
    for( auto component : singularityDebugSpheres )
    {
        component->RemoveFromRoot();
        component->UnregisterComponent();
        component->DestroyComponent();
    }
    
    singularityDebugSpheres.Empty();
}

void URegionGrowingComponent::addSingularity(const FVector& vertex, float k)
{
    if( !_didInit && !_initMeshInterface() )
        return;
    
    int32 index;
    if( !_meshInterface->indexOfVertex(vertex, index) )
        return;
    
    singularities.Emplace(index, k);
    
    _updateDirectionField();
    _updateDrawDirectionField();
}

void URegionGrowingComponent::removeSingularityAt(unsigned int index)
{
    singularities.RemoveAt(index);
    
    if( !_didInit && !_initMeshInterface() )
        return;
    
    _updateDirectionField();
    _updateDrawDirectionField();
}

// -----------------------------------------------------------------------------
/// - Seeds
// -----------------------------------------------------------------------------

void URegionGrowingComponent::addSeed(const FVector& localPoint)
{
    seeds.Emplace(localPoint);
    
    _updateSeeds();
}

void URegionGrowingComponent::removeSeedAt(unsigned int index)
{
    seeds.RemoveAt(index);
    
    _updateSeeds();
}

UInstancedStaticMeshComponent* URegionGrowingComponent::seedPointsInstancedMeshComponent()
{
    if( _seedPointsInstancedMeshComponent == nullptr )
    {
        UStaticMeshComponent * meshComponent = staticMeshComponent();
        
        if( meshComponent == nullptr )
            return nullptr;
        
        const auto localTransform = meshComponent->GetRelativeTransform();
        const auto bounds = meshComponent->Bounds;
        
        auto spherePath = FName(TEXT("/Engine/BasicShapes/Sphere"));
        UStaticMesh * sphere = LoadObjectFromPath<UStaticMesh>(spherePath);
        
        auto materialPath = FName(TEXT("/Engine/BasicShapes/BasicShapeMaterial"));
        UMaterial * material = LoadObjectFromPath<UMaterial>(materialPath);
        UMaterialInstanceDynamic * materialInstance = UMaterialInstanceDynamic::Create(material, this);
        
        materialInstance->SetVectorParameterValue(FName(TEXT("Color")), FLinearColor(0.0f, 1.0f, 0.0f, 1.0f));
        
        const float pointSize = (bounds.SphereRadius / 50.0f) / 100.0f; // sphere radius is 100, and we'd like about 100 spheres to fit across the mesh
        
        _seedPointsInstancedMeshComponent = NewObject<UInstancedStaticMeshComponent>(GetOwner());
        
        _seedPointsInstancedMeshComponent->SetStaticMesh(sphere);
        _seedPointsInstancedMeshComponent->SetMobility(EComponentMobility::Static);
        _seedPointsInstancedMeshComponent->SetCollisionProfileName(TEXT("NoCollision"));
		_seedPointsInstancedMeshComponent->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );        
        _seedPointsInstancedMeshComponent->SetMaterial(0, materialInstance);

		_seedPointsInstancedMeshComponent->RegisterComponent();
    }
    
    return _seedPointsInstancedMeshComponent;
}

void URegionGrowingComponent::_updateSeeds()
{
    if( drawDistanceField )
        _drawSeeds();
    else
        _destroySeeds();
}

void URegionGrowingComponent::_drawSeeds()
{
    _destroySeeds();
    
    UStaticMeshComponent * meshComponent = staticMeshComponent();
    if( !meshComponent )
        return;
    
    auto instancedMesh = seedPointsInstancedMeshComponent();
    if( !instancedMesh )
        return;
    
    instancedMesh->ClearInstances();
    
    // sphere radius is 100, and we'd like about 100 spheres to fit across the mesh
    const auto bounds = meshComponent->Bounds;
    const float pointSize = (bounds.SphereRadius / 50.0f) / 100.0f; // sphere radius is 100, and we'd like about 100 spheres to fit across the mesh
    
    for( auto& seed : seeds )
    {
        FTransform instanceTransform;
        instanceTransform.SetLocation(seed);
        instanceTransform.SetScale3D(FVector(pointSize));
        
        instancedMesh->AddInstance(instanceTransform);
    }
}

void URegionGrowingComponent::_destroySeeds()
{
    auto instance = seedPointsInstancedMeshComponent();
    if( !instance )
        return;
    
    instance->ClearInstances();
}


// -----------------------------------------------------------------------------
// - Half Edge Drawing
// -----------------------------------------------------------------------------

void URegionGrowingComponent::_updateGenerationLimits()
{
	if(!_boundsBox || drawBounds == false)
		return;

	FVector centre = _boundsBox->GetComponentLocation();
	FVector extents = _boundsBox->GetScaledBoxExtent();

	FVector min = centre - extents;
	FVector max = centre + extents;

	generationLimitsMin = min;
	generationLimitsMax = max;
}

void URegionGrowingComponent::_updateDrawBounds()
{
	if(drawBounds)
	{
		_drawBounds();
		_updateGenerationLimits();
	}
	else
		_undrawBounds();
}

FBox URegionGrowingComponent::generationBoundsWorld()
{
	return FBox( generationLimitsMin, generationLimitsMax );
}



void URegionGrowingComponent::_drawBounds()
{
	if(_boundsBox == nullptr)
	{
		_boundsBox = NewObject<UBoxComponent>( GetOwner() );
		
		AActor * actor = GetOwner();

		USceneComponent * root = actor->GetRootComponent();

		_boundsBox->AttachToComponent( root, FAttachmentTransformRules::KeepRelativeTransform );
		actor->AddInstanceComponent( _boundsBox );
		_boundsBox->RegisterComponent();
	}

	if(!_boundsBox->IsVisible())
		_boundsBox->SetVisibility( true );

	FBox box = FBox( generationLimitsMin, generationLimitsMax );

	_boundsBox->SetWorldScale3D( FVector( 1.0f ) );
	_boundsBox->SetWorldRotation( FQuat::Identity, false, nullptr, ETeleportType::TeleportPhysics );
	_boundsBox->SetWorldLocation( box.GetCenter() );
	_boundsBox->SetBoxExtent( box.GetExtent() );
}

void URegionGrowingComponent::_undrawBounds()
{
	if(_boundsBox)
		_boundsBox->SetVisibility( false );
}

// -----------------------------------------------------------------------------
// - Half Edge Drawing
// -----------------------------------------------------------------------------

void URegionGrowingComponent::_updateDrawHalfEdge()
{
	if(drawHalfEdge)
		_drawHalfEdge();
	else
		_undrawHalfEdge();
}

void URegionGrowingComponent::_drawHalfEdge()
{
	using namespace tcods;
	using namespace std;
	using namespace DDG;

	UStaticMeshComponent * meshComponent = staticMeshComponent();

	if(meshComponent == nullptr)
		return;


	const auto transform = meshComponent->GetComponentTransform();
	const auto invTransform = transform.Inverse();
	const auto bounds = meshComponent->Bounds;
	float scale = (bounds.SphereRadius / 100.0f) * 1.2f;

	TArray<FVector> outVertices;
	TArray<int32> outIndices;
	int32 currentIndex = 0;

	for (auto& pair : _meshInterface->_sections)
	{
		Mesh& mesh = *pair.second.get();

		for (HalfEdge& edge : mesh.halfedges)
		{
			Vector v0 = edge.from->position;
			Vector v1 = edge.next->from->position;

			Face& face = *edge.face;

			Vector e1, e2; face.frame(e1, e2);
			Vector n = cross(e2, e1);

			// draw an arrow
			{
				Vector mid = (v0 + v1) / 2.0f;



				Vector b = face.barycenter();

				Vector lift = n * scale * .1;
				Vector d = (v1 - v0);

				Vector p0 = mid + (b - mid) * .1 + lift;
				Vector p1 = mid;
				Vector p2 = mid + d * 0.1;


				// The vertices that we get from the above operations are already in world
				// space, but because our procedural mesh component is a child of our rootComponent
				// we need to apply our inverse transform to them.
				const FVector p0_ = invTransform.TransformPosition(unreal(p0));
				const FVector p1_ = invTransform.TransformPosition(unreal(p1));
				const FVector p2_ = invTransform.TransformPosition(unreal(p2));

				// two sided
				outVertices.Add(p0_);
				outVertices.Add(p1_);
				outVertices.Add(p2_);

				outVertices.Add(p0_);
				outVertices.Add(p2_);
				outVertices.Add(p1_);

				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);

				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
			}

			// draw a wire-frame
			{
				Vector v2 = edge.next->next->from->position;

				Vector d = (v2 - v0);
				Vector t = cross((v1 - v0), n);

				t = t.unit() * 5.0;

				const FVector p0 = unreal(v0);
				const FVector p1 = unreal(v0 + t);
				const FVector p2 = unreal(v1 + t);
				const FVector p3 = unreal(v0 + t);

				const FVector p0_ = invTransform.TransformPosition(p0);
				const FVector p1_ = invTransform.TransformPosition(p1);
				const FVector p2_ = invTransform.TransformPosition(p2);
				const FVector p3_ = invTransform.TransformPosition(p3);


				outVertices.Add(p0_);
				outVertices.Add(p1_);
				outVertices.Add(p2_);
				outVertices.Add(p0_);
				outVertices.Add(p2_);
				outVertices.Add(p1_);


				outVertices.Add(p0_);
				outVertices.Add(p2_);
				outVertices.Add(p3_);
				outVertices.Add(p0_);
				outVertices.Add(p3_);
				outVertices.Add(p2_);

				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
				outIndices.Add(currentIndex++);
			}
		}
	}


	URuntimeMeshComponent * procedural = debugMeshComponent();
	procedural->CreateMeshSection( MeshSections::HalfEdge, outVertices, outIndices, TArray<FVector>() /* normals */, TArray<FVector2D>() /* uvs*/, TArray<FColor>() /* colors */, TArray<FRuntimeMeshTangent>(), false );
	procedural->UpdateBounds();
}

void URegionGrowingComponent::_undrawHalfEdge()
{
	URuntimeMeshComponent * procedural = debugMeshComponent();
	procedural->ClearMeshSection( MeshSections::HalfEdge );
}

void URegionGrowingComponent::updateMeshInterface()
{
	if (meshInterfaceMode == EMeshInterfaceMode::ChunkedMesh)
	{
		// We don't do a dynamic_cast because it is disabled by default in Unreal.
		// It can be turned on with bUseRTTI = true in the build.cs file.
		tcodsChunkedGridMeshInterface * chunkInterface = StaticCast<tcodsChunkedGridMeshInterface*>(_meshInterface.get());

		if ( !chunkedStaticMeshActor)
			return;

		UChunkedVolumeComponent * volumeComponent = chunkedStaticMeshActor->FindComponentByClass<UChunkedVolumeComponent>();

		ChunkGrid<float>& grid = volumeComponent->grid();

		chunkInterface->buildMesh(grid, chunkedStaticMeshActor->ActorToWorld(), volumeComponent->isoLevel, generationBoundsWorld(), meshInterfaceRuntimeMesh);

		// set materials
		for (auto& section : chunkInterface->_sections)
		{
			meshInterfaceRuntimeMesh->SetMaterial(section.first, volumeComponent->material);
			meshInterfaceRuntimeMesh->SetMeshSectionVisible(section.first, drawSurface);
		}

		if (drawSurface)
		{
			URuntimeMeshComponent * volumeMesh = chunkedStaticMeshActor->FindComponentByClass<URuntimeMeshComponent>();

			if (volumeMesh)
				volumeMesh->SetVisibility(false);
		}

		_updateDirectionField();
		_updateDrawDirectionField();
	}
}

ASimulationSnapshotActor* URegionGrowingComponent::createGraphicalSnapshotActor(UWorld * targetWorld)
{
	ASimulationSnapshotActor * newActor = targetWorld->SpawnActor<ASimulationSnapshotActor>();

	USceneComponent * newScene = NewObject<USceneComponent>( newActor );
	{
		AActor * originalActor = GetOwner();

		USceneComponent * originalScene = originalActor->GetRootComponent();

		FTransform newTransform = originalScene->GetComponentTransform();

		newActor->SetRootComponent( newScene );
		newActor->AddInstanceComponent( newScene );

		newScene->SetWorldTransform( newTransform );
	}
    
    
    for( auto& typeInstancePair : _instanceManager.elementTypesToInstancedMeshes )
    {
        UInstancedStaticMeshComponent * instance = typeInstancePair.second;

		// I tried DuplicateObject but the ISMC just wasn't showing up after the duplication, so we'll do it all manually
        UInstancedStaticMeshComponent * newInstance = NewObject<UInstancedStaticMeshComponent>( newActor );
		
		newInstance->SetStaticMesh( instance->GetStaticMesh() );

		for(int32 mi = 0; mi < instance->GetNumMaterials(); mi++)
			newInstance->SetMaterial( mi, instance->GetMaterial( mi ) );
        
		newInstance->SetCollisionEnabled( ECollisionEnabled::NoCollision );
		newInstance->SetCollisionProfileName( TEXT( "NoCollision" ) );
		newInstance->AttachToComponent( newScene, FAttachmentTransformRules::KeepRelativeTransform );
        
		newActor->AddInstanceComponent( newInstance );

		newInstance->RegisterComponent();

		// copy the instances
		int32 n = instance->GetInstanceCount();
		for(int32 i = 0; i < n; i++)
		{
			FTransform transform;

			instance->GetInstanceTransform( i, transform, true );

			newInstance->AddInstanceWorldSpace( transform );
		}
    }

    return newActor;
}

OutputDomainExport URegionGrowingComponent::exportOutputDomain()
{
	// wait for work to be finished
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	OutputDomainExport domainExportInfo;

	Domain& outputDomain = algorithm->output();

	FTransform toWorld = FTransform::Identity; // we are already in world space
	domainExportInfo.read(outputDomain.graph, toWorld);

	domainExportInfo.meshInterface = _meshInterface;

	return domainExportInfo;
}

void URegionGrowingComponent::_loadTypeDescriptionsFromNewElementsAndSetTypes( const std::vector<FGraphNodeHandle>& elements, FGraph& graph )
{
	// we need to map to the exemplar
	std::map<ElementTypeDescription, int, InferredElementTypeDescriptionComparator> types;

	int nextType = 0;
	for(auto& pair : _elementTypesToDescriptions)
	{
		types.emplace( pair.second, pair.first );

		if(pair.first > nextType)
			nextType = pair.first;
	}

	nextType++;

	bool newTypesFound = false;

	// look for new types
	const UScriptStruct * meshGraphObjectScriptStruct = FGraphMesh::StaticStruct();
	for(FGraphNodeHandle  e : elements)
	{
		ElementTypeDescription typeDescription;

		ElementTuple tuple(e, graph);


		if (tuple.node.hasComponent<FGraphMesh>())
		{
			FGraphMesh& mesh = graph.component<FGraphMesh>(e);

			typeDescription.mesh = mesh.staticMesh;
			typeDescription.material = mesh.material;

			if (mesh.material)
				typeDescription.materials.Add(mesh.material);
		}

		auto foundType = types.find( typeDescription );

		// new type time!
		if(foundType == types.end())
		{
			newTypesFound = true;

			types[typeDescription] = nextType;

			// HACK
			// set the element type
			tuple.element.type = nextType;

			nextType++;
		}
	}

	// update the _elementTypesToDesciptions tables
	_elementTypesToDescriptions.clear();

	for(auto& pair : types)
	{
		_elementTypesToDescriptions[pair.second] = pair.first;
	}

	// update the ISMCs
	_loadInstancedStaticMeshes();
}

void URegionGrowingComponent::_reassignFaceIndices(std::vector<FGraphNodeHandle>& elements, FGraph& graph)
{
	Element::eachTuple(elements, graph, [&](ElementTuple tuple)
	{
		if (!tuple.element.surfaceIndex.isOnSurface())
			return;

		FVector p = tuple.node.position;

		auto nearestPoint = _meshInterface->nearestPointOnMesh(p);

		tuple.element.surfaceIndex = nearestPoint.surfaceIndex;
		tuple.node.position = nearestPoint.point;
	});
}

void URegionGrowingComponent::loadElementDomain(FGraphSnapshot& toLoad )
{
	// wait for work to be finished
	_generationWorker.push( []( int threadID ) mutable {} ).get();

	// clear everything
	_privateClear();

	// update our ElementTypeDescriptions and faces
	{
		auto toLoadHandles = toLoad.graph.nodeHandles();

		_loadTypeDescriptionsFromNewElementsAndSetTypes( toLoadHandles, toLoad.graph );

		_reassignFaceIndices( toLoadHandles, toLoad.graph );
	}

	_readyToSubmitNextGenerateCall = false;

	_generationWorker.push( [this, toLoad]( int threadID ) mutable
	{
		// load the result into the algorithm
		algorithm->loadOutputElements( toLoad );

		auto loadedElements = algorithm->output().elements();

		AlgorithmResult result;

		result.generated = loadedElements;

		_generationWorkDone = [this,result]()
		{
			_loadResult( result );
		};

		_readyToSubmitNextGenerateCall = true;
	} );

}

// -----------------------------------------------------------------------------
// - Show matching visualization
// -----------------------------------------------------------------------------

UInstancedStaticMeshComponent& URegionGrowingComponent::debug_matchingISMC()
{
	if (_debug_matchingISMC == nullptr)
	{
		_debug_matchingISMC = NewObject<UInstancedStaticMeshComponent>(GetOwner());

		_debug_matchingISMC->SetStaticMesh(debug_matchingPointMesh);
		_debug_matchingISMC->SetMobility(EComponentMobility::Static);
		_debug_matchingISMC->SetCollisionProfileName(TEXT("NoCollision"));
		_debug_matchingISMC->AttachToComponent(this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
		_debug_matchingISMC->SetMaterial(0, debug_matchingPointMaterial);

		_debug_matchingISMC->RegisterComponent();
	}

	return *_debug_matchingISMC;
}

void URegionGrowingComponent::_updateShowMatching()
{
	UInstancedStaticMeshComponent& ismc = debug_matchingISMC();

	ismc.ClearInstances();

	if (!showMatchingExamples)
		return;

	for (auto& elementObject : algorithm->output())
	{
		std::vector<FGraphNodeHandle> matchingElements = algorithm->sourceExamples(elementObject.nodeHandle());

		Element::eachNode(matchingElements, algorithm->exemplar().graph, [&](FGraphNode& matchingNode) {
			FTransform transform(matchingNode.orientation, matchingNode.position, FVector(debug_matchingPointSize));

			if (showMatchingExamplesBounds.IsInside(matchingNode.position))
			{
				ismc.AddInstanceWorldSpace(transform);
			}
		});
	}
}