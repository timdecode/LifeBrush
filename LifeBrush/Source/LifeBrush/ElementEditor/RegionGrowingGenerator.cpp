//
//  LifeBrush
//
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ShipEditorSimulation/MeshSimulation.h"

#include "RegionGrowingGenerator.h"
#include "Simulation/FlexElements.h"
#include "Simulation/Aggregates.h"

void URegionGrowingGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager_hack)
{
	Super::attach(context, simulationManager_hack);

	if (!context || !simulationManager_hack)
		return;

	_generationWorker.resize(1);

	_outputContext = context;
	
	// make sure the kd-tree is updated, other generators might not do this
	_outputContext->domain.rebalance();

	if( !_outputSimulationManager)
		_outputSimulationManager = NewObject<UGraphSimulationManager>(this);

	_outputSimulationManager->init(context->graph(), *simulationManager_hack->actor());
	_outputSimulationManager->registerSimulation<UMeshSimulation>();
	_outputSimulationManager->attachSimulations();

	if (!_algorithm)
		_algorithm = std::make_unique<Algorithm_RegionGrowing>(*_outputContext);

	if( !_exemplarSimulationManager)
		_exemplarSimulationManager = NewObject<UGraphSimulationManager>(this);

	_exemplarSimulationManager->init(_algorithm->exemplarGraph(), *simulationManager_hack->actor());
	_exemplarSimulationManager->attachSimulations();


	loadParameters();

	_algorithm->reloadContext();

	syncExemplarFromElementActors();
}

void URegionGrowingGenerator::detach()
{
	Super::detach();

	// sync the algorithm and execute any work queued for the main thread
	_generationWorker.push([](int threadID) mutable {}).get();

	if (_generationWorkerMainThreadWork)
	{
		_generationWorkerMainThreadWork();
		_generationWorkerMainThreadWork = nullptr;
	}


	// hack set the static positions
	auto& staticPositions = _outputContext->graph().componentStorage<FStaticPositionObject>();

	for (auto& sp : staticPositions)
	{
		FGraphNode& node = _outputContext->graph().node(sp.nodeHandle());

		sp.position = node.position;
	}

	_outputSimulationManager->detachSimulations();
	_exemplarSimulationManager->detachSimulations();
}

void URegionGrowingGenerator::tick(float deltaT)
{
	if (!_outputContext || !_outputSimulationManager || !_algorithm)
		return;

	if (!_generationWorkFinished)
		return;

	_generationWorkFinished = false;

	if (_generationWorkerMainThreadWork)
	{
		_generationWorkerMainThreadWork();
		_generationWorkerMainThreadWork = nullptr;
	}

	_outputContext->domain.graph.pushTransactionContext();
	{
		// tick the simulation manager
		if (pauseSynthesis)
		{
			_outputSimulationManager->tick_paused(deltaT);
		}
		else
		{
			_outputSimulationManager->tick(deltaT);
		}

		// the exemplar is always paused
		_exemplarSimulationManager->tick_paused(deltaT);

	}
	_outputContext->domain.graph.popTransactionContext();

	// short circuit
	if (pauseSynthesis)
		return;
	
	Algorithm::AABB limits;
	limits.aabb.min() = eigen(_outputContext->limits.Min);
	limits.aabb.max() = eigen(_outputContext->limits.Max);

	// copy the toRemove
	auto toRemove = _toRemove;
	_toRemove.clear();

	// push the context here, we pop it in the _generationWorkerMainThreadWork lambda (on the main thread)
	_outputContext->domain.graph.pushTransactionContext();
	_generationWorker.push([this, limits, toRemove](int threadID) mutable
	{

		// remove points
		AlgorithmResult removalResult;
		{
			std::set<FGraphNodeHandle> removedElements;
			for (auto& p : toRemove)
			{
				auto erasedThisTime = _algorithm->eraseAt(p);

				removedElements.insert(erasedThisTime.begin(), erasedThisTime.end());
			}

			removalResult.removed.insert(removalResult.removed.end(), removedElements.begin(), removedElements.end());
		}

		// and generate
		AlgorithmResult generationResult;

		std::vector<PositionFace> dummyNoPositions;

		_algorithm->generate(dummyNoPositions, 10000.0f, limits);

		generationResult.append(removalResult);

		// timings
		// we do the calculations here, as the roundSummary duration could change in the meantime
		RoundRecord round = _algorithm->roundSummary().back();

		double roundDuration = round.duration();
		double instantaneousGenerationRate = double(round.elementsGenerated) / round.duration();
		double totalDuration = _algorithm->roundSummary().duration();

		double totalGenerationRate = double(round.elementsTotal) / totalDuration;

		// put some results back on the main thread
		_generationWorkerMainThreadWork = 
			[this,
			generationResult,
			round,
			instantaneousGenerationRate,
			totalGenerationRate,
			totalDuration,
			roundDuration] ()
			{
				if (round.elementsGenerated)
					UE_LOG(LogTemp, Warning, TEXT("instantaneous generation rate %f"), instantaneousGenerationRate);

				if (round.elementsTotal)
					UE_LOG(LogTemp, Warning, TEXT("total generation rate %f total time %f"), totalGenerationRate, totalDuration);

				// HERE! pop the context, on the main thread, yay
				_outputContext->domain.graph.popTransactionContext();
			};

		_generationWorkFinished = true;
	});
}

void URegionGrowingGenerator::start()
{
	pauseSynthesis = false;
}

void URegionGrowingGenerator::stop()
{
	// sync the algorithm and execute any work queued for the main thread
	_generationWorker.push([](int threadID) mutable {}).get();

	pauseSynthesis = true;
}

void URegionGrowingGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex)
{
	// hack for now
// we should do something more intelligent, like checking if the brush stroke has completed growing
	_generationWorker.push([this](int threadID) mutable
	{
		_algorithm->clearBrushPoints();
	});

	_toRemove.clear();
}

void URegionGrowingGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	Eigen::Vector3f eigen_point = eigen(point);

	_generationWorker.push([this, eigen_point, radius, surfaceIndex](int threadID) mutable
	{
		PositionRadiusFace brushPoint(eigen_point, radius, surfaceIndex);
		_algorithm->addBrushPoint(brushPoint);
		_algorithm->setCurrentBrushPoint(brushPoint);
	});
}

void URegionGrowingGenerator::endBrushPath()
{
	_generationWorker.push([this](int threadID) mutable
	{
		_algorithm->clearCurrentBrushPoint();
	});
}

void URegionGrowingGenerator::beginEraserPath(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	_generationWorker.push([this](int threadID) mutable
	{
		_algorithm->clearBrushPoints();
	});

	_toRemove.clear();
}

void URegionGrowingGenerator::eraseInRadiusAt(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	// space the points apart
	if (_toRemove.size())
	{
		auto& last = _toRemove.back();

		const float distSqrd = (eigen(point) - last.position).squaredNorm();
		if (distSqrd < last.radius * 0.2)
			return;
	}

	_toRemove.emplace_back(eigen(point), brushSize, surfaceIndex);
}

void URegionGrowingGenerator::endEraserPath()
{

}







void URegionGrowingGenerator::setSelection(std::vector<AElementActor*> elementActors)
{
	if (!_algorithm)
		return;

	// join the worker
	_generationWorker.push([](int threadID) mutable {}).get();

	// build the selection from the _actorToElement map
	std::vector<FGraphNodeHandle> selection;
	selection.reserve(elementActors.size());

	for (auto actor : elementActors)
	{
		FGraphNodeHandle handle = _actorToElement[actor];

		if (handle)
			selection.push_back(handle);
	}

	// resync
	syncExemplarFromElementActors();

	_algorithm->updateExampleSelection(selection);
}

void URegionGrowingGenerator::setGenerationMode(EGenerationMode generationMode)
{
	_generationMode = generationMode;

	_generationWorker.push([this, generationMode](int threadID) mutable
	{
		_algorithm->generationMode = generationMode;
	});
}

EGenerationMode URegionGrowingGenerator::generationMode()
{
	return _generationMode;
}

void URegionGrowingGenerator::syncExemplarFromElementActors()
{
	if (!_algorithm || !exemplar)
		return;

	TArray<USceneComponent*> exemplarSceneComponents;
	exemplar->GetRootComponent()->GetChildrenComponents(false, exemplarSceneComponents);

	auto& exemplarDomain = _algorithm->exemplar();

	exemplarDomain.graph.beginTransaction();

	std::map<URegionGrowingComponent::ElementTypeDescription, int, URegionGrowingComponent::InferredElementTypeDescriptionComparator> actorToTypes;

	// infer our element types
	int nextIndex = 0;
	for (auto sceneComponent : exemplarSceneComponents)
	{
		AElementActor * elementActor = Cast<AElementActor>(sceneComponent->GetOwner());

		if (!elementActor)
			continue;

		URegionGrowingComponent::ElementTypeDescription typeDescription(elementActor);

		if (actorToTypes.find(typeDescription) == actorToTypes.end())
		{
			actorToTypes[typeDescription] = nextIndex;
			nextIndex++;
		}
	}

	// build the elements
	TMap<FGraphNodeHandle, AElementActor*> elementToActor;

	for (auto sceneComponent : exemplarSceneComponents)
	{
		AElementActor * elementActor = Cast<AElementActor>(sceneComponent->GetOwner());

		if (!elementActor)
			continue;

		int elementID = actorToTypes[elementActor];

		FVector position = sceneComponent->GetComponentLocation();
		FQuat orientation = sceneComponent->GetComponentRotation().Quaternion();
		float scale = sceneComponent->GetComponentScale().X;

		// we already have this guy, update him
		if (_actorToElement.find(elementActor) != _actorToElement.end())
		{
			auto h = _actorToElement[elementActor];

			ElementTuple element(h, exemplarDomain.graph);

			element.node.position = position;
			element.node.orientation = orientation;
			element.node.scale = scale;

			// first nuke the components that aren't an element
			static const ComponentType ElementType = componentType<FElementObject>();

			auto components = element.node.components;
			for (ComponentType component : components)
			{
				if( component == ElementType ) continue;

				element.node.removeComponent(exemplarDomain.graph, component);
			}

			elementActor->writeToElement(element, exemplarDomain.graph);

			elementToActor.Add(h, elementActor);
		}
		// create a new element
		else
		{
			auto h = exemplarDomain.insert(position, orientation, scale);
			ElementTuple element(h, exemplarDomain.graph);

			element.element.type = elementID;

			_actorToElement[elementActor] = h;

			UStaticMeshComponent * mesh = (UStaticMeshComponent*)elementActor->GetComponentByClass(UStaticMeshComponent::StaticClass());

			if (mesh == nullptr || mesh->GetStaticMesh() == nullptr)
				continue;

			elementActor->writeToElement(element, exemplarDomain.graph);

			elementToActor.Add(h, elementActor);
		}
	}

	// nuke element that are gone from the exemplar
	{
		// get rid of the nodes
		for (FGraphNode& node : exemplarDomain.graph.allNodes)
		{
			if (!node.isValid()) continue;

			if (!elementToActor.Contains(node.handle()))
			{
				exemplarDomain.graph.removeNode(node.handle());
			}
		}

		// get rid of old references in the table
		for (auto& pair : _actorToElement)
		{
			if (!elementToActor.Contains(pair.second))
				_actorToElement.erase(pair.first);
		}
	}

	exemplarDomain.graph.endTransaction();
}


















































FString URegionGrowingGenerator::parametersString()
{
	FString params = "";

	params += FString::Printf(TEXT("perElementParameters %d\n"), perElementParameters);

	params += FString::Printf(TEXT("disableOptimization %d\n"), disableOptimization);
	params += FString::Printf(TEXT("optimizationRounds %d\n"), optimizationRounds);
	params += FString::Printf(TEXT("useGlobalOptimization %d\n"), useGlobalOptimization);

	params += FString::Printf(TEXT("enableRoundSummaries %d\n"), enableRoundSummaries);
	params += FString::Printf(TEXT("calculate_kCoherenceEnergy %d\n"), calculate_kCoherenceEnergy);
	params += FString::Printf(TEXT("calculate_bruteForceEnergy %d\n"), calculate_bruteForceEnergy);

	{
		const UEnum * GenerationMode = FindObject<UEnum>(ANY_PACKAGE, TEXT("EGenerationMode"));
		int32 int32_generationMode = (int32)_generationMode;

		params += FString::Printf(TEXT("generationMode: %s\n"), *(GenerationMode->GetDisplayNameTextByValue(int32_generationMode).ToString()));
	}

	params += FString::Printf(TEXT("sketchBasedForceTangentPlane %d\n"), sketchBasedForceTangentPlane);

	params += FString::Printf(TEXT("generationParameters\n"));
	params += _optimizationParametersString(generationParameters);

	params += FString::Printf(TEXT("optimizationParameters\n"));
	params += _optimizationParametersString(optimizationParameters);

	params += FString::Printf(TEXT("useTypeVoting %d\n"), useTypeVoting);
	params += FString::Printf(TEXT("removeOverlaps %d\n"), removeOverlaps);
	params += FString::Printf(TEXT("enableReassignment %d\n"), enableReassignment);
	params += FString::Printf(TEXT("generationInnerRadius %f\n"), generationRadius);
	params += FString::Printf(TEXT("kCoherence %d\n"), kCoherence);
	params += FString::Printf(TEXT("minAssignmentDistance %f\n"), minAssignmentDistance);
	params += FString::Printf(TEXT("freespaceRadius %f\n"), freespaceRadius);
	params += FString::Printf(TEXT("typeCost %f\n"), typeCost);
	params += FString::Printf(TEXT("gradientTerm %f\n"), gradientTerm);
	params += FString::Printf(TEXT("relaxation %f\n"), relaxation);
	params += FString::Printf(TEXT("flipSurfaceNormals %d\n"), flipSurfaceNormals);
	params += FString::Printf(TEXT("ignoreBadSuggestions %d\n"), ignoreBadSuggestions);
	params += FString::Printf(TEXT("ignoreBadSuggestionsDistanceFactor %f\n"), ignoreBadSuggestionsDistanceFactor);

	params += FString::Printf(TEXT("kCoherenceClusteringPenalty %f\n"), kCoherenceClusteringPenalty);
	params += FString::Printf(TEXT("kCoherenceClusteringPenaltyRadius %f\n"), kCoherenceClusteringPenaltyRadius);

	params += FString::Printf(TEXT("sourceHistogramRadius %f\n"), sourceHistogramRadius);
	params += FString::Printf(TEXT("sourceHistogramWeight %f\n"), sourceHistogramWeight);
	params += FString::Printf(TEXT("samplingDistanceWeight %f\n"), samplingDistanceWeight);

	params += FString::Printf(TEXT("brushSize %f\n"), brushSize);

	return params;
}

void URegionGrowingGenerator::loadParameters()
{
	if (!_algorithm)
		return;

	_algorithm->perElementParameters = perElementParameters;

	_algorithm->enableRoundSummaries = enableRoundSummaries;
	_algorithm->calculate_kCoherenceEnergy = calculate_kCoherenceEnergy;
	_algorithm->calculate_bruteForceEnergy = calculate_bruteForceEnergy;
	_algorithm->cellExtents = cinpactCellSize;
	_algorithm->useCinpactEnergy = useCinpactEnergy;

	_algorithm->generationMode = _generationMode;
	_algorithm->sketchBasedForceTangentPlane = sketchBasedForceTangentPlane;
	_algorithm->generationParameters = generationParameters;
	_algorithm->optimizationParameters = optimizationParameters;
	_algorithm->useTypeVoting = useTypeVoting;
	_algorithm->removeOverlaps = removeOverlaps;
	_algorithm->enableReassignment = enableReassignment;
	_algorithm->generationInnerRadius = generationRadius;
	_algorithm->frozenElementRadius = frozenElementRadius;
	_algorithm->kCoherence = kCoherence;
	_algorithm->minAssignmentDistance = minAssignmentDistance;
	_algorithm->freespaceRadius = freespaceRadius;
	_algorithm->typeCost = typeCost;
	_algorithm->gradientTerm = gradientTerm;
	_algorithm->relaxation = relaxation;
	_algorithm->flipSurfaceNormals = flipSurfaceNormals;
	_algorithm->ignoreBadSuggestions = ignoreBadSuggestions;
	_algorithm->ignoreBadSuggestionsDistanceFactor = ignoreBadSuggestionsDistanceFactor;

	_algorithm->kCoherenceClusteringPenalty = kCoherenceClusteringPenalty;
	_algorithm->kCoherenceClusteringPenaltyRadius = kCoherenceClusteringPenaltyRadius;

	_algorithm->enableVolumeSurfaceInteraction = enableVolumeSurfaceInteraction;

	_algorithm->disableOptimization = disableOptimization;
	_algorithm->optimizationRounds = optimizationRounds;
	_algorithm->useGlobalOptimization = useGlobalOptimization;

	_algorithm->forceRotation = false;
	_algorithm->forcedRotation = FQuat::Identity;


	_algorithm->sourceHistogramRadius = sourceHistogramRadius;
	_algorithm->sourceHistogramWeight = sourceHistogramWeight;
	_algorithm->samplingDistanceWeight = samplingDistanceWeight;
}

FString URegionGrowingGenerator::_optimizationParametersString(FNeighbourhoodParameters& params)
{
	FString result = "";

	result += FString::Printf(TEXT("   radius %f\n"), params.radius);
	result += FString::Printf(TEXT("   kNearest %f\n"), params.kNearest);

	return result;
}
