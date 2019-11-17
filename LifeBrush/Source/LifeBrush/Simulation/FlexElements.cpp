// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "Visualization/Timeline.h"
#include "SimulationSnapshotActor.h"
#include "MolecularLego/MolecularLego_Relaxation.h"

#include "Utility.h"
#include "SimulationSnapshotActor.h"
#include "VolumeComponent.h"
#include "ElementEditor/SwarmGenerator.h"




UGraphComponent::UGraphComponent()
{
	bWantsInitializeComponent = true;
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.TickGroup = TG_EndPhysics;
	bIsActive = true;

	simulationManager = CreateDefaultSubobject<UGraphSimulationManager>(TEXT("GraphSimulationManager"));
}

void UGraphComponent::InitializeComponent()
{
	Super::InitializeComponent();

	initGraph();
}

void UGraphComponent::initGraph()
{
	if (_didInitGraph)
		return;

	graph.init();

	simulationManager->init(graph, *GetOwner(), true);

	_didInitGraph = true;
}







void UFlexSimulationComponent::playSimulation()
{
	if (!flexSimulation() || !_didInitOnce )
		return;

	flexSimulation()->updateFlexState();

	flexSimulation()->begin();

	flexSimulation()->play();
}

void UFlexSimulationComponent::pauseSimulation()
{
	if (!flexSimulation() )
		return;

	flexSimulation()->pause();
}

void UFlexSimulationComponent::initWithCameraComponent(UCameraComponent * camera)
{
	if (_didInitOnce)
		return;

	FTransform meshInterfaceToWorld = FTransform::Identity;

	init(meshInterfaceToWorld, camera);

	_didInitOnce = true;
}

bool UFlexSimulationComponent::isSimulationPlaying()
{
	return flexSimulation()->isPlaying();
}

UFlexSimulationComponent::UFlexSimulationComponent()
{
	bWantsInitializeComponent = true;
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.TickGroup = TG_EndPhysics;
	bIsActive = true;
	
	bTickInEditor = true;
}

void UFlexSimulationComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	_initContext();

	if (drawLimits && !_limitsBoxComponent )
		_drawLimits();

	_updateLimits();

	if( !_flexSimulation)
		initFlexSimulationObject();

	if( tickFlexSimulation )
		_flexSimulation->tick(timeStep);
}

void UFlexSimulationComponent::BeginDestroy()
{
	if( _flexSimulation )
		_flexSimulation->uninit();

	Super::BeginDestroy();
}

void UFlexSimulationComponent::InitializeComponent()
{
	Super::InitializeComponent();

	if (!rulesActor)
	{
		if (softMLRulesActor.IsValid() && softMLRulesActor.TryLoad())
			rulesActor = Cast<AActor>(softMLRulesActor.ResolveObject());
	}

	if (!swarmGrammarRulesActor)
	{
		if (softSwarmGrammarRulesActor.IsValid() && softSwarmGrammarRulesActor.TryLoad())
			swarmGrammarRulesActor = Cast<AActor>(softSwarmGrammarRulesActor.ResolveObject());
	}

	initFlexSimulationObject();
}

void UFlexSimulationComponent::BeginPlay()
{
	Super::BeginPlay();

	_cacheMeshesForMeshInterface();
}

void UFlexSimulationComponent::initFlexSimulationObject()
{
	if (_flexSimulation)
		return;

	initGraph();

	_updateLimits();
	_initContext();

	_flexSimulation = std::make_unique<FFlexSimulation>(
		graph,
		*simulationManager,
		flexParams
	);

	_readMLRules();
	_readSGRules();

	_flexSimulation->setInstanceManagerBounds(limits);
	_flexSimulation->initSimulationManager(GetOwner());
}

void UFlexSimulationComponent::_initContext()
{
	if (_context)
		return;

	_context = std::make_unique<SynthesisContext>(graph);
	_context->limits = limits;
}

void UFlexSimulationComponent::_readMLRules()
{
	if (!rulesActor) return;

	UMLElementSimulation * ruleSimulation = _flexSimulation->simulationManager.registerSimulation<UMLElementSimulation>();

	ruleSimulation->ruleGraph.init();
	ruleSimulation->ruleGraph.clear();

	TArray<USceneComponent*> childSceneComponents;
	rulesActor->GetRootComponent()->GetChildrenComponents(false, childSceneComponents);

	for (USceneComponent * child : childSceneComponents)
	{
		AMLRuleActor * ruleActor = Cast<AMLRuleActor>(child->GetOwner());

		if( !ruleActor ) continue;

		ruleActor->writeToGraph(ruleSimulation->ruleGraph);
	}
}

void UFlexSimulationComponent::_readSGRules()
{
	if (!swarmGrammarRulesActor) return;

	USwarmSimulation * ruleSimulation = _flexSimulation->simulationManager.registerSimulation<USwarmSimulation>();

	ruleSimulation->ruleGraph.init();
	ruleSimulation->ruleGraph.clear();

	TArray<USceneComponent*> childSceneComponents;
	swarmGrammarRulesActor->GetRootComponent()->GetChildrenComponents(false, childSceneComponents);

	FGraph * ruleGraph = &ruleSimulation->ruleGraph;

	for (USceneComponent * child : childSceneComponents)
	{
		AElementActor * ruleActor = Cast<AElementActor>(child->GetOwner());

		if (!ruleActor) continue;

		FGraphNodeHandle newHandle = FGraphNodeHandle(ruleGraph->addNode(
			ruleActor->GetActorLocation(),
			ruleActor->GetActorQuat(),
			ruleActor->GetActorScale3D().GetMax()));

		ruleGraph->node(newHandle).addComponent<FElementObject>(*ruleGraph);

		ElementTuple tuple = ElementTuple(newHandle, *ruleGraph);

		ruleActor->writeToElement(tuple, ruleSimulation->ruleGraph);
	}
}

auto UFlexSimulationComponent::init(
	FTransform meshInterfaceToWorld,
	UCameraComponent * camera) -> void
{
	_initMeshInterface(_flexSimulation->meshInterface());

	simulationManager->camera = camera;
	simulationManager->init(graph, *GetOwner(), true);
}

void UFlexSimulationComponent::_initMeshInterface(std::shared_ptr<tcodsMeshInterface> meshInterfacePtr)
{
	auto& meshInterface = *meshInterfacePtr.get();

	for (UStaticMeshComponent * meshComponent : _meshesForMeshInterface)
	{
		// build the mesh interface up
		if (meshComponent && meshComponent->GetStaticMesh())
		{
			auto meshes = tcodsUStaticMeshBuilder::buildMeshes(*meshComponent->GetStaticMesh(), meshComponent->GetComponentToWorld());

			for (auto& mesh : meshes)
			{
				auto section = meshInterface.addMesh(std::move(mesh));
			}
		}
	}

	meshInterface.setWorldLimits(_flexSimulation->instanceManagerBounds());
	meshInterface.rebuild();

	_context->meshInterface = meshInterfacePtr;
}

void UFlexSimulationComponent::_cacheMeshesForMeshInterface()
{
	AActor * actor = GetOwner();
	USceneComponent * root = actor->GetRootComponent();

	TArray<UStaticMeshComponent*> components;
	actor->GetComponents<UStaticMeshComponent>(components);

	// create our RMC
	_meshInterfaceRMC = NewObject<URuntimeMeshComponent>(actor);
	_meshInterfaceRMC->AttachToComponent(root, FAttachmentTransformRules::KeepRelativeTransform);
	_meshInterfaceRMC->RegisterComponent();
	actor->AddInstanceComponent(_meshInterfaceRMC);

	// find the static meshes
	for (USceneComponent * sceneComponent : components)
	{
		if (Cast<UInstancedStaticMeshComponent>(sceneComponent) != nullptr) continue;

		UStaticMeshComponent * meshComponent = Cast<UStaticMeshComponent>(sceneComponent);

		if( meshComponent && meshComponent->GetStaticMesh() )
			_meshesForMeshInterface.Add(meshComponent);
	}
}

void UFlexSimulationComponent::updateMeshInterface(UChunkedVolumeComponent * chunkVolume)
{
	checkfSlow(chunkVolume, TEXT("nullptr chunkVolume"));
	checkfSlow(context->meshInterface, TEXT("The context is missing a meshInterface"));

	auto& meshInterface = *_flexSimulation->meshInterface().get();

	ChunkGrid<float>& grid = chunkVolume->grid();

	// clear our old sections
	for (auto section : _chunkSections)
	{
		meshInterface.removeMesh(section);

		if( _meshInterfaceRMC) _meshInterfaceRMC->ClearMeshSection(section);
	}

	_chunkSections.clear();

	// add new sections
	auto meshes = tcodsChunkedGridMeshBuilder::buildMeshes(grid, chunkVolume->GetOwner()->ActorToWorld(), chunkVolume->isoLevel);

	for (auto& mesh : meshes)
	{
		auto section = meshInterface.addMesh(std::move(mesh));

		_chunkSections.push_back(section);
	}

	meshInterface.rebuild();

	FBox limits = _flexSimulation->instanceManagerBounds();

	if (_meshInterfaceRMC)
	{
		const FTransform transform = _meshInterfaceRMC->GetComponentTransform().Inverse();

		for (auto section : _chunkSections)
		{
			meshInterface.buildRuntimeMeshSection(_meshInterfaceRMC, meshInterface.mesh(section), section, transform, true, true, limits);
			_meshInterfaceRMC->SetMaterial(section, meshInterfaceMaterial);
		}
	}
}



ASimulationSnapshotActor* UFlexSimulationComponent::createGraphicalSnapshotActor(UWorld * world)
{
	ASimulationSnapshotActor * newActor = world->SpawnActor<ASimulationSnapshotActor>();

	USceneComponent * newScene = NewObject<USceneComponent>(newActor);
	{
		AActor * originalActor = GetOwner();

		USceneComponent * originalScene = originalActor->GetRootComponent();

		FTransform newTransform = originalScene->GetComponentTransform();

		newActor->SetRootComponent(newScene);
		newActor->AddInstanceComponent(newScene);

		newScene->SetWorldTransform(newTransform);
	}

	for (UObjectSimulation * simulation : _flexSimulation->simulationManager.simulations)
	{
		simulation->snapshotToActor(newActor);
	}

	return newActor;
}

AActor * UFlexSimulationComponent::snapshotSimulationStateToWorld(UWorld * world)
{
	AActor * newActor = world->SpawnActor<AActor>();

	AActor * originalActor = GetOwner();

	// create scene root
	USceneComponent * newScene = NewObject<USceneComponent>(newActor);
	{
		USceneComponent * originalScene = originalActor->GetRootComponent();

		FTransform newTransform = originalScene->GetComponentTransform();

		newActor->SetRootComponent(newScene);
		newActor->AddInstanceComponent(newScene);

		newScene->SetWorldTransform(newTransform);
	}

	// duplicate the simulation
	{
		UFlexSimulationComponent * newElements = DuplicateObject(this, newActor);

		newActor->AddInstanceComponent(newElements);
	}

	TArray<UStaticMeshComponent*> components;
	originalActor->GetComponents<UStaticMeshComponent>(components);

	// duplicate the static meshes
	{
		for (USceneComponent * sceneComponent : components)
		{
			if( Cast<UInstancedStaticMeshComponent>(sceneComponent) ) continue;

			UStaticMeshComponent * meshComponent = Cast<UStaticMeshComponent>(sceneComponent);

			if( !meshComponent) continue;

			UStaticMeshComponent * newMeshComponent = DuplicateObject(meshComponent, newActor);

			newActor->AddInstanceComponent(newMeshComponent);
		}
	}

	// duplicate the _meshInterfaceRMC
	if (_meshInterfaceRMC && _meshInterfaceRMC->GetNumSections() > 0)
	{
		// export the visible mesh
		URuntimeMeshComponent * newRMC = Utility::duplicateRuntimeMeshComponentToActor(_meshInterfaceRMC, newActor);

		if (newRMC)
			newActor->AddInstanceComponent(newRMC);
	}

	{
		// export the collision mesh too (we need this for tcods)
		URuntimeMeshComponent * tcodsMeshInterfaceRMC = NewObject<URuntimeMeshComponent>(newActor);
		tcodsMeshInterfaceRMC->SetCollisionEnabled(ECollisionEnabled::NoCollision);
		tcodsMeshInterfaceRMC->AttachToComponent(newActor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
		newActor->AddInstanceComponent(tcodsMeshInterfaceRMC);
		tcodsMeshInterfaceRMC->RegisterComponent();

		URuntimeMeshComponent * tcodsVisibleRMC = NewObject<URuntimeMeshComponent>(newActor);
		tcodsVisibleRMC->SetCollisionEnabled(ECollisionEnabled::NoCollision);
		tcodsVisibleRMC->AttachToComponent(newActor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
		newActor->AddInstanceComponent(tcodsVisibleRMC);
		tcodsVisibleRMC->RegisterComponent();

		const FTransform transform = tcodsMeshInterfaceRMC->GetComponentTransform().Inverse();

		auto& meshInterface = *_flexSimulation->meshInterface().get();

		for (auto section : _chunkSections)
		{
			meshInterface.buildRuntimeMeshSection(tcodsMeshInterfaceRMC, meshInterface.mesh(section), section, transform, false /* clip to limits */, false /* double sided */);		
			//tcodsMeshInterfaceRMC->SetMaterial(section, meshInterfaceMaterial);

			meshInterface.buildRuntimeMeshSection(tcodsVisibleRMC, meshInterface.mesh(section), section, transform, true /* clip to limits */, true /* double sided */);
			tcodsVisibleRMC->SetMaterial(section, meshInterfaceMaterial);
		}
	}
	


	return newActor;
}

void UFlexSimulationComponent::_updateLimits()
{
	if (!_limitsBoxComponent || drawLimits == false)
		return;

	FTransform relativeTransform = _limitsBoxComponent->GetRelativeTransform();

	FVector centre = relativeTransform.GetLocation();
	FVector extents = _limitsBoxComponent->GetUnscaledBoxExtent();

	extents *= relativeTransform.GetScale3D();

	limits.Min = centre - extents;
	limits.Max = centre + extents;
}

void UFlexSimulationComponent::_updateDrawLimits()
{
	if (drawLimits)
	{
		_drawLimits();
		_updateLimits();
	}
	else
		_hideLimits();
}

void UFlexSimulationComponent::_drawLimits()
{
	if (_limitsBoxComponent == nullptr)
	{
		_limitsBoxComponent = NewObject<UBoxComponent>(GetOwner());

		AActor * actor = GetOwner();

		USceneComponent * root = actor->GetRootComponent();

		_limitsBoxComponent->SetCollisionEnabled(ECollisionEnabled::NoCollision);
		_limitsBoxComponent->SetCollisionProfileName(TEXT("NoCollision"));

		_limitsBoxComponent->AttachToComponent(root, FAttachmentTransformRules::KeepRelativeTransform);
		actor->AddInstanceComponent(_limitsBoxComponent);
		_limitsBoxComponent->RegisterComponent();


	}

	if (!_limitsBoxComponent->IsVisible())
		_limitsBoxComponent->SetVisibility(true);

	_limitsBoxComponent->SetWorldScale3D(FVector(1.0f));
	_limitsBoxComponent->SetWorldRotation(FQuat::Identity, false, nullptr, ETeleportType::TeleportPhysics);
	_limitsBoxComponent->SetWorldLocation(limits.GetCenter());
	_limitsBoxComponent->SetBoxExtent(limits.GetExtent());
}

void UFlexSimulationComponent::_hideLimits()
{
	if (_limitsBoxComponent)
		_limitsBoxComponent->SetVisibility(false);
}















auto FFlexSimulation::tick(float deltaT) -> void
{
	if (!_playing)
	{
		_tickGraph_paused(deltaT);

		return;
	}

	// do we need to load state into flex?
	bool needsFirstInit = false;

	if (_needsFlexReset)
	{
		if (_didInitFlex)
			uninit();

		_initFlex();
		_didInitFlex = true;

		_needsFlexReset = false;

		needsFirstInit = true;

		_hackInitChannels();
	}

	if (!_didInitFlex)
		return;

	// map
	FVector4 * particles = (FVector4*)NvFlexMap(_particleBuffer, eNvFlexMapWait);
	FVector * velocities = (FVector*)NvFlexMap(_velocityBuffer, eNvFlexMapWait);
	int * phases = (int*)NvFlexMap(_phaseBuffer, eNvFlexMapWait);
	int * active = (int*)NvFlexMap(_activeBuffer, eNvFlexMapWait);

	NvFlexCollisionGeometry * geometry = (NvFlexCollisionGeometry*)NvFlexMap(_geometryBuffer, eNvFlexMapWait);
	FVector4 * shapePositions = (FVector4*)NvFlexMap(_shapePositionsBuffer, eNvFlexMapWait);
	FQuat * shapeRotations = (FQuat*)NvFlexMap(_shapeRotationsBuffer, eNvFlexMapWait);
	int * shapeFlags = (int*)NvFlexMap(_shapeFlagsBuffer, eNvFlexMapWait);

	_springs->map();

	_neighborsIndicesBuffer->map();
	_neighborsCountsBuffer->map();
	_neighborsAPIToInternal->map();
	_neighborsInternalToAPI->map();

	bool debug_didInitSprings = false;
	if (needsFirstInit)
	{
		// make sure we have created our buffers
		_springs->init();

		debug_didInitSprings = true;
	}



	// spawn
	if (!_didSpawnShapesAndMesh)
	{
		_spawnShapes(geometry, shapePositions, shapeRotations, shapeFlags);
		_spawnMesh(geometry, shapePositions, shapeRotations, shapeFlags);

		_didSpawnShapesAndMesh = true;
	}

	// springs and positions
	if (needsFirstInit)
	{
		_writeFlexState(particles, velocities, phases, active);
		_writeSpringState();
	}

	// update the push-interaction sphere
	(FVector4&)shapePositions[0] = _spherePosition;
	shapePositions[0].W = 1.0f;
	geometry[0].sphere.radius = _sphereRadius;

	// Read Flex state into the simulation
	_readFlexState(particles, velocities, phases);

	// cache this, as we need to access it after the step too (during which it can be modified)

	// rigids
	if (needsFirstInit)
	{
		_rigids->init();
		_rigids->mapAll();
		_loadRigids();
	}
	else
	{
		if (_flexRigidsDirty)
		{
			_rigids->mapAll();
			_loadRigids();
		}
		else
			_rigids->mapRotationsTranslations();

		_readRigidRotations();
	}

	// simulate
	_preTick(deltaT);

	{
		// do tick work
		for (auto& work : _tickWork)
			work();

		_tickWork.clear();

		// do flex tick work
		for (auto& work : _flexTickWork)
		{
			work(deltaT,
				*_neighborsIndicesBuffer.get(),
				*_neighborsCountsBuffer.get(),
				*_neighborsAPIToInternal.get(),
				*_neighborsInternalToAPI.get(),
				_solverDescription.maxParticles
			);
		}

		// graph tick
		_tickGraph(deltaT);

		_flexTick(deltaT);

		// simulate particle rotation outside of Flex
		_integrateRotations(deltaT);
	}
	//timeline->endFrame();
	_postTick(deltaT);

	// Write simulation state back into Flex
	_writeFlexState(particles, velocities, phases, active);

	_writeSpringState();



	bool needsLoadRigids = _flexRigidsDirty;

	if (needsLoadRigids)
	{
		_rigids->init();
		_rigids->mapRemaining();
		_loadRigids();
	}

	// unmap
	NvFlexUnmap(_particleBuffer); particles = nullptr;
	NvFlexUnmap(_velocityBuffer); velocities = nullptr;
	NvFlexUnmap(_phaseBuffer); phases = nullptr;
	NvFlexUnmap(_activeBuffer); active = nullptr;

	NvFlexUnmap(_geometryBuffer);
	NvFlexUnmap(_shapePositionsBuffer);
	NvFlexUnmap(_shapeRotationsBuffer);
	NvFlexUnmap(_shapeFlagsBuffer);

	_springs->unmap();

	if (needsFirstInit)
		_rigids->unmapAll();
	else if (needsLoadRigids)
		_rigids->unmapAll();
	else
		_rigids->unmapRotationsTranslations();

	_neighborsIndicesBuffer->unmap();
	_neighborsCountsBuffer->unmap();
	_neighborsAPIToInternal->unmap();
	_neighborsInternalToAPI->unmap();

	// write to device (async)
	NvFlexSetParticles(_solver, _particleBuffer, nullptr);
	NvFlexSetVelocities(_solver, _velocityBuffer, nullptr);
	NvFlexSetPhases(_solver, _phaseBuffer, nullptr);
	NvFlexSetActive(_solver, _activeBuffer, nullptr);
	NvFlexSetActiveCount(_solver, nParticles);

	if (!_springs->indices.empty())
	{
		NvFlexSetSprings(_solver, _springs->indices.buffer, _springs->lengths.buffer, _springs->coefficients.buffer, _springs->coefficients.size());
	}

	if (needsFirstInit || needsLoadRigids)
	{
		NvFlexSetRigids(
			_solver,
			_rigids->offsets.buffer,
			_rigids->indices.buffer,
			_rigids->localRestPositions.buffer,
			_rigids->localRestNormals.buffer,
			_rigids->stiffness.buffer,
			_rigids->plasticThresholds.buffer,
			_rigids->plasticCreeps.buffer,
			_rigids->bodyRotations.buffer,
			_rigids->bodyTranslations.buffer,
			_rigids->offsets.size() - 1, // numRigids
			_rigids->indices.size() // numIndices
		);
	}


	NvFlexSetShapes(_solver, _geometryBuffer, _shapePositionsBuffer, _shapeRotationsBuffer, nullptr, nullptr, _shapeFlagsBuffer, nGeometries);

	// tick
	NvFlexSetParams(_solver, &_params);

	NvFlexUpdateSolver(_solver, deltaT, 1, false);

	// read back (async)
	NvFlexGetParticles(_solver, _particleBuffer, nullptr);
	NvFlexGetVelocities(_solver, _velocityBuffer, nullptr);
	NvFlexGetPhases(_solver, _phaseBuffer, nullptr);
	NvFlexGetActive(_solver, _activeBuffer, nullptr);
	NvFlexGetNeighbors(
		_solver,
		_neighborsIndicesBuffer->buffer,
		_neighborsCountsBuffer->buffer,
		_neighborsAPIToInternal->buffer,
		_neighborsInternalToAPI->buffer
	);
	NvFlexGetRigids(
		_solver,
		nullptr,
		nullptr,
		nullptr,
		nullptr,
		nullptr,
		nullptr,
		nullptr,
		_rigids->bodyRotations.buffer,
		_rigids->bodyTranslations.buffer
	);
}

auto FFlexSimulation::uninit() -> void
{
	// we already released
	if(!_solver)
		return;

	NvFlexFreeBuffer( _particleBuffer ); _particleBuffer = nullptr;
	NvFlexFreeBuffer( _velocityBuffer ); _velocityBuffer = nullptr;
	NvFlexFreeBuffer( _phaseBuffer ); _phaseBuffer = nullptr;
	NvFlexFreeBuffer( _activeBuffer ); _activeBuffer = nullptr;

	NvFlexFreeBuffer( _geometryBuffer ); _geometryBuffer = nullptr;
	NvFlexFreeBuffer( _shapePositionsBuffer ); _shapePositionsBuffer = nullptr;
	NvFlexFreeBuffer( _shapeRotationsBuffer ); _shapeRotationsBuffer = nullptr;
	NvFlexFreeBuffer( _shapeFlagsBuffer ); _shapeFlagsBuffer = nullptr;

	_mesh0_indices = nullptr;
	_mesh0_positions = nullptr;

	_springs = nullptr;

	_rigids = nullptr;

	_neighborsIndicesBuffer = nullptr;
	_neighborsCountsBuffer = nullptr;
	_neighborsAPIToInternal = nullptr;
	_neighborsInternalToAPI = nullptr;

	NvFlexDestroySolver( _solver ); _solver = nullptr;
	NvFlexShutdown( _library ); _library = nullptr;

	_didSpawnShapesAndMesh = false;

	graphSimulation.removeComponentListener<FFlexParticleObject>((ComponentListener*)this);
	graphSimulation.removeComponentListener<FFlexRigidBodyObject>((ComponentListener*)this);
	graphSimulation.removeEdgeObjectListener<FFlexRigidBodyConnection>((EdgeObjectListener*)this);
}

void FFlexSimulation::_preTick(float DeltaTime)
{
	for (auto simulation : simulationManager.simulations)
	{
		IFlexGraphSimulation * flexSim = Cast<IFlexGraphSimulation>(simulation);

		if (flexSim)
			flexSim->preTick(DeltaTime);
	}
}

void FFlexSimulation::_tickGraph( float DeltaTime )
{
	if(_didInitSImulationManager)
	{
		simulationManager.tick(DeltaTime);
	}
}

void FFlexSimulation::_tickGraph_paused(float DeltaTime)
{
	if (_didInitSImulationManager)
	{
		simulationManager.tick_paused(DeltaTime);
	}
}

void FFlexSimulation::_flexTick(float DeltaTime)
{
	for (auto simulation : simulationManager.simulations)
	{
		IFlexGraphSimulation * flexSim = Cast<IFlexGraphSimulation>(simulation);

		if (flexSim)
		{
			flexSim->flexTick(
				DeltaTime,
				*_neighborsIndicesBuffer.get(),
				*_neighborsCountsBuffer.get(),
				*_neighborsAPIToInternal.get(),
				*_neighborsInternalToAPI.get(),
				_solverDescription.maxParticles
			);
		}
	}


}

void FFlexSimulation::_postTick(float DeltaTime)
{
	for (auto simulation : simulationManager.simulations)
	{
		IFlexGraphSimulation * flexSim = Cast<IFlexGraphSimulation>(simulation);

		if (flexSim)
			flexSim->postTick(DeltaTime);
	}
}

auto FFlexSimulation::begin() -> void
{
	simulationManager.begin();
}

auto FFlexSimulation::_initFlex() -> void
{
	/** Descriptor used to initialize Flex

	struct NvFlexInitDesc
	{
		int deviceIndex;				//!< The GPU device index that should be used, if there is already a CUDA context on the calling thread then this parameter will be ignored and the active CUDA context used. Otherwise a new context will be created using the suggested device ordinal.
		bool enableExtensions;			//!< Enable or disable NVIDIA/AMD extensions in DirectX, can lead to improved performance.
		void* renderDevice;				//!< Direct3D device to use for simulation, if none is specified a new device and context will be created.
		void* renderContext;			//!< Direct3D context that the app is using for rendering. In DirectX 12 this should be a ID3D12CommandQueue pointer.
		void* computeContext;           //!< Direct3D context to use for simulation, if none is specified a new context will be created, in DirectX 12 this should be a pointer to the ID3D12CommandQueue where compute operations will take place.
		bool runOnRenderContext;		//!< If true, run Flex on D3D11 render context, or D3D12 direct queue. If false, run on a D3D12 compute queue, or vendor specific D3D11 compute queue, allowing compute and graphics to run in parallel on some GPUs.

		NvFlexComputeType computeType;	//!< Set to eNvFlexD3D11 if DirectX 11 should be used, eNvFlexD3D12 for DirectX 12, this must match the libraries used to link the application
	};
	*/

	_initDesc.deviceIndex = 1;
	_initDesc.enableExtensions;
	_initDesc.runOnRenderContext = false;

	// _library = NvFlexInit(NV_FLEX_VERSION, 0, &_initDesc);

	_library = NvFlexInit(NV_FLEX_VERSION, 0, 0);


	NvFlexSetSolverDescDefaults(&_solverDescription);
	_solverDescription.maxParticles = _maxParticles;
	_solverDescription.maxNeighborsPerParticle = _maxNeighbors;
	_solverDescription.maxDiffuseParticles = 0;

	_solver = NvFlexCreateSolver(_library, &_solverDescription);

	_particleBuffer = NvFlexAllocBuffer(_library, _maxParticles, sizeof(FVector4), eNvFlexBufferHost);
	_velocityBuffer = NvFlexAllocBuffer(_library, _maxParticles, sizeof(FVector4), eNvFlexBufferHost);
	_phaseBuffer = NvFlexAllocBuffer(_library, _maxParticles, sizeof(int), eNvFlexBufferHost);
	_activeBuffer = NvFlexAllocBuffer(_library, _maxParticles, sizeof(int), eNvFlexBufferHost);

	_geometryBuffer = NvFlexAllocBuffer(_library, nGeometries, sizeof(NvFlexCollisionGeometry), eNvFlexBufferHost);
	_shapePositionsBuffer = NvFlexAllocBuffer(_library, nGeometries, sizeof(FVector4), eNvFlexBufferHost);
	_shapeRotationsBuffer = NvFlexAllocBuffer(_library, nGeometries, sizeof(FVector4), eNvFlexBufferHost);
	_shapeFlagsBuffer = NvFlexAllocBuffer(_library, nGeometries, sizeof(int), eNvFlexBufferHost);

	_springs = std::make_unique<Springs>(_library);

	_rigids = std::make_unique<Rigids>(_library);

	{
		const int neighboursSize = _solverDescription.maxParticles * _solverDescription.maxNeighborsPerParticle;

		_neighborsIndicesBuffer = std::make_unique<NvFlexVector<int>>(_library, neighboursSize);
		_neighborsCountsBuffer =  std::make_unique<NvFlexVector<int>>(_library, _solverDescription.maxParticles);
		_neighborsAPIToInternal = std::make_unique<NvFlexVector<int>>(_library, neighboursSize);
		_neighborsInternalToAPI = std::make_unique<NvFlexVector<int>>(_library, neighboursSize);
	}

	if (_meshInterface)
		_loadTcodsMesh(*_meshInterface.get());


	_initParams();
	NvFlexSetParams(_solver, &_params);
}


auto FFlexSimulation::clear() -> void
{
	graphSimulation.clear();

	if( _meshInterface )
		_needsFlexReset = true;
	
	_didInitFlex = false;
	_didSpawnShapesAndMesh = false;
}

auto FFlexSimulation::pause() -> void
{
	_playing = false;


}

auto FFlexSimulation::play() -> void
{
	_playing = true;
}

auto FFlexSimulation::updateSphereWorldSpace( FVector position, float radius ) -> void
{
	FTransform transform = owner->GetRootComponent()->GetComponentTransform();

	_spherePosition = transform.InverseTransformPosition( position );
	_sphereRadius = radius / transform.GetScale3D().GetMax();
}

auto FFlexSimulation::updateFlexState() -> void
{
	graphSimulation.beginTransaction();

	FTransform toWorld = FTransform::Identity;

	auto& velocities = graphSimulation.componentStorage<FVelocityGraphObject>();
	auto& elements = graphSimulation.componentStorage<FElementObject>();
	auto& particles = graphSimulation.componentStorage<FFlexParticleObject>();

	for (FGraphNode& node : graphSimulation.allNodes)
	{
		if( !node.isValid() )
			continue;


		// velocity
		if (!node.hasComponent<FVelocityGraphObject>())
		{
			node.addComponent<FVelocityGraphObject>(graphSimulation);
		}

		bool onSurface = false;

		// constrain objects to the surface if they have a face
		FElementObject * elementObject = elements.componentPtrForNode(node.handle());
		if (elementObject && elementObject->surfaceIndex.isOnSurface())
		{
			FSurfaceBoundGraphObject& boundObject = node.addComponent<FSurfaceBoundGraphObject>(graphSimulation);

			auto rotationAndNormal = _meshInterface->rotationAndNormalAtIndex(elementObject->surfaceIndex);

			boundObject.lastSurfaceRotation = rotationAndNormal.first;
			boundObject.surfaceIndex = elementObject->surfaceIndex;

			onSurface = true;
		}

		// set particle
		if (FFlexParticleObject * particle = particles.componentPtrForNode(node.handle()))
		{
			// don't interact with the surface
			particle->channel = onSurface ? eNvFlexPhaseShapeChannel1 : eNvFlexPhaseShapeChannel0;
		}

	}

	graphSimulation.endTransaction();

	URandomWalkSimulation * randomWalkSim = simulationManager.simulation<URandomWalkSimulation>();

	randomWalkSim->meshInterface = _meshInterface;

	_needsFlexReset = true;
}

auto FFlexSimulation::exportElementDomain() -> FGraphSnapshot
{
	FGraphSnapshot result;

	result.snapshot(graphSimulation);

	// transform to world space
	FTransform toWorld = owner->GetRootComponent()->GetComponentTransform();

	for (FGraphNode& node : result.graph.allNodes)
	{
		if( !node.isValid() )
			continue;

		node.position = toWorld.TransformPosition(node.position);
		node.orientation = toWorld.TransformRotation(node.orientation);
		node.scale = toWorld.GetScale3D().X * node.scale;
	}

	return result;
}



auto FFlexSimulation::setInstanceManagerBounds( FBox bounds ) -> void
{
	for(int i = 0; i < 3; ++i)
	{
		FVector n_up( 0.0f );
		n_up[i] = 1.0f;

		FVector n_down( 0.0f );
		n_down[i] = -1.0f;

		float w_up = -(bounds.Min[i] * n_up[i]);
		float w_down = -(bounds.Max[i] * n_down[i]);

		flexParams.planes[i * 2 + 0] = FVector4( n_up.X, n_up.Y, n_up.Z, w_up );
		flexParams.planes[i * 2 + 1] = FVector4( n_down.X, n_down.Y, n_down.Z, w_down );
	}

	flexParams.numPlanes = 6;

	_bounds = bounds;
}

auto FFlexSimulation::addTickWork(std::function<void()> work) -> void
{
	_tickWork.push_back(work);
}

auto FFlexSimulation::addFlexTickWork(FlexTickWork_t work) -> void
{
	_flexTickWork.push_back(work);
}

auto FFlexSimulation::initSimulationManager(AActor * owner_in) -> void
{
	_FFlexParticleObjectType = FGraphObject::componentType( FFlexParticleObject::StaticStruct() );
	_FFlexRigidObjectType = FGraphObject::componentType(FFlexRigidBodyObject::StaticStruct());

	owner = owner_in;

	graphSimulation.tickCount = 0;

	// the order is important
	// disconnect all listeners
	// instantiate the simulations
	// init the graph
	// init the simulation manager

	graphSimulation.resetListeners();


	FTransform toWorld = owner->GetRootComponent()->GetComponentTransform();

	UVisualization_AgentPathLines * agentPathLines = simulationManager.registerSimulation<UVisualization_AgentPathLines>();
	agentPathLines->toWorld = toWorld;

	URandomWalkSimulation * randomWalkSim = simulationManager.registerSimulation<URandomWalkSimulation>();
	randomWalkSim->toWorld = toWorld;
	randomWalkSim->meshInterface = _meshInterface;

	// setup listeners
	graphSimulation.addComponentListener<FFlexParticleObject>((ComponentListener*) this );
	graphSimulation.addComponentListener<FFlexRigidBodyObject>((ComponentListener*)this);

	graphSimulation.addEdgeObjectListener<FFlexConnection>(this);
	graphSimulation.addEdgeObjectListener<FFlexRigidBodyConnection>(this);

	simulationManager.attachAllSimulations();

	_didInitSImulationManager = true;

	_needsFlexReset = true;
	_didInitFlex = false;
	_didSpawnShapesAndMesh = false;
}

auto FFlexSimulation::setCamera(UCameraComponent * camera) -> void
{
	_camera = camera;

	if (UVisualization_AgentPathLines * agentPathLines = simulationManager.simulation<UVisualization_AgentPathLines>() )
	{
		agentPathLines->camera = camera;
	}
}

auto FFlexSimulation::_initParams() -> void
{
	// sim params
	_params = flexParams.toNvFlexParams();

	if(_params.solidRestDistance == 0.0f)
		_params.solidRestDistance = _params.radius;

	// if fluid present then we assume solid particles have the same radius
	if(_params.fluidRestDistance > 0.0f)
		_params.solidRestDistance = _params.fluidRestDistance;

	// set collision distance automatically based on rest distance if not already set
	if(_params.collisionDistance == 0.0f)
		_params.collisionDistance = std::max( _params.solidRestDistance, _params.fluidRestDistance )*0.5f;

	// default particle friction to 10% of shape friction
	if(_params.particleFriction == 0.0f)
		_params.particleFriction = _params.dynamicFriction*0.1f;

	// add a margin for detecting contacts between particles and shapes
	if(_params.shapeCollisionMargin == 0.0f)
		_params.shapeCollisionMargin = _params.collisionDistance*0.5f;
}

auto FFlexSimulation::_instanceCount() -> size_t
{
	size_t sum = 0;
	for(auto& pair : _instanceManager.elementTypesToInstancedMeshes)
	{
		sum += pair.second->GetInstanceCount();
	}

	return sum;
}

void FFlexSimulation::_hackInitChannels()
{
	auto& flexParticles = graphSimulation.componentStorage<FFlexParticleObject>();

	for (FFlexParticleObject& particleObject : flexParticles)
	{
		FGraphNode& node = graphSimulation.node(particleObject.nodeIndex);

		bool onSurface = node.hasComponent<FSurfaceBoundGraphObject>();

		// don't interact with the surface
		particleObject.channel = onSurface ? eNvFlexPhaseShapeChannel1 : eNvFlexPhaseShapeChannel0;
	}
}

void FFlexSimulation::_writeSpringState()
{
	auto connections = graphSimulation.edgeView<FFlexConnection>();

	int n = 0;
	connections.each([&](FFlexConnection& connection, FGraphEdge& edge)
	{
		n++;
	});

	auto& indices = _springs->indices;
	auto& lengths = _springs->lengths;
	auto& coefficients = _springs->coefficients;

	indices.resize(n * 2);
	lengths.resize(n);
	coefficients.resize(n);


	int i = 0;
	connections.each([&](FFlexConnection& connection, FGraphEdge& edge)
	{
		indices[(i * 2) + 0] = edge.a;
		indices[(i * 2) + 1] = edge.b;

		lengths[i] = connection.length;
		coefficients[i] = connection.coefficient;

		++i;
	});
}

void FFlexSimulation::_loadRigids()
{
	TypedImmutableEdgeStorage<FFlexRigidBodyConnection> rigidEdges = graphSimulation.edgeStorage<FFlexRigidBodyConnection>();

	auto& velocities = graphSimulation.componentStorage<FVelocityGraphObject>();
	auto& particles = graphSimulation.componentStorage<FFlexParticleObject>();
	auto& baseRotations = graphSimulation.componentStorage<FFlexBaseRotation>();
	auto& members = graphSimulation.componentStorage<FFlexRigidMember>();

	auto rigidN = graphSimulation.componentStorage<FFlexRigidBodyObject>().size();

	_rigids->offsets.reserve(rigidEdges.validSize() + 1);
	_rigids->indices.reserve(rigidEdges.validSize());
	_rigids->stiffness.reserve(rigidN);
	_rigids->plasticThresholds.reserve(rigidN);
	_rigids->plasticCreeps.reserve(rigidN);
	_rigids->bodyRotations.reserve(rigidN);
	_rigids->bodyTranslations.reserve(rigidN);
	_rigids->localRestPositions.reserve(rigidN);
	_rigids->localRestNormals.reserve(rigidN);


	// Flex: first entry has to be 0, for some reason
	if( _rigids->indices.empty() )
		_rigids->offsets.push_back(0);

	int rigidIndex = 0;
	graphSimulation.view<FFlexRigidBodyObject>().each([&](FGraphNode& node, FFlexRigidBodyObject& rigid) 
	{
		rigid.flexRigidIndex = rigidIndex;

		node.position = rigid.calculateCenterOfMass(graphSimulation);

		checkfSlow(!node.hasComponent<FFlexParticleObject>(), TEXT("The node with the FFlexRigidBodyObject must not have a particle."));


		FVelocityGraphObject nodeVelocity = node.component<FVelocityGraphObject>(graphSimulation);

		FQuat inverseRotation = node.orientation.Inverse();

		// - set the rest positions and normals
		// - find the centre of mass
		for (auto edgeIndex : node.edges)
		{
			// do we have a rigid connection?
			FGraphEdgeHandle edgeHandle(edgeIndex);
			FFlexRigidBodyConnection * rigidConnection = rigidEdges.objectPtr(edgeHandle);

			if (!rigidConnection) continue;

			FGraphNodeHandle other = edgeHandle(graphSimulation).other(node.handle());

			FFlexBaseRotation& otherBaseRotation = other(graphSimulation).addComponent<FFlexBaseRotation>(graphSimulation);
			otherBaseRotation.rotation = inverseRotation * other(graphSimulation).orientation;
			otherBaseRotation.flexRigidBodyIndex = rigidIndex;

			// load the particle at the other end
			if (FFlexParticleObject * otherParticle = particles.componentPtrForNode(other))
			{
				if (FVelocityGraphObject * otherVelocity = velocities.componentPtrForNode(other))
				{
					otherVelocity->linearVelocity = nodeVelocity.linearVelocity;
					otherVelocity->angularVelocity = nodeVelocity.angularVelocity;
				}

				otherParticle->group = node.handle().index;
				otherParticle->selfCollide = false;

				_rigids->indices.push_back(other.index);

				const FVector restPosition = inverseRotation.RotateVector(other(graphSimulation).position - node.position);
				// The restNormal.w is supposed to be the signed rest distance inside of the shape.
				// HACK, we'll just say it is at the surface.
				const FVector4 restNormal = FVector4(restPosition.GetSafeNormal(), 0.0f);

				_rigids->localRestPositions.push_back(restPosition);
				_rigids->localRestNormals.push_back(restNormal);
			}

			// load the member at the other end
			if (FFlexRigidMember * member = members.componentPtrForNode(other))
			{
				member->relativeOffset = inverseRotation.RotateVector(other(graphSimulation).position - node.position);
				member->relativeOrientation = inverseRotation * other(graphSimulation).orientation;
				member->flexRigidIndex = rigidIndex;
			}
		}

		_rigids->stiffness.push_back(rigid.stiffness);
		_rigids->plasticThresholds.push_back(rigid.plasticDeformationThreshold);
		_rigids->plasticCreeps.push_back(rigid.plasticDeformationCreep);

		_rigids->bodyRotations.push_back(node.orientation);
		_rigids->bodyTranslations.push_back(node.position);

		// Flex uses offsets like this:
		// body_i start = offsets[i]
		// body_i end = offsets[i + 1] // exclusive
		// Therefore, offsets must be n + 1 in size, since the last rigid body, needs an end index
		// This is also why we do the push_back here.
		_rigids->offsets.push_back(_rigids->indices.size());

		rigidIndex++;
	});

	_flexRigidsDirty = false;
}

void FFlexSimulation::_readRigidRotations()
{
	// read the base rotations to the rigid body sub-particles
	graphSimulation.each_node_object<FFlexBaseRotation>([&](FGraphNode& node, FFlexBaseRotation& baseRotation) {
		FQuat rotation = _rigids->bodyRotations[baseRotation.flexRigidBodyIndex];

		node.orientation = rotation * baseRotation.rotation;
	});

	// read our flex rigid body object position and orientation
	graphSimulation.each_node_object<FFlexRigidBodyObject>([&](FGraphNode& node, FFlexRigidBodyObject& rigid) {
		node.position = _rigids->bodyTranslations[rigid.flexRigidIndex];
		node.orientation = _rigids->bodyRotations[rigid.flexRigidIndex];
	});

	// read FFlexRigidMembers
	graphSimulation.each_node_object<FFlexRigidMember>([&](FGraphNode& node, FFlexRigidMember& member) {
		FQuat rotation = _rigids->bodyRotations[member.flexRigidIndex];
		FVector position = _rigids->bodyTranslations[member.flexRigidIndex];

		node.position = rotation.RotateVector(member.relativeOffset) + position;
		node.orientation = rotation * member.relativeOrientation;
	});
}

void FFlexSimulation::_spawnShapes( NvFlexCollisionGeometry * geometry, FVector4 * shapePositions, FQuat * shapeRotations, int * shapeFlags )
{
	shapeFlags[0] = NvFlexMakeShapeFlagsWithChannels( eNvFlexShapeSphere, true, eNvFlexPhaseShapeChannel0 );
	geometry[0].sphere.radius = _sphereRadius;

	(FVector4&)shapePositions[0] = _spherePosition; 	
	shapePositions[0].W = 1.0f;

	shapeRotations[0] = FQuat::Identity;
}

void FFlexSimulation::_spawnMesh( NvFlexCollisionGeometry * geometry, FVector4 * shapePositions, FQuat * shapeRotations, int * shapeFlags )
{
	const size_t meshIndex = 1;

	if (!_mesh0_positions || _mesh0_positions->empty())
		return;

	_mesh0_Id = NvFlexCreateTriangleMesh( _library );

	NvFlexUpdateTriangleMesh(
		_library,
		_mesh0_Id,
		_mesh0_positions->buffer,
		_mesh0_indices->buffer,
		_mesh0_positions->size(),
		_mesh0_indices->size() / 3, // num triangles?
		&_mesh0_bounds.Min[0],
		&_mesh0_bounds.Max[0]
	);

	// shapeFlags[meshIndex] = NvFlexMakeShapeFlagsWithChannels(eNvFlexShapeTriangleMesh, false, eNvFlexPhaseShapeChannel0 | eNvFlexPhaseShapeChannel1);
	// only interact with channel 0
	// channel 1 is for surface bound objects
	shapeFlags[meshIndex] = NvFlexMakeShapeFlagsWithChannels( eNvFlexShapeTriangleMesh, false, eNvFlexPhaseShapeChannel0 );

	geometry[meshIndex].triMesh.mesh = _mesh0_Id;

	geometry[meshIndex].triMesh.scale[0] = 1.0f;
	geometry[meshIndex].triMesh.scale[1] = 1.0f;
	geometry[meshIndex].triMesh.scale[2] = 1.0f;

	shapePositions[meshIndex] = FVector4( 0.0f );
	shapeRotations[meshIndex] = FQuat::Identity;


}

auto FFlexSimulation::_readFlexState( FVector4 * positions, FVector * velocities, int * phases ) -> void
{
	// read positions
	auto& particles = graphSimulation.componentStorage<FFlexParticleObject>();

	for(FFlexParticleObject& particle : particles)
	{
		if (!particle.isValid()) continue;

		FGraphNode& node = graphSimulation.node( particle.nodeIndex );

		FVector4& p4 = positions[node.id];
		FVector p( p4.X, p4.Y, p4.Z );

		particle.inverseMass = p4.W;

		node.position = p;
	}

	// read velocities
	auto& velocityObjects = graphSimulation.componentStorage<FVelocityGraphObject>();

	for(FVelocityGraphObject& velocityObject : velocityObjects)
	{
		if (!velocityObject.isValid()) continue;

		auto id = velocityObject.nodeIndex;

		FVector v = velocities[id];

		velocityObject.linearVelocity = v;
	}
}




auto FFlexSimulation::_writeFlexState( FVector4 * positions, FVector * velocities, int * phases, int * active ) -> void
{
	// write positions
	auto& particles = graphSimulation.componentStorage<FFlexParticleObject>();

	nParticles = 0;

	int activeIndex = 0;

	for(FFlexParticleObject& particle : particles)
	{
		if (!particle.isValid()) continue;

		FGraphNode& node = graphSimulation.node( particle.nodeIndex );

		positions[node.id] = node.position;
		positions[node.id].W = particle.inverseMass;

		phases[node.id] = NvFlexMakePhaseWithChannels(
			particle.group, 
			(particle.selfCollide ? eNvFlexPhaseSelfCollide : 0) | (particle.isFlud ? eNvFlexPhaseFluid : 0), 
			particle.channel);

		active[activeIndex] = particle.nodeIndex;

		++activeIndex;
		++nParticles;
	}

	// write velocities
	auto& velocityObjects = graphSimulation.componentStorage<FVelocityGraphObject>();

	for(FVelocityGraphObject& velocityObject : velocityObjects)
	{
		if( !velocityObject.isValid() ) continue;

		auto nodeID = velocityObject.nodeIndex;

		velocities[nodeID] = velocityObject.linearVelocity;
	}
}

auto FFlexSimulation::_integrateRotations(float deltaT) -> void
{
	auto& velocityObjects = graphSimulation.componentStorage<FVelocityGraphObject>();

	for(FVelocityGraphObject& velocityObject : velocityObjects)
	{
		if (!velocityObject.isValid()) continue;

		// https://gamedev.stackexchange.com/questions/108920/applying-angular-velocity-to-quaternion
		// q = q_0 + dt/2 * w * q_0

		FGraphNodeHandle node = FGraphNodeHandle(velocityObject.nodeIndex);

		FVector angularVelocity = velocityObject.angularVelocity;

		FQuat w = FQuat( angularVelocity.X, angularVelocity.Y, angularVelocity.Z, 0.0f );

		FQuat& rotation = node(graphSimulation).orientation;

		w *= 0.5f * deltaT;
		w *= rotation;

		rotation = rotation + (w);

		rotation.Normalize();
	}
}

auto FFlexSimulation::_loadTcodsMesh(tcodsMeshInterface& tcodsMesh) -> void
{
	FTransform worldToUsTransform = owner->GetRootComponent()->GetComponentTransform().Inverse();

	FTransform toUsTransform = _meshToWorld * worldToUsTransform;

	// data to be consumed by flex
	std::vector<FVector4> flexPositions;
	std::vector<int> flexIndices;

	_mesh0_bounds = FBox(EForceInit::ForceInitToZero);

	

	for (auto& ptr : tcodsMesh._sections)
	{
		using namespace tcods;

		const int32 sectionOffset = flexPositions.size();

		Mesh& mesh = *ptr.second.get();

		for (auto& vertex : mesh.vertices)
		{
			FVector v = unreal(vertex.position);

			v = toUsTransform.TransformPosition(v);

			_mesh0_bounds += v;

			flexPositions.push_back(FVector4(v, 0.0f));
		}

		for (auto& face : mesh.faces)
		{
			if (face.isBoundary())
				continue;

			auto he = face.he;

			// add a face
			std::array<int32, 3> triangle;
			int i = 0;
			do 
			{
				triangle[i++] = he->from->index + sectionOffset;

				he = he->next;
			} while ( he != face.he);

			assert(i == 3);

			flexIndices.push_back(triangle[0]);
			flexIndices.push_back(triangle[1]);
			flexIndices.push_back(triangle[2]);

			// add the reverse face (double sided)
			flexIndices.push_back(triangle[0]);
			flexIndices.push_back(triangle[2]);
			flexIndices.push_back(triangle[1]);
		}
	}

	_mesh0_positions = std::make_unique< NvFlexVector<FVector4> >(_library, &flexPositions[0], flexPositions.size());
	_mesh0_indices = std::make_unique< NvFlexVector<int> >(_library, &flexIndices[0], flexIndices.size());
}

auto FFlexSimulation::_loadMesh( UStaticMeshComponent * meshComponent ) -> void
{
	FTransform worldToUsTransform = owner->GetRootComponent()->GetComponentTransform().Inverse();


	FTransform toUsTransform = meshComponent->GetComponentTransform();
	toUsTransform = toUsTransform * worldToUsTransform;

	UStaticMesh& uMesh = *meshComponent->GetStaticMesh();

	if(uMesh.RenderData == nullptr) return;
	if(uMesh.RenderData->LODResources.Num() == 0) return;

	// data to be consumed by flex
	std::vector<FVector4> flexPositions;
	std::vector<int> flexIndices;

	FStaticMeshLODResources& resource = uMesh.RenderData->LODResources[meshComponent->PreviousLODLevel];

	FPositionVertexBuffer& positionBuffer = resource.PositionVertexBuffer;
	FRawStaticIndexBuffer& indexBuffer = resource.IndexBuffer;
	FStaticMeshVertexBuffer& vertexBuffer = resource.VertexBuffer;

	uint32 vCount = positionBuffer.GetNumVertices();

	// Unreal stores duplicate vertices, we need to resolve this with a map
	std::unordered_map<FVector, int32> vertexLookup;

	std::vector<uint32> vertexIndices;
	vertexIndices.resize( vCount );

	_mesh0_bounds = FBox();

	for(uint32 i = 0; i < vCount; ++i)
	{
		const FVector& v = positionBuffer.VertexPosition( i );


		auto found = vertexLookup.find( v );
		if(found == vertexLookup.end())
		{
			vertexLookup[v] = flexPositions.size();

			vertexIndices[i] = flexPositions.size();

			FVector vTransformed = toUsTransform.TransformPosition( v );

			_mesh0_bounds += vTransformed;

			flexPositions.push_back( FVector4( vTransformed, 0.0f) );
		}
		else
			vertexIndices[i] = found->second;
	}

	FIndexArrayView indexArrayView = indexBuffer.GetArrayView();
	int32 iCount = indexArrayView.Num();
	for(int32 i = 0; i + 2 < iCount; i += 3)
	{
		int i0 = vertexIndices[indexArrayView[i + 0]];
		int i1 = vertexIndices[indexArrayView[i + 1]];
		int i2 = vertexIndices[indexArrayView[i + 2]];

		// double-sided geometry
		flexIndices.push_back( i0 );
		flexIndices.push_back( i1 );
		flexIndices.push_back( i2 );

		flexIndices.push_back( i0 );
		flexIndices.push_back( i2 );
		flexIndices.push_back( i1 );
	}

	_mesh0_positions = std::make_unique< NvFlexVector<FVector4> >( _library, &flexPositions[0], flexPositions.size() );
	_mesh0_indices = std::make_unique< NvFlexVector<int> >( _library, &flexIndices[0], flexIndices.size() );
}



void FFlexSimulation::componentAdded( FGraphNodeHandle handle, ComponentType type )
{
	const static ComponentType FFlexRigidMemberType = componentType<FFlexRigidMember>();
	const static ComponentType FFlexRigidBodyObjectType = componentType<FFlexRigidBodyObject>();


	if (type == _FFlexParticleObjectType)
	{
		FGraphNode& node = graphSimulation.node(handle);

		FFlexParticleObject& particle = node.component<FFlexParticleObject>(graphSimulation);

		bool onSurface = node.hasComponent<FSurfaceBoundGraphObject>();

		// don't interact with the surface
		particle.channel = onSurface ? eNvFlexPhaseShapeChannel1 : eNvFlexPhaseShapeChannel0;
	}
	else if (type == FFlexRigidBodyObjectType)
	{
		_flexRigidsDirty = true;
	}
	else if (type == FFlexRigidMemberType)
	{
		_flexRigidsDirty = true;
	}
}

void FFlexSimulation::componentRemoved( FGraphNodeHandle node, ComponentType type )
{
	const static ComponentType FFlexRigidBodyObjectType = componentType<FFlexRigidBodyObject>();

	if (type != FFlexRigidBodyObjectType)
	{
	}
	else if (type == _FFlexRigidObjectType)
	{
		flexRigidBodyObjectRemoved(node);
	}
}


void FFlexSimulation::flexRigidBodyObjectRemoved(FGraphNodeHandle handle)
{
	_flexRigidsDirty = true;

	// remove the attached connections
	FGraphNode& node = handle(graphSimulation);

	auto rigidBodyEdgeStorage = graphSimulation.edgeStorage<FFlexRigidBodyConnection>();

	for (auto ei : node.edges)
	{
		auto edgeObject = rigidBodyEdgeStorage.objectPtr(FGraphEdgeHandle(ei));

		if( !edgeObject || !rigidBodyEdgeStorage.isValid(*edgeObject) ) continue;

		graphSimulation.removeEdgeObject<FFlexRigidBodyConnection>(FGraphEdgeHandle(ei));
	}
}

void FFlexSimulation::connectionAdded( int32 edgeIndex, FGraphEdge& edge, ComponentType type )
{
}

void FFlexSimulation::connectionRemoved( int32 edgeIndex, FGraphEdge& oldEdge, ComponentType type )
{
}

void FFlexSimulation::edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type)
{
	if (type == edgeType<FFlexConnection>())
	{

	}
	else if (type == edgeType<FFlexRigidBodyConnection>())
	{
		_flexRigidsDirty = true;
	}
}

void FFlexSimulation::edgeObjectRemoved(FGraphEdgeHandle handle, EdgeObjectType type)
{
	if (type == edgeType<FFlexConnection>())
	{

	}
	else if (type == edgeType<FFlexRigidBodyConnection>())
	{
		_flexRigidsDirty = true;
	}
}




































void URandomWalkSimulation::attach()
{

}

void URandomWalkSimulation::tick( float deltaT )
{
	auto& walkers = graph->componentStorage<FRandomWalkGraphObject>();

	FRandomStream rand;
	rand.GenerateNewSeed();

	auto& surfaceDwellers = graph->componentStorage<FSurfaceBoundGraphObject>();
	auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();

	for(FSurfaceBoundGraphObject& dweller : surfaceDwellers)
	{
		FGraphNode& node = graph->node( dweller.nodeIndex );

		FVector p = toWorld.TransformPosition( node.position );

		auto nearest = meshInterface->nearestPointOnMesh( p );

		if( !nearest.surfaceIndex.isOnSurface() )
			continue;

		auto rotationAndNormal = meshInterface->rotationAndNormalAtIndex( nearest.surfaceIndex );

		FQuat surfaceRotation = toWorld.InverseTransformRotation(rotationAndNormal.first);

		FQuat rotation = surfaceRotation * dweller.lastSurfaceRotation.Inverse();

		node.position = toWorld.InverseTransformPosition( nearest.point );
		//node.orientation = rotation * node.orientation;

		dweller.lastSurfaceRotation = surfaceRotation;
		dweller.surfaceIndex = nearest.surfaceIndex;

		// project linear velocity back onto the surface
		if(node.hasComponent<FVelocityGraphObject>())
		{
			FVelocityGraphObject& velocityObject = node.component<FVelocityGraphObject>(*graph);

			FVector& v = velocityObject.linearVelocity;

			// rotate the vector onto the surface
			v = rotation.RotateVector( v );

			// trip any residual vertical motion

			FVector normal = toWorld.InverseTransformVectorNoScale(rotationAndNormal.second);

			v = v - (FVector::DotProduct( v, normal ) * normal);
		}
	}

	for(FRandomWalkGraphObject& walker : walkers)
	{
		walker.timeLeft -= deltaT;

		// times' up, choose a random direction
		if(walker.timeLeft >= 0.0f )
			continue;

		FGraphNode& node = graph->node( walker.nodeIndex );

		if(!node.hasComponent<FVelocityGraphObject>())
			continue;

		FVelocityGraphObject& velocityObject = node.component<FVelocityGraphObject>( *graph );

		walker.timeLeft = rand.GetFraction() * 0.3f;


		// surface walker
		if(node.hasComponent<FSurfaceBoundGraphObject>())
		{
			FSurfaceBoundGraphObject& boundObject = node.component<FSurfaceBoundGraphObject>( *graph );

			// calculate random unit circle vector
			float t = rand.GetFraction() * 2.0f * M_PI;

			FVector v = FVector::ZeroVector;

			FMath::SinCos( &v.Y, &v.X, t );

			v = v *  (walker.baseVelocity + walker.maxVelocityOffset * rand.GetFraction());

			// now, we need to rotate it from surface space into sim space
			auto rotationAndNormal = meshInterface->rotationAndNormalAtIndex( boundObject.surfaceIndex );

			FQuat r = toWorld.InverseTransformRotation(rotationAndNormal.first);

			v = r.RotateVector( v );

			velocityObject.linearVelocity = v;
		}
		else // space walker
		{
			FVector v = rand.GetUnitVector() * (walker.baseVelocity + walker.maxVelocityOffset * rand.GetFraction());

			velocityObject.linearVelocity += v;
		}
	}

	auto& spinners = graph->componentStorage<FSpinnerGraphObject>();

	for(FSpinnerGraphObject& spinner : spinners)
	{
		FGraphNodeHandle node( spinner.nodeIndex );

		if(!node( graph ).hasComponent<FVelocityGraphObject>())
			continue;

 		FVelocityGraphObject& velocityObject = node( graph ).component<FVelocityGraphObject>( *graph );


		FQuat& orientation = node( graph ).orientation;

		FVector up = FVector::UpVector;

		up = orientation.RotateVector( up );

		up *= spinner.angularVelocityMagnitude;

		velocityObject.angularVelocity = up;
	}
}




void UATPSynthaseSimulation::attach()
{
	_g = a_max * (r_min * r_min);
}

void UATPSynthaseSimulation::tick( float deltaT )
{

}

void UATPSynthaseSimulation::flexTick( 
	float deltaT, 
	NvFlexVector<int>& neighbourIndices, 
	NvFlexVector<int>& neighbourCounts, 
	NvFlexVector<int>& apiToInternal, 
	NvFlexVector<int>& internalToAPI, 
	int maxParticles 
)
{
	auto& synthases = graph->componentStorage<FATPSynthaseGraphObject>();
	auto& hydrogens = graph->componentStorage<FHydrogenGraphObject>();
	auto& particles = graph->componentStorage<FFlexParticleObject>();

	const int stride = maxParticles;

	// If we have hydrogen below us, spin and teleport it across the membrane
	// If we are spinning and have ADP above us, destroy it and produce ADP

	// we'll be adding and removing components, so start a transaction
	graph->beginTransaction();

	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	for(FATPSynthaseGraphObject& synthase : synthases)
	{
		if( !synthase.isValid() ) continue;

		auto nodeIndex = synthase.nodeIndex;
		FGraphNode& atpSynthaseNode = graph->node(nodeIndex);

		if(!atpSynthaseNode.hasComponent<FFlexParticleObject>())
			continue;

		FVector atpSynthaseUp = FVector::UpVector;
		atpSynthaseUp = atpSynthaseNode.orientation.RotateVector( atpSynthaseUp );

		const float interactionRadiusSqrd = synthase.hyrogenInteractionRadius * synthase.hyrogenInteractionRadius;

		int nodeIndex_flexInternal = apiToInternal[nodeIndex];
		int neighborCount = neighbourCounts[nodeIndex_flexInternal];

		// use hydrogen to spin
		if(synthase.timerH < 0.0f)
		{
			// check neighbours for hydrogen

			float nearestSqrdDistance = std::numeric_limits<float>::max();
			FGraphNodeHandle nearestHydrogen;

			int aboveCount = 0;
			int belowCount = 0;

			for (int ni = 0; ni < neighborCount; ++ni)
			{
				int neighborIndex = internalToAPI[neighbourIndices[ni*stride + nodeIndex_flexInternal]];

				FGraphNodeHandle handle(neighborIndex);

				if (FHydrogenGraphObject * hydrogen = hydrogens.componentPtrForNode(handle))
				{
					FGraphNode& node = graph->node(handle);

					FVector direction = node.position - atpSynthaseNode.position;

					bool isAbove = FVector::DotProduct(direction, atpSynthaseUp) > 0;

					aboveCount += isAbove ? 1 : 0;
					belowCount += isAbove ? 0 : 1;
				}
			}

			for(int ni = 0; ni < neighborCount; ++ni)
			{
				int neighborIndex = internalToAPI[neighbourIndices[ ni*stride + nodeIndex_flexInternal ]];

				FGraphNode& hydrogenNode = graph->node( neighborIndex );

				if( hydrogens.componentPtrForNode(hydrogenNode.handle()) == nullptr) continue;

				// check for a hydrogen gradient
				// this is naive, it just checks the local neighborhood
				// we want about 50% more hydrogen below the synthase than above it (the exit point)
				if (aboveCount * 1.5 > belowCount) continue;

				FVector direction = hydrogenNode.position - atpSynthaseNode.position;

				float dot = FVector::DotProduct(direction.GetSafeNormal(), atpSynthaseUp.GetSafeNormal());

				float angle = FMath::Acos(-dot);

				bool isAbove = dot > 0;

				const float maxAngle = FMath::DegreesToRadians(45.0f);


				// hydrogen enters through the bottom of ATPSynthase, through it and out the top
				// so, only interact with hydrogen below
				if(isAbove || angle > maxAngle)
					continue;

				float distanceSqrd = direction.SizeSquared();

				if(distanceSqrd > interactionRadiusSqrd)
					continue;

				if (distanceSqrd < nearestSqrdDistance)
				{
					nearestSqrdDistance = distanceSqrd;
					nearestHydrogen = FGraphNodeHandle(hydrogenNode);
				}

				break;
			}

			if (nearestHydrogen)
			{
				FGraphNode& hydrogenNode = nearestHydrogen(graph);

				// spin
				UEvent_SpinATPSynthase& spinEvent = timeline->recordEvent<UEvent_SpinATPSynthase>(hydrogenNode.position);
				spinEvent.triggeringAgent = FGraphNodeHandle(hydrogenNode);
				spinEvent.otherAgents.Add(FGraphNodeHandle(atpSynthaseNode));

				// setup spin state
				synthase.state = EATPSynthaseState::Spinning;

				synthase.timerSpin = synthase.timerSpinDuration;
				synthase.timerH = synthase.refractoryPeriodH;

				if (!atpSynthaseNode.hasComponent<FSpinnerGraphObject>())
				{
					atpSynthaseNode.addComponent<FSpinnerGraphObject>(*graph);
				}

				FSpinnerGraphObject& spinner = atpSynthaseNode.component<FSpinnerGraphObject>(*graph);

				spinner.angularVelocityMagnitude = 5.0f;

				// teleport hydrogen
				F3PointPathFollower& follower = hydrogenNode.addComponent<F3PointPathFollower>(*graph);

				follower.points[0] = hydrogenNode.position;
				follower.points[1] = atpSynthaseNode.position - atpSynthaseUp * .75f;
				follower.points[2] = atpSynthaseNode.position + atpSynthaseUp * 4.0f;

				follower.terminalVelocity = atpSynthaseUp * 4.0f;

				follower.duration = 1.0f;

				// disable particle interaction
				if (FFlexParticleObject * particle = particles.componentPtrForNode(nearestHydrogen))
				{
					follower.restoreChannel = particle->channel;
					particle->channel = 1;
				}

				if (hydrogenNode.hasComponent<FRandomWalkGraphObject>())
					hydrogenNode.component<FRandomWalkGraphObject>(*graph).timeLeft = 2.0f;
			}
		}

		// convert ADP to ATP
		if( synthase.state == EATPSynthaseState::Spinning && synthase.timerADP < 0.0f )
		{
			int neighborCount = neighbourCounts[nodeIndex_flexInternal];

			FGraphNodeHandle nearestADP;
			float nearestDistanceSqrd = std::numeric_limits<float>::max();

			for(int ni = 0; ni < neighborCount; ++ni)
			{
				int adpIndex = internalToAPI[neighbourIndices[ni*stride + nodeIndex_flexInternal]];

				FGraphNode& adpNode = graph->node( adpIndex );

				if(adpNode.hasComponent<FADPGraphObject>() && adpNode.hasComponent<FFlexParticleObject>() )
				{
					FVector direction = adpNode.position - atpSynthaseNode.position;

					bool isAbove = FVector::DotProduct( direction, atpSynthaseUp ) > 0;

					if(!isAbove)
						continue;

					float distanceSqrd = direction.SizeSquared();

					if(distanceSqrd > interactionRadiusSqrd)
						continue;

					if (distanceSqrd < nearestDistanceSqrd)
					{
						nearestDistanceSqrd = distanceSqrd;
						nearestADP = FGraphNodeHandle(adpNode);
					}
				}


			}

			if (nearestADP)
			{
				FGraphNode& adpNode = nearestADP(graph);

				// Spawn ATP
				// ---------
				const FVector atpPostion = atpSynthaseNode.position + atpSynthaseUp * 2.0;

				FGraphNode& atpNode = _spawnATP(atpPostion, adpNode.orientation, adpNode.scale);

				const FVector eventPosition = atpSynthaseNode.position + (adpNode.position - atpSynthaseNode.position).GetSafeNormal() * 3.0f;

				UEvent_SpawnATP& spawnATPEvent = timeline->recordEvent<UEvent_SpawnATP>(eventPosition);
				spawnATPEvent.triggeringAgent = FGraphNodeHandle(atpNode);
				spawnATPEvent.otherAgents.Add(FGraphNodeHandle(atpSynthaseNode));
				spawnATPEvent.otherAgents.Add(FGraphNodeHandle(adpNode));

				// launch the atp out
				atpNode.addComponent<FVelocityGraphObject>(*graph).linearVelocity = atpSynthaseUp * 3.0f;

				if (atpNode.hasComponent<FFlexParticleObject>())
				{
					FFlexParticleObject& atpParticle = atpNode.component<FFlexParticleObject>(*graph);
					FFlexParticleObject& adpParticle = adpNode.component<FFlexParticleObject>(*graph);

					atpParticle.channel = adpParticle.channel;
					atpParticle.group = adpParticle.group;
				}

				if (atpNode.hasComponent<FRandomWalkGraphObject>())
					atpNode.component<FRandomWalkGraphObject>(*graph).timeLeft = 2.0f;

				// kill the adp
				// TIM TODO We should remove ADP from the simulation, but for the event system, I have to keep it around as a zombie. We
				// should correctly handle this, perhaps tyeing into the simulation snapshots
				// graph->removeNode(adpNode.id);
				{
					adpNode.component<FGraphMesh>(graph).visible = false;

					// we're mutating, keep a copy
					auto componentsCached = adpNode.components;

					// remove all of the components, turn this guy into a zombie
					for (ComponentType type : componentsCached)
					{
						// we'll keep the graph mesh
						if (type == componentType<FGraphMesh>())
							continue;

						// but kill everything else
						adpNode.removeComponent(*graph, type);
					}

					// finally, hack the position to the event position
					adpNode.position = eventPosition;
				}

				synthase.timerADP = synthase.refractoryPeriodADP;
			}
		}

		if(synthase.timerSpin < 0.0f)
		{
			synthase.state = EATPSynthaseState::Inactive;

			if(atpSynthaseNode.hasComponent<FSpinnerGraphObject>())
			{
				atpSynthaseNode.component<FSpinnerGraphObject>(*graph).angularVelocityMagnitude = 0.0f;
			}

			if(atpSynthaseNode.hasComponent<FVelocityGraphObject>())
			{
				atpSynthaseNode.component<FVelocityGraphObject>( *graph ).angularVelocity = FVector::ZeroVector;
			}
		}

		// update timers
		synthase.timerH -= deltaT;
		synthase.timerADP -= deltaT;
		synthase.timerSpin -= deltaT;
	}

	_tickProtonPumps( deltaT, neighbourIndices, neighbourCounts, apiToInternal, internalToAPI, maxParticles );
	_tickFollowPath(deltaT);

	graph->endTransaction();
}



void UATPSynthaseSimulation::_tickFollowPath(float deltaT)
{
	auto& pathFollowers = graph->componentStorage<F3PointPathFollower>();
	auto& particles = graph->componentStorage<FFlexParticleObject>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	for (F3PointPathFollower& follower : pathFollowers)
	{
		if( !follower.isValid() ) continue;

		FGraphNode& node = graph->node(follower.nodeHandle());


		
		const int n = 3;

		const float t = follower.time / follower.duration;

		const int a_i = int(t * (n - 1));
		const int b_i = a_i + 1;

		// remainder
		const float u = t * (n - 1) - float(a_i);

		const FVector& a = follower.points[a_i];
		const FVector& b = follower.points[b_i];

		FVector last = node.position;
		node.position = a * (1.0f - u) + b * u;


		// the path is done, kill it
		if (follower.time > follower.duration)
		{
			if (FFlexParticleObject * particle = particles.componentPtrForNode(follower.nodeHandle()))
				particle->channel = follower.restoreChannel;

			if (FVelocityGraphObject * velocity = velocities.componentPtrForNode(follower.nodeHandle()))
				velocity->linearVelocity = follower.terminalVelocity;

			node.removeComponent<F3PointPathFollower>(*graph);
			continue;
		}




		follower.time += deltaT;
	}

}

void UATPSynthaseSimulation::_tickProtonPumps( 
	float deltaT,
	NvFlexVector<int>& neighbourIndices,
	NvFlexVector<int>& neighbourCounts,
	NvFlexVector<int>& apiToInternal,
	NvFlexVector<int>& internalToAPI,
	int maxParticles )
{
	auto& protonPumps = graph->componentStorage<FProtonPumpGraphObject>();
	auto& particles = graph->componentStorage<FFlexParticleObject>();

	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	const int stride = maxParticles;

	const float interactionRadiusSqrd = std::pow(hyrogenInteractionRadius, 2.0f);

	for(FProtonPumpGraphObject& pump : protonPumps)
	{
		if( !pump.isValid() ) continue;

		auto nodeIndex = pump.nodeIndex;
		FGraphNode& pumpNode = graph->node( nodeIndex );

		if(!pumpNode.hasComponent<FFlexParticleObject>())
			continue;

		FVector pumpUp = FVector::UpVector;
		pumpUp = pumpNode.orientation.RotateVector( pumpUp );


		int nodeIndex_flexInternal = apiToInternal[nodeIndex];
		int neighborCount = neighbourCounts[nodeIndex_flexInternal];

		if(pump.timerH < 0.0f)
		{
			// check neighbours for hydrogen

			FGraphNodeHandle nearestHydrogen;
			float nearestDistanceSqrd = std::numeric_limits<float>::max();

			for(int ni = 0; ni < neighborCount; ++ni)
			{
				int neighborIndex = internalToAPI[neighbourIndices[ni*stride + nodeIndex_flexInternal]];

				FGraphNode& hydrogenNode = graph->node( neighborIndex );

				if(!hydrogenNode.hasComponent<FHydrogenGraphObject>()) continue;
				if(!hydrogenNode.hasComponent<FVelocityGraphObject>()) continue;
				// if it's already pumped, don't pump it
				if(hydrogenNode.hasComponent<F3PointPathFollower>()) continue;

				FVector direction = hydrogenNode.position - pumpNode.position;

				float dot = FVector::DotProduct(direction.GetSafeNormal(), pumpUp.GetSafeNormal());

				float angle = FMath::Acos(dot);

				bool isBelow = dot < 0;

				const float maxAngle = FMath::DegreesToRadians(60.0f);

				// hydrogen through the top
				if(isBelow || angle > maxAngle)
					continue;

				float distanceSqrd = direction.SizeSquared();

				if(distanceSqrd > interactionRadiusSqrd)
					continue;

				if (distanceSqrd < nearestDistanceSqrd)
				{
					nearestDistanceSqrd = distanceSqrd;
					nearestHydrogen = FGraphNodeHandle(hydrogenNode.id);
				}
			}

			if (nearestHydrogen)
			{
				FGraphNode& hydrogenNode = graph->node(nearestHydrogen);

				// record the event before we mess with the position
				UEvent_ProtonPumped& protonPumpedEvent = timeline->recordEvent<UEvent_ProtonPumped>(hydrogenNode.position);

				protonPumpedEvent.triggeringAgent = FGraphNodeHandle(hydrogenNode);
				protonPumpedEvent.otherAgents.Add(FGraphNodeHandle(pumpNode));

				// teleport hydrogen
				F3PointPathFollower& follower = hydrogenNode.addComponent<F3PointPathFollower>(*graph);

				follower.points[0] = hydrogenNode.position;
				follower.points[1] = pumpNode.position + pumpUp * 2.5f;
				follower.points[2] = pumpNode.position - pumpUp * 2.5f;

				follower.terminalVelocity = -pumpUp * 3.0f;

				follower.duration = 1.0f;

				// disable particle interaction
				if (FFlexParticleObject * particle = particles.componentPtrForNode(nearestHydrogen))
				{
					follower.restoreChannel = particle->channel;
					particle->channel = 1;
				}

				if (hydrogenNode.hasComponent<FRandomWalkGraphObject>())
					hydrogenNode.component<FRandomWalkGraphObject>(*graph).timeLeft = 2.0f;


				pump.timerH = pump.refractoryPeriodH;
			}
		}

		// update timers
		pump.timerH -= deltaT;
	}
}

FGraphNode& UATPSynthaseSimulation::_spawnATP( FVector position, FQuat orientation, float scale )
{
	FGraphNode& node = graph->node(graph->addNode( position, orientation, scale ));

	// spawn our components
	for(FTimStructBox& box : atpTemplate)
	{
		if(!box.IsValid())
			continue;

		// create the graph object in the simulation
		ComponentType type = FGraphObject::componentType( box.scriptStruct );

		FGraphObject * object = node.addComponent( *graph, type );

		// then copy the memory from the element
		box.scriptStruct->CopyScriptStruct( object, box.structMemory );

		object->nodeIndex = node.id;
	}

	return node;
}

void FFlexRigidBodyObject::_ensureFlexParticle(FGraphNode& node, FGraph& graph)
{
	if (!node.hasComponent<FFlexParticleObject>())
		node.addComponent<FFlexParticleObject>(graph);

	if (!node.hasComponent<FVelocityGraphObject>())
		node.addComponent<FVelocityGraphObject>(graph);
}


FVector FFlexRigidBodyObject::calculateCenterOfMass(FGraph& graph)
{
	FGraphNode& node = graph.node(nodeHandle());

	auto rigidEdgeStorage = graph.edgeStorage<FFlexRigidBodyConnection>();
	auto& particles = graph.componentStorage<FFlexParticleObject>();

	FVector sum = FVector::ZeroVector;

	int32 count = 0;
	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		FFlexRigidBodyConnection * rigidConnection = rigidEdgeStorage.objectPtr(edgeHandle);

		if (!rigidConnection) continue;

		auto other = graph.edge(edgeHandle).other(node.handle());

		if (!particles.componentPtrForNode(other)) continue;

		sum += other.node(graph).position;

		count++;
	}

	if (count) sum /= float(count);

	return sum;
}

void FFlexRigidBodyObject::applyRotationTranslation(const FQuat& rotation, const FVector& translation, FGraph& graph)
{
	FGraphNode& node = graph.node(nodeHandle());

	const FVector pr = node.position; // original node position
	const FVector pr_ = pr + translation; // new node position

	node.position = pr_;
	node.orientation = (rotation * node.orientation).GetNormalized();

	nodeHandle().node(graph).each<FFlexRigidBodyConnection>(graph, [&](FGraphNodeHandle subRigidHandle, FFlexRigidBodyConnection& rigidConnection) {
		FGraphNode& subNode = subRigidHandle(graph);

		subNode.position = rotation.RotateVector(subNode.position - pr) + pr_;

		FQuat qsub = rotation * subNode.orientation;
		subNode.orientation = qsub.GetNormalized();
	});
}

// So, passing FQuat by value can cause a compiler bug: https://answers.unrealengine.com/questions/125784/doing-a-fquatfquat-when-one-of-them-is-a-parameter.html
// The quat passed in is not aligned, and things blow up in the SSE instructions.
void FFlexRigidBodyObject::applyRotationTranslationInferVelocity(const FQuat& rotation, const FVector& translation, FGraph& graph)
{
	auto& velocities = graph.componentStorage<FVelocityGraphObject>();

	FGraphNode& node = graph.node(nodeHandle());

	const FVector pr = node.position; // original node position
	const FVector pr_ = pr + translation; // new node position


	node.position = pr_;
	node.orientation = (rotation * node.orientation).GetNormalized();

	if (FVelocityGraphObject * nodeVelocity = velocities.componentPtrForNode(nodeHandle())) nodeVelocity->linearVelocity += translation;

	nodeHandle().node(graph).each<FFlexRigidBodyConnection>(graph, [&](FGraphNodeHandle subRigidHandle, FFlexRigidBodyConnection& rigidConnection) {
		FGraphNode& subNode = subRigidHandle(graph);

		subNode.position = rotation.RotateVector(subNode.position - pr) + pr_;

		FQuat qsub = rotation * subNode.orientation;
		subNode.orientation = qsub.GetNormalized();

		FVelocityGraphObject * velocity = velocities.componentPtrForNode(subRigidHandle);

		if( velocity ) velocity->linearVelocity += translation;
	});
}

void FFlexRigidBodyObject::setRotationPosition(const FQuat rotation_unrealBug, const FVector position, FGraph& graph)
{
	FGraphNode & node = graph.node(nodeHandle());

	FVector deltaT = position - node.position;
	const FQuat rotation_fixedUnrealBug = rotation_unrealBug;
	FQuat deltaQ = rotation_fixedUnrealBug * node.orientation.Inverse();

	applyRotationTranslation(deltaQ, deltaT, graph);
}

void FFlexRigidBodyObject::setRotationPositionInferVelocity(const FQuat rotation_unrealBug, const FVector position, FGraph& graph)
{
	FGraphNode & node = graph.node(nodeHandle());

	FVector deltaT = position - node.position;
	const FQuat rotation_fixedUnrealBug = rotation_unrealBug;
	FQuat deltaQ = rotation_fixedUnrealBug * node.orientation.Inverse();

	applyRotationTranslationInferVelocity(deltaQ, deltaT, graph);
}

FFlexRigidBodyObject& FFlexRigidBodyObject::createRigidBody(FGraph& graph, TArray<FGraphNodeHandle>& particlesInBody, FGraphNodeHandle rigidNodeHandle)
{
	checkfSlow(particlesInBody.Num() > 2, TEXT("There must be more than two particles in the body."));

	graph.beginTransaction();

	if (!rigidNodeHandle)
		rigidNodeHandle = FGraphNodeHandle(graph.addNode(FVector::ZeroVector));

	checkfSlow(!rigidNodeHandle(graph).hasComponent<FFlexParticleObject>(), TEXT("The rigidNodeHandle must not have a particle."));

	FVector sum = FVector::ZeroVector;

	auto& particles = graph.componentStorage<FFlexRigidBodyObject>();

	for (auto& other : particlesInBody)
	{
		graph.connectNodes<FFlexRigidBodyConnection>(rigidNodeHandle, other);

		_ensureFlexParticle(other(graph), graph);

		other(graph).component<FFlexParticleObject>(graph).group = rigidNodeHandle.index;

		sum += other(graph).position;
	}

	// configure the rigid body node
	FVector com = sum / float(particlesInBody.Num());

	FGraphNode& rigidNode = rigidNodeHandle(graph);

	FFlexRigidBodyObject& body = rigidNode.addComponent<FFlexRigidBodyObject>(graph);
	body.stiffness = 0.9f;
	body.plasticDeformationCreep = 0.01f;
	body.plasticDeformationThreshold = 5.0f;

	FVelocityGraphObject& velocity = rigidNode.addComponent<FVelocityGraphObject>(graph);

	rigidNode.position = com;

	graph.endTransaction();

	return body;
}

void FFlexRigidBodyObject::edgesAndNodes(FGraph& graph, TArray<FGraphNodeHandle>& nodes_out, TArray<FGraphEdgeHandle>& edges_out)
{
	auto rigidConnections = graph.edgeStorage<FFlexRigidBodyConnection>();

	FGraphNode& node = graph.node(nodeHandle());

	nodes_out.Add(nodeHandle());

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if (!rigidConnections.objectPtr(edgeHandle)) continue;

		edges_out.Add(edgeHandle);

		FGraphEdge& edge = graph.edge(edgeHandle);

		nodes_out.Add(edge.other(nodeHandle()));
	}
}


FGraphNodeHandle FFlexRigidBodyObject::getRigidBodyHandle(FGraph& graph, FGraphNodeHandle subNodeHandle)
{
	auto& node = graph.node(subNodeHandle);

	if (node.hasComponent<FFlexRigidBodyObject>())
		return subNodeHandle;

	auto edgeStorage = graph.edgeStorage<FFlexRigidBodyConnection>();

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if (edgeStorage.isValid(edgeHandle) && edgeStorage.objectPtr(edgeHandle))
		{
			FGraphEdge& edge = graph.edge(edgeHandle);

			if( edge.isValid() )
				return edge.other(subNodeHandle);
		}
	}

	return FGraphNodeHandle::null;
}

void UStaticPositionSimulation::attach()
{

}

void UStaticPositionSimulation::tick(float deltaT)
{
	_tickStatic(deltaT);
	_tickStabalized(deltaT);
}

void UStaticPositionSimulation::tick_paused(float deltaT)
{
	tick(deltaT);
}

void UStaticPositionSimulation::_tickStabalized(float deltaT)
{
	auto& staticPositions = graph->componentStorage<FStabalizedPosition>();
	auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	// initialize position
	for (FStabalizedPosition& staticObject : staticPositions)
	{
		if (!staticObject.isValid()) continue;
		if (staticObject.didLoad) continue;

		FGraphNode& node = graph->node(staticObject.nodeHandle());

		staticObject.position = node.position;
		staticObject.orientation = node.orientation;

		staticObject.didLoad = true;
	}

	// set the static positions
	for (FStabalizedPosition& staticObject : staticPositions)
	{
		if (!staticObject.isValid()) continue;

		FGraphNode& node = graph->node(staticObject.nodeHandle());

		if (!node.isValid()) continue;

		FVector position = FMath::Lerp(node.position, staticObject.position, deltaT * staticObject.strength);
		FQuat orientation = FMath::Lerp(node.orientation, staticObject.orientation, deltaT * staticObject.strength);


		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, staticObject.nodeHandle());

		if (rigidHandle)
		{
			FQuat dq = orientation * node.orientation.Inverse();
			FVector dp = position - node.position;

			FFlexRigidBodyObject * rigid = rigids.componentPtrForNode(rigidHandle);

			rigid->applyRotationTranslationInferVelocity(dq, dp, *graph);
		}
		else
		{
			FVector dv = position - node.position;
			FQuat dq = orientation * node.orientation.Inverse();
		

			node.position = position;
			node.orientation = orientation;

			if (FVelocityGraphObject * velocity = velocities.componentPtrForNode(staticObject.nodeHandle()))
			{
				if (deltaT == 0.0f) deltaT = 0.001f;
				velocity->linearVelocity += dv / deltaT;
			}
		}
	}

	
}

void UStaticPositionSimulation::_tickStatic(float deltaT)
{
	auto& staticPositions = graph->componentStorage<FStaticPositionObject>();
	auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	// initialize position
	for (FStaticPositionObject& staticObject : staticPositions)
	{
		if (!staticObject.isValid()) continue;
		if (staticObject.didLoad) continue;

		FGraphNode& node = graph->node(staticObject.nodeHandle());

		staticObject.position = node.position;
		staticObject.orientation = node.orientation;

		staticObject.didLoad = true;
	}

	// set the static positions
	for (FStaticPositionObject& staticObject : staticPositions)
	{
		if (!staticObject.isValid()) continue;

		FGraphNode& node = graph->node(staticObject.nodeHandle());

		if (!node.isValid()) continue;


		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, staticObject.nodeHandle());

		if (rigidHandle)
		{
			FQuat dq = staticObject.orientation * node.orientation.Inverse();
			FVector dp = staticObject.position - node.position;

			FFlexRigidBodyObject * rigid = rigids.componentPtrForNode(rigidHandle);

			rigid->applyRotationTranslation(dq, dp, *graph);
		}
		else
		{
			FVector dv = staticObject.position - node.position;
			
			node.position = staticObject.position;
			node.orientation = staticObject.orientation;

			node.orientation = node.orientation.GetNormalized();

			// we add velocity for flex connected objects, otherwise they keep pulling this node
			if (FVelocityGraphObject * velocity = velocities.componentPtrForNode(staticObject.nodeHandle()))
			{
				if (deltaT == 0.0f) deltaT = 0.001f;

				velocity->linearVelocity += dv / deltaT;
			}
		}
	}
}
