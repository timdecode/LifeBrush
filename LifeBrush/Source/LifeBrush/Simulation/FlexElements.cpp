// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "Visualization/Timeline.h"
#include "SimulationSnapshotActor.h"

#include "Utility.h"
#include "SimulationSnapshotActor.h"




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

	if( !_flexSimulation)
		initFlexSimulationObject();

	_flexSimulation->tick(DeltaTime);
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

	initFlexSimulationObject();
}

void UFlexSimulationComponent::BeginPlay()
{
	Super::BeginPlay();


}

void UFlexSimulationComponent::initFlexSimulationObject()
{
	if (_flexSimulation)
		return;

	initGraph();

	_flexSimulation = std::make_unique<FFlexSimulation>(
		graph,
		*simulationManager,
		flexParams
	);

	_readRules();

	_flexSimulation->initSimulationManager(GetOwner());
}

void UFlexSimulationComponent::_readRules()
{
	//if (!rulesActor) return;

	//UMLElementSimulation * ruleSimulation = _flexSimulation->simulationManager.registerSimulation<UMLElementSimulation>();

	//ruleSimulation->ruleGraph.init();
	//ruleSimulation->ruleGraph.clear();

	//TArray<USceneComponent*> childSceneComponents;
	//rulesActor->GetRootComponent()->GetChildrenComponents(false, childSceneComponents);

	//for (USceneComponent * child : childSceneComponents)
	//{
	//	AMLRuleActor * ruleActor = Cast<AMLRuleActor>(child->GetOwner());

	//	if( !ruleActor ) continue;

	//	ruleActor->writeToGraph(ruleSimulation->ruleGraph);
	//}
}

auto UFlexSimulationComponent::init(
	std::shared_ptr<tcodsMeshInterfaceBase> meshInterface,
	FTransform meshInterfaceToWorld,
	UCameraComponent * camera) -> void
{

	_flexSimulation->initMeshInterface(meshInterface, meshInterfaceToWorld, camera);
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

	// create scene root
	USceneComponent * newScene = NewObject<USceneComponent>(newActor);
	{
		AActor * originalActor = GetOwner();

		USceneComponent * originalScene = originalActor->GetRootComponent();

		FTransform newTransform = originalScene->GetComponentTransform();

		newActor->SetRootComponent(newScene);
		newActor->AddInstanceComponent(newScene);

		newScene->SetWorldTransform(newTransform);
	}

	// duplicate the simulation
	UFlexSimulationComponent * newElements = DuplicateObject(this, newActor);

	newActor->AddInstanceComponent(newElements);

	return newActor;
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
		_loadPositionsAndVelocities(particles, velocities, phases, active);

		// make sure we have created our buffers
		_springs->init();
	}

	_loadSprings();

	// update the push-interaction sphere
	(FVector4&)shapePositions[0] = _spherePosition;
	shapePositions[0].W = 1.0f;
	geometry[0].sphere.radius = _sphereRadius;

	// Read Flex state into the simulation
	_readFlexState(particles, velocities, phases);
	_spawnDespawnSprings();

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
		_rigids->mapRotationsTranslations();

		_readRigidRotations();
	}

	// simulate
	_preTick(deltaT);

	//UTimelineSimulation * timeline = simulationManager.simulation<UTimelineSimulation>();

	//timeline->beginFrame();
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

	// called immediately after the simulation, to generate new particles
	_spawnDespawnParticles(particles, velocities, phases, active);

	// Write simulation state back into Flex
	_writeFlexState(particles, velocities, phases, active);



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
		NvFlexSetSprings(_solver, _springs->indices.buffer, _springs->lengths.buffer, _springs->coefficients.buffer, _springs->lengths.size());

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

auto FFlexSimulation::initMeshInterface(
	std::shared_ptr<tcodsMeshInterfaceBase> meshInterface,
	FTransform meshInterfaceToWorld,
	UCameraComponent * camera) -> void
{
	checkfSlow(_didInitSImulationManager, TEXT("initSimulationManager must be called before initMeshInterface."));

	_meshInterface = meshInterface;

	_meshToWorld = meshInterfaceToWorld;

	_needsFlexReset = true;
	_didInitFlex = false;
	_didSpawnShapesAndMesh = false;

	URandomWalkSimulation * randomWalkSim = simulationManager.simulation<URandomWalkSimulation>();

	randomWalkSim->meshInterface = _meshInterface;

	UVisualization_AgentPathLines * pathLines = simulationManager.simulation<UVisualization_AgentPathLines>();

	pathLines->camera = camera;
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

auto FFlexSimulation::loadExportedDomainInfo( OutputDomainExport& exportInfo ) -> void
{
	graphSimulation.clear();

	graphSimulation.beginTransaction();

	_meshInterface = exportInfo.meshInterface;

	FTransform worldToUsTransform = owner->GetRootComponent()->GetComponentTransform().Inverse();

	FTransform toWorld = owner->GetRootComponent()->GetComponentTransform();

	exportInfo.write(graphSimulation, worldToUsTransform);

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
		FElementObject * elementObject = node.hasComponent<FElementObject>() ? &node.component<FElementObject>(graphSimulation) : nullptr;
		if (elementObject && elementObject->surfaceIndex.isOnSurface())
		{
			FSurfaceBoundGraphObject& boundObject = node.addComponent<FSurfaceBoundGraphObject>(graphSimulation);

			FVector p = toWorld.TransformPosition(node.position);

			auto rotationAndNormal = _meshInterface->rotationAndNormalAtIndex(elementObject->surfaceIndex);

			boundObject.lastSurfaceRotation = rotationAndNormal.first;
			boundObject.surfaceIndex = elementObject->surfaceIndex;

			onSurface = true;
		}

		// set particle
		if (node.hasComponent<FFlexParticleObject>())
		{
			FFlexParticleObject& particleObject = node.component<FFlexParticleObject>(graphSimulation);

			// don't interact with the surface
			particleObject.channel = onSurface ? eNvFlexPhaseShapeChannel1 : eNvFlexPhaseShapeChannel0;
		}

		if (node.hasComponent<FStaticPositionObject>())
		{
			node.component<FStaticPositionObject>(graphSimulation).position = node.position;
		}
	}

	// setup the positions

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
	
	// setup listeners
	graphSimulation.addComponentListener<FFlexParticleObject>((ComponentListener*) this );
	graphSimulation.addComponentListener<FFlexRigidBodyObject>((ComponentListener*)this);

	graphSimulation.addEdgeObjectListener<FFlexConnection>(this);
	graphSimulation.addEdgeObjectListener<FFlexRigidBodyConnection>(this);

	simulationManager.attachSimulations();

	_didInitSImulationManager = true;
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

auto FFlexSimulation::_spawnDespawnParticles( FVector4 * positions, FVector * velocities, int * phases, int * active ) -> void
{
	float r = _params.radius;

	nParticles = graphSimulation.componentStorage<FFlexParticleObject>().size();

	if(_flexParticlesToAdd.empty() && _flexParticlesToRemove.empty())
		return;

	// spawn particles
	for(FGraphNodeHandle handle : _flexParticlesToAdd)
	{
		FGraphNode& node = handle(graphSimulation);

		if(!node.isValid())
			continue;

		if(	!node.hasComponent<FVelocityGraphObject>() ||
			!node.hasComponent<FFlexParticleObject>() )
			continue;

		FFlexParticleObject& particleObject = node.component<FFlexParticleObject>(graphSimulation);
		FVelocityGraphObject& velocityObject = node.component<FVelocityGraphObject>(graphSimulation);

		FVector p = node.position;
		FVector v = velocityObject.linearVelocity;

		positions[node.id] = p;
		positions[node.id].W = 0.125;
		
		velocities[node.id] = v;
		
		phases[node.id] = NvFlexMakePhaseWithChannels( particleObject.group, particleObject.selfCollide ? eNvFlexPhaseSelfCollide : 0, particleObject.channel );
	}

	// We don't need to remove particles, just deactivate them.
	// We do this by update the active set for the entire simulation---brute force, 
	// but we do more expensive things like reading every particle position.
	{
		auto& flexParticles = graphSimulation.componentStorage<FFlexParticleObject>();

		size_t i = 0;
		for( FFlexParticleObject& particle : flexParticles )
		{
			active[i] = particle.nodeIndex;

			++i;
		}

		// We don't have to worry about setting inactive particles, because NvFlexSetActiveCount will truncate the active array
		// to just the indices we set above.
	}

	_flexParticlesToAdd.clear();
	_flexParticlesToRemove.clear(); // we don't need to parse this guy, it's handled implicitly by the above loop
}

auto FFlexSimulation::_spawnDespawnSprings() -> void
{
	if (_flexSpringsToAdd.empty() && _flexSpringsToRemove.empty())
		return;

	// hack for now, just load everything again
	_loadSprings();

	_flexSpringsToAdd.clear();
	_flexSpringsToRemove.clear();
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

void FFlexSimulation::_loadPositionsAndVelocities(FVector4 * positions, FVector * velocities, int * phases, int * active)
{
	auto& flexParticles = graphSimulation.componentStorage<FFlexParticleObject>();

	nParticles = flexParticles.size();

	for (FFlexParticleObject& particleObject : flexParticles)
	{
		phases[particleObject.nodeIndex] = NvFlexMakePhaseWithChannels(particleObject.group, particleObject.selfCollide ? eNvFlexPhaseSelfCollide : 0, particleObject.channel);

		positions[particleObject.nodeIndex] = graphSimulation.node(particleObject.nodeIndex).position;
		positions[particleObject.nodeIndex].W = particleObject.inverseMass;
	}

	// activate
	size_t i = 0;
	for (FFlexParticleObject& particle : flexParticles)
	{
		active[i] = particle.nodeIndex;

		++i;
	}

	auto& flexVelocities = graphSimulation.componentStorage<FVelocityGraphObject>();

	for (FVelocityGraphObject& velocityObject : flexVelocities)
		velocities[velocityObject.nodeIndex] = velocityObject.linearVelocity;
}




void FFlexSimulation::_loadSprings()
{
	auto n = graphSimulation.edgeStorage<FFlexConnection>().validSize();

	auto& indices = _springs->indices;
	auto& lengths = _springs->lengths;
	auto& coefficients = _springs->coefficients;

	indices.resize(n * 2);
	lengths.resize(n);
	coefficients.resize(n);

	auto connections = graphSimulation.edgeView<FFlexConnection>();

	int i = 0;
	connections.each([&](FFlexConnection& connection, FGraphEdge& edge)
	{
		indices[(i * 2) + 0] = edge.a;
		indices[(i * 2) + 1] = edge.b;

		FGraphNode& a = graphSimulation.node(edge.a);
		FGraphNode& b = graphSimulation.node(edge.b);

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

	shapeFlags[meshIndex] = NvFlexMakeShapeFlagsWithChannels( eNvFlexShapeTriangleMesh, false, eNvFlexPhaseShapeChannel0 | eNvFlexPhaseShapeChannel1 );

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
		auto id = velocityObject.nodeIndex;

		FVector v = velocities[id];

		velocityObject.linearVelocity = v;
	}
}




auto FFlexSimulation::_writeFlexState( FVector4 * positions, FVector * velocities, int * phases, int * active ) -> void
{
	// write positions
	auto& particles = graphSimulation.componentStorage<FFlexParticleObject>();

	size_t i = 0;

	for(FFlexParticleObject& particle : particles)
	{
		FGraphNode& node = graphSimulation.node( particle.nodeIndex );

		FVector4& p = positions[node.id];

		float w = p.W;

		p = FVector4( node.position, w );

		phases[i] = NvFlexMakePhaseWithChannels(particle.group, particle.selfCollide ? eNvFlexPhaseSelfCollide : 0, particle.channel);

		active[i] = particle.nodeIndex;

		++i;
	}

	// write velocities
	auto& velocityObjects = graphSimulation.componentStorage<FVelocityGraphObject>();

	for(FVelocityGraphObject& velocityObject : velocityObjects)
	{
		auto id = velocityObject.nodeIndex;

		velocities[id] = velocityObject.linearVelocity;
	}
}

auto FFlexSimulation::_integrateRotations(float deltaT) -> void
{
	auto& velocityObjects = graphSimulation.componentStorage<FVelocityGraphObject>();

	for(FVelocityGraphObject& velocityObject : velocityObjects)
	{
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

auto FFlexSimulation::_loadTcodsMesh(tcodsMeshInterfaceBase& tcodsMesh) -> void
{
	FTransform worldToUsTransform = owner->GetRootComponent()->GetComponentTransform().Inverse();

	FTransform toUsTransform = _meshToWorld * worldToUsTransform;

	// data to be consumed by flex
	std::vector<FVector4> flexPositions;
	std::vector<int> flexIndices;

	_mesh0_bounds = FBox(EForceInit::ForceInitToZero);

	size_t sectionOffset = 0;

	for (auto& ptr : tcodsMesh._sections)
	{
		using namespace tcods;

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
			flexIndices.push_back(triangle[2]);
			flexIndices.push_back(triangle[1]);
			flexIndices.push_back(triangle[0]);
		}

		sectionOffset += flexPositions.size();
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

	if (type == _FFlexParticleObjectType)
	{
		_flexParticlesToAdd.insert(handle);
		_flexParticlesToRemove.erase(handle);

		FGraphNode& node = graphSimulation.node(handle);

		FFlexParticleObject& particle = node.component<FFlexParticleObject>(graphSimulation);

		bool onSurface = node.hasComponent<FSurfaceBoundGraphObject>();

		// don't interact with the surface
		particle.channel = onSurface ? eNvFlexPhaseShapeChannel1 : eNvFlexPhaseShapeChannel0;
	}
	else if (type == _FFlexRigidObjectType)
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
	if (type != _FFlexParticleObjectType)
	{
		_flexParticlesToRemove.insert(node);
		_flexParticlesToAdd.erase(node);
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
		_flexSpringsToAdd.insert(handle);
		_flexSpringsToRemove.erase(handle);
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
		_flexSpringsToRemove.insert(handle);
		_flexSpringsToAdd.erase(handle);
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


	for(FSurfaceBoundGraphObject& dweller : surfaceDwellers)
	{
		FGraphNode& node = graph->node( dweller.nodeIndex );

		FVector p = toWorld.TransformPosition( node.position );

		auto nearest = meshInterface->nearestPointOnMesh( p );

		auto rotationAndNormal = meshInterface->rotationAndNormalAtIndex( nearest.surfaceIndex );

		FQuat surfaceRotation = toWorld.InverseTransformRotation(rotationAndNormal.first);

		FQuat rotation = surfaceRotation * dweller.lastSurfaceRotation.Inverse();

		node.position = toWorld.InverseTransformPosition( nearest.point );
		node.orientation = rotation * node.orientation;

		dweller.lastSurfaceRotation = surfaceRotation;

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

	auto& staticPositions = graph->componentStorage<FStaticPositionObject>();

	for(FStaticPositionObject& staticObject : staticPositions)
	{
		FGraphNode& node = graph->node( staticObject.nodeIndex );

		node.position = staticObject.position;

		if(node.hasComponent<FVelocityGraphObject>())
			node.component<FVelocityGraphObject>(*graph ).linearVelocity = FVector::ZeroVector;
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

	const int stride = maxParticles;

	// If we have hydrogen below us, spin and teleport it across the membrane
	// If we are spinning and have ADP above us, destroy it and produce ADP

	// we'll be adding and removing components, so start a transaction
	graph->beginTransaction();

	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	for(FATPSynthaseGraphObject& synthase : synthases)
	{
		auto nodeIndex = synthase.nodeIndex;
		FGraphNode& atpSynthaseNode = graph->node(nodeIndex);

		if(!atpSynthaseNode.hasComponent<FFlexParticleObject>())
			continue;

		FVector atpSynthaseUp = FVector::UpVector;
		atpSynthaseUp = atpSynthaseNode.orientation.RotateVector( atpSynthaseUp );

		const float interactionRadiusSqrd = synthase.hyrogenInteractionRadius * synthase.hyrogenInteractionRadius;

		int nodeIndex_flexInternal = apiToInternal[nodeIndex];

		// use hydrogen to spin
		if(synthase.timerH < 0.0f)
		{
			// check neighbours for hydrogen
			int neighborCount = neighbourCounts[nodeIndex_flexInternal];

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

				bool isAbove = FVector::DotProduct( direction, atpSynthaseUp ) > 0;

				// hydrogen enters through the bottom of ATPSynthase, through it and out the top
				// so, only interact with hydrogen below
				if(isAbove)
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

				// teleport hydrogen and spin
				// ---
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
				hydrogenNode.position = atpSynthaseNode.position + atpSynthaseUp * 4.0;

				if (hydrogenNode.hasComponent<FVelocityGraphObject>())
					hydrogenNode.component<FVelocityGraphObject>(*graph).linearVelocity = atpSynthaseUp * 4.0f;

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

	graph->endTransaction();
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

	UTimelineSimulation * timeline = simulationManager->simulation<UTimelineSimulation>();

	const int stride = maxParticles;

	for(FProtonPumpGraphObject& pump : protonPumps)
	{
		auto nodeIndex = pump.nodeIndex;
		FGraphNode& pumpNode = graph->node( nodeIndex );

		if(!pumpNode.hasComponent<FFlexParticleObject>())
			continue;

		FVector pumpUp = FVector::UpVector;
		pumpUp = pumpNode.orientation.RotateVector( pumpUp );

		const float interactionRadiusSqrd = pump.hyrogenInteractionRadius * pump.hyrogenInteractionRadius;

		int nodeIndex_flexInternal = apiToInternal[nodeIndex];

		// use hydrogen to spin
		if(pump.timerH < 0.0f)
		{
			// check neighbours for hydrogen
			int neighborCount = neighbourCounts[nodeIndex_flexInternal];

			FGraphNodeHandle nearest;
			float nearestDistanceSqrd = std::numeric_limits<float>::max();

			for(int ni = 0; ni < neighborCount; ++ni)
			{
				int neighborIndex = internalToAPI[neighbourIndices[ni*stride + nodeIndex_flexInternal]];

				FGraphNode& hydrogenNode = graph->node( neighborIndex );

				if(!hydrogenNode.hasComponent<FHydrogenGraphObject>())
					continue;

				if (!hydrogenNode.hasComponent<FVelocityGraphObject>())
					continue;

				FVector direction = hydrogenNode.position - pumpNode.position;

				bool isAbove = FVector::DotProduct( direction, pumpUp ) > 0;

				// hydrogen through the top
				if(!isAbove)
					continue;

				// apply acceleration
				//FVelocityGraphObject& velocityObject = hydrogenNode.component<FVelocityGraphObject>(*graph);

				float distanceSqrd = direction.SizeSquared();

				//float a = _g / distanceSqrd;

				//velocityObject.linearVelocity = velocityObject.linearVelocity + (direction.GetSafeNormal() * -a * deltaT);


				if(distanceSqrd > interactionRadiusSqrd)
					continue;

				if (distanceSqrd < nearestDistanceSqrd)
				{
					nearestDistanceSqrd = distanceSqrd;
					nearest = FGraphNodeHandle(hydrogenNode.id);
				}
			}

			if (nearest)
			{
				FGraphNode& hydrogenNode = nearest(graph);

				// record the event before we mess with the position
				UEvent_ProtonPumped& protonPumpedEvent = timeline->recordEvent<UEvent_ProtonPumped>(hydrogenNode.position);

				protonPumpedEvent.triggeringAgent = FGraphNodeHandle(hydrogenNode);
				protonPumpedEvent.otherAgents.Add(FGraphNodeHandle(pumpNode));

				// teleport hydrogen
				hydrogenNode.position = pumpNode.position - pumpUp * 0.1;

				FVelocityGraphObject& velocityObject = hydrogenNode.component<FVelocityGraphObject>(*graph);

				velocityObject.linearVelocity = -pumpUp * 4.0f;

				if (hydrogenNode.hasComponent<FRandomWalkGraphObject>())
					hydrogenNode.component<FRandomWalkGraphObject>(*graph).timeLeft = 1.0f;

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

	// cleaner to read, these will be inlined by the compiler
	auto translate = [rotation, pr, pr_](FVector p) -> FVector {
		// find the original direction vector, rotate and add it to the new position
		return rotation.RotateVector(p - pr) + pr_;
	};

	auto rotate = [rotation](FQuat q) -> FQuat {
		return rotation * q;
	};

	node.position = translate(node.position);
	node.orientation = rotate(node.orientation);

	nodeHandle().node(graph).each<FFlexRigidBodyConnection>(graph, [&](FGraphNodeHandle subRigidHandle, FFlexRigidBodyConnection& rigidConnection) {
		FGraphNode& subNode = subRigidHandle(graph);

		subNode.position = translate(subNode.position);
		subNode.orientation = rotate(subNode.orientation);
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
	node.orientation = rotation * node.orientation;

	if (FVelocityGraphObject * nodeVelocity = velocities.componentPtrForNode(nodeHandle())) nodeVelocity->linearVelocity += translation;

	nodeHandle().node(graph).each<FFlexRigidBodyConnection>(graph, [&](FGraphNodeHandle subRigidHandle, FFlexRigidBodyConnection& rigidConnection) {
		FGraphNode& subNode = subRigidHandle(graph);

		subNode.position = rotation.RotateVector(subNode.position - pr) + pr_;

		FQuat qsub = rotation * subNode.orientation;
		subNode.orientation = qsub;

		FVelocityGraphObject * velocity = velocities.componentPtrForNode(subRigidHandle);

		if( velocity ) velocity->linearVelocity += translation;
	});
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

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if (graph.edgeObjectPtr<FFlexRigidBodyConnection>(edgeHandle))
		{
			FGraphEdge& edge = graph.edge(edgeHandle);

			return edge.other(subNodeHandle);
		}
	}

	return FGraphNodeHandle::null;
}

