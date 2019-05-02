//
//  DiscreteElementEditorComponent.cpp
//  LifeBrush
//
//  Created by Timothy Davison on 2018-12-21.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementGenerator.h"

#include "DiscreteElementEditorComponent.h"
#include "ShipEditorSimulation/GraphSnapshot.h"
#include "tcodsMeshInterface.h"

#include "RegionGrowingComponent.h"
#include "VolumeComponent.h"
#include "RegionGrowingGenerator.h"

UDiscreteElementEditorComponent::UDiscreteElementEditorComponent()
{
	bWantsInitializeComponent = true;
	PrimaryComponentTick.bCanEverTick = true;
	PrimaryComponentTick.TickGroup = ETickingGroup::TG_PrePhysics;
	bIsActive = true;
	bTickInEditor = true;

	PostPhysicsTick.TickGroup = ETickingGroup::TG_PostPhysics;
}

void UDiscreteElementEditorComponent::InitializeComponent()
{
	Super::InitializeComponent();

	generators.Remove(nullptr);

	_initContextAndGraph();

	if (!_didInit)
		_init();
}

void UDiscreteElementEditorComponent::BeginPlay()
{
	Super::BeginPlay();

	context->limits = limits;

	// We can't init flex until after initialize component, other guys are going to set our camera.
	// Hence this code doesn't go into _init().
	_flexSimulation = std::make_unique<FFlexSimulation>(graph, *simulationManager, flexParams);
	_flexSimulation->initSimulationManager(GetOwner());
	_flexSimulation->initMeshInterface(context->meshInterface, FTransform::Identity, camera);
	_flexSimulation->setInstanceManagerBounds(limits);

	_flexSimulation->begin();
	_flexSimulation->play();
}

void UDiscreteElementEditorComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	// in case we tick in the editor, we fire this off. it uses conditional initialization
	_initContextAndGraph();

	if (ThisTickFunction == &PrimaryComponentTick)
		preTickPhysics(DeltaTime);
	else if (ThisTickFunction == &PostPhysicsTick)
		postTickPhysics(DeltaTime);
}

void UDiscreteElementEditorComponent::_initContextAndGraph()
{
	if (context)
		return;

	graph.init();

	context = std::make_unique<SynthesisContext>(graph);
}

void UDiscreteElementEditorComponent::preTickPhysics(float DeltaTime)
{
	_updateLimits();

#if WITH_EDITOR
	// don't tick, in the editor, past this point
	bool inEditor = GetWorld()->WorldType == EWorldType::Editor;

	if (inEditor)
		return;
#endif

	if (_generator && !_paused)
	{
		// The generator wants a flex simulation wrapped around it.
		// So, add its work to the tick-work queue and tick flex.
		if (_generator->wantsFlex())
		{
			_flexSimulation->addTickWork([&, DeltaTime]() {
				_generator->tick(DeltaTime);
			});

			IFlexGraphSimulation * flexTicker = Cast<IFlexGraphSimulation>(_generator);

			if (flexTicker)
			{
				_flexSimulation->addFlexTickWork([&, DeltaTime, flexTicker](
					float dt,
					NvFlexVector<int>& neighbourIndices,
					NvFlexVector<int>& neighbourCounts,
					NvFlexVector<int>& apiToInternal,
					NvFlexVector<int>& internalToAPI,
					int maxParticles)
				{
					flexTicker->flexTick(dt,
						neighbourIndices,
						neighbourCounts,
						apiToInternal,
						internalToAPI,
						maxParticles);
				});
			}

			_flexSimulation->tick(DeltaTime);
		}
		// This generator doesn't need our flex simulation
		else
		{
			_generator->tick(DeltaTime);
		}
	}
}

void UDiscreteElementEditorComponent::postTickPhysics(float DeltaTime)
{

}

void UDiscreteElementEditorComponent::setCurrentGenerator(UElementGenerator * generator)
{
	if (generator == _generator)
		return;

	if (_generator)
	{
		_generator->detach();
		_generator->flexSimulation = nullptr;
	}

	_generator = generator;

	if (_generator)
	{
		if (_generator->wantsFlex())
			_generator->flexSimulation = _flexSimulation.get();

		_generator->attach(context.get(), simulationManager);

		if (!_paused)
			_generator->start();
	}
}

UElementGenerator* UDiscreteElementEditorComponent::currentGenerator()
{
	return _generator;
}

AActor * UDiscreteElementEditorComponent::exemplarActor()
{
	if (URegionGrowingGenerator * rgc = generator<URegionGrowingGenerator>())
		return rgc->exemplar;
	else
		return nullptr;
}

void UDiscreteElementEditorComponent::stop()
{
	_paused = true;

	if( _generator )
		_generator->stop();
}

void UDiscreteElementEditorComponent::start()
{
	_paused = false;

	if (_generator)
		_generator->start();

	_flexSimulation->begin();
}

OutputDomainExport UDiscreteElementEditorComponent::exportOutputDomain()
{
	OutputDomainExport result;

	FTransform toWorld = GetOwner()->GetRootComponent()->GetComponentTransform();

	result.read(graph);

	// transform to world space

	for (FGraphNode& node : result.graph.allNodes)
	{
		if (!node.isValid())
			continue;

		node.position = toWorld.TransformPosition(node.position);
		node.orientation = toWorld.TransformRotation(node.orientation);
		node.scale = toWorld.GetScale3D().X * node.scale;
	}

	result.meshInterface = context->meshInterface;

	return result;
}

void UDiscreteElementEditorComponent::loadElementDomain(FGraphSnapshot& toLoad)
{
	auto cachcedGenerator = currentGenerator();

	if (cachcedGenerator)
		setCurrentGenerator(nullptr);

	toLoad.restore(graph);

	if (cachcedGenerator)
		setCurrentGenerator(cachcedGenerator);
}


void UDiscreteElementEditorComponent::_init()
{
	if (!simulationManager)
		return;
		
	context->meshInterface = std::make_shared<tcodsMeshInterface>();

	context->limits = limits;

	simulationManager->camera = camera;
	simulationManager->init(graph, *GetOwner(), false);

	_didInit = true;
}


void UDiscreteElementEditorComponent::RegisterComponentTickFunctions(bool bRegister)
{
	if (bRegister)
	{
		if (SetupActorComponentTickFunction(&PrePhysicsTick))
		{
			PrePhysicsTick.Target = this;
		}

		if (SetupActorComponentTickFunction(&PostPhysicsTick))
		{
			PostPhysicsTick.Target = this;
		}
	}
	else
	{
		if (PrePhysicsTick.IsTickFunctionRegistered())
		{
			PrePhysicsTick.UnRegisterTickFunction();
		}

		if (PostPhysicsTick.IsTickFunctionRegistered())
		{
			PostPhysicsTick.UnRegisterTickFunction();
		}
	}

	Super::RegisterComponentTickFunctions(bRegister);

}

// ------------------------------------------------------------
// Limits Drawing
// ------------------------------------------------------------

void UDiscreteElementEditorComponent::_updateLimits()
{
	if (!_limitsBoxComponent || drawLimits == false)
		return;

	FVector centre = _limitsBoxComponent->GetComponentLocation();
	FVector extents = _limitsBoxComponent->GetScaledBoxExtent();

	limits.Min = centre - extents;
	limits.Max = centre + extents;
}

void UDiscreteElementEditorComponent::_updateDrawLimits()
{
	if (drawLimits)
	{
		_drawLimits();
		_updateLimits();
	}
	else
		_hideLimits();
}

void UDiscreteElementEditorComponent::_drawLimits()
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

void UDiscreteElementEditorComponent::_hideLimits()
{
	if (_limitsBoxComponent)
		_limitsBoxComponent->SetVisibility(false);
}
