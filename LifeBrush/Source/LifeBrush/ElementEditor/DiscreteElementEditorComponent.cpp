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

	if (!_didInit)
		_init();
}

void UDiscreteElementEditorComponent::BeginPlay()
{
	Super::BeginPlay();
}

void UDiscreteElementEditorComponent::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);


	if (ThisTickFunction == &PrimaryComponentTick)
		preTickPhysics(DeltaTime);
	else if (ThisTickFunction == &PostPhysicsTick)
		postTickPhysics(DeltaTime);
}

void UDiscreteElementEditorComponent::preTickPhysics(float DeltaTime)
{

#if WITH_EDITOR
	// don't tick, in the editor, past this point
	bool inEditor = GetWorld()->WorldType == EWorldType::Editor;

	if (inEditor)
		return;
#endif

	if (!_generator || !flexSimulation() )
		return;


	// The generator wants a flex simulation wrapped around it.
	// So, add its work to the tick-work queue and tick flex.
	if (_generator->wantsFlex() && flexSimulation()->isPlaying())
	{
		flexSimulation()->addTickWork([&, DeltaTime]() {
			_generator->tick(DeltaTime);
		});
	}
	// This generator doesn't need our flex simulation
	else
	{
		_generator->tick(DeltaTime);
	}

	if (!flexSimulation()->isPlaying())
	{
		_generator->tickPaused(DeltaTime);
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
	}

	_generator = generator;

	if (_generator)
	{
		_generator->flexSimulation = flexSimulation();
		_generator->flexSimulationComponent = flexSimulationComponent();
		_generator->elementEditorComponent = this;

		_generator->attach(context(), &flexSimulation()->simulationManager);

		if (!_paused)
			_generator->start();
	}
}

UElementGenerator* UDiscreteElementEditorComponent::currentGenerator()
{
	return _generator;
}

UFlexSimulationComponent * UDiscreteElementEditorComponent::flexSimulationComponent()
{
	AActor * host = simulationActor ? simulationActor : GetOwner();

	UFlexSimulationComponent * flexComponent = host->FindComponentByClass<UFlexSimulationComponent>();

	return flexComponent;
}

FFlexSimulation * UDiscreteElementEditorComponent::flexSimulation()
{
	AActor * host = simulationActor ? simulationActor : GetOwner();

	UFlexSimulationComponent * flexComponent = host->FindComponentByClass<UFlexSimulationComponent>();

	if (flexComponent)
		return flexComponent->flexSimulation();
	else
		return nullptr;
}

SynthesisContext * UDiscreteElementEditorComponent::context()
{
	AActor * host = simulationActor ? simulationActor : GetOwner();

	UFlexSimulationComponent * flexComponent = host->FindComponentByClass<UFlexSimulationComponent>();

	if (flexComponent)
		return flexComponent->context();
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

	flexSimulation()->begin();
}

OutputDomainExport UDiscreteElementEditorComponent::exportOutputDomain()
{
	OutputDomainExport result;

	FTransform toWorld = GetOwner()->GetRootComponent()->GetComponentTransform();

	//result.read(graph);

	// transform to world space

	for (FGraphNode& node : result.graph.allNodes)
	{
		if (!node.isValid())
			continue;

		node.position = toWorld.TransformPosition(node.position);
		node.orientation = toWorld.TransformRotation(node.orientation);
		node.scale = toWorld.GetScale3D().X * node.scale;
	}

	result.meshInterface = context()->meshInterface;

	return result;
}

void UDiscreteElementEditorComponent::loadElementDomain(FGraphSnapshot& toLoad)
{
	auto cachcedGenerator = currentGenerator();

	if (cachcedGenerator)
		setCurrentGenerator(nullptr);

	//toLoad.restore(graph);

	if (cachcedGenerator)
		setCurrentGenerator(cachcedGenerator);
}


void UDiscreteElementEditorComponent::_init()
{
	if (!flexSimulation()) return;

	_autoInstantiateGenerators();

	for (auto gen : generators)
	{
		gen->elementEditorComponent = this;
		gen->flexSimulation = flexSimulation();
		gen->flexSimulationComponent = flexSimulationComponent();
		gen->exemplarActor = exemplarActor;
	}

	TMap<UClass*, UElementGenerator*> structToGenerator;

	for (auto gen : generators)
	{
		auto theClass = gen->GetClass();

		structToGenerator.Add(theClass, gen);
	}

	TSet<UElementGenerator*> initialized;

	for (auto gen : generators)
	{
		_initDependencies(gen, initialized, structToGenerator);
	}


	_didInit = true;
}


void UDiscreteElementEditorComponent::_autoInstantiateGenerators()
{
	// nuke classes that might not be there anymore (for example after code changes)
	generators.Remove(nullptr);

	for (TObjectIterator<UClass> it; it; ++it)
	{
		if (!it->IsChildOf(UElementGenerator::StaticClass()) || *it == UElementGenerator::StaticClass())
			continue;

		_registerGenerator(*it);
	}
}

UElementGenerator * UDiscreteElementEditorComponent::_registerGenerator(TSubclassOf<UElementGenerator> generatorClass)
{
	// see if we already have the generator
	for (auto generator : generators)
	{
		if (generator->GetClass() == generatorClass)
			return generator;
	}

	// create a new one
	UElementGenerator * generator = NewObject<UElementGenerator>(this, generatorClass);

	generators.Add(generator);

	return generator;
}

void UDiscreteElementEditorComponent::_initDependencies(
	UElementGenerator* generator, 
	TSet<UElementGenerator*>& initialized, 
	TMap<UClass*, UElementGenerator*>& classToGenerator)
{
	if (initialized.Contains(generator))
		return;

	auto dependencies = generator->dependencies();

	for (UClass * dependencyClass : dependencies)
	{
		auto dependency = classToGenerator.Find(dependencyClass);

		checkf(dependency != nullptr, TEXT("We don't have this dependency registered."));

		_initDependencies(*dependency, initialized, classToGenerator);
	}

	generator->init(context(), &flexSimulation()->simulationManager);
	initialized.Add(generator);
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

