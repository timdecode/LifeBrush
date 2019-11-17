// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "ObjectSimulation.h"

void UGraphSimulationManager::init(FGraph& graph_in, AActor& actor, bool autoInstantiateSimulations /*= false*/)
{
	this->_graph = &graph_in;
	_actor = &actor;

	if (autoInstantiateSimulations)
		_autoInstantiateSimulations();

	for (auto simulation : simulations)
	{
		simulation->graph = _graph;
		simulation->actor = _actor;
		simulation->simulationManager = this;
	}
}

void UGraphSimulationManager::attachAllSimulations()
{
	for (auto simulation : simulations)
		_attachSimulation(simulation);
}

void UGraphSimulationManager::detachAllSimulations()
{
	auto cachedAttachedSimulations = _attachedSimulations;

	for (auto simulation : cachedAttachedSimulations)
	{
		_detachSimulation(simulation);
	}
}

void UGraphSimulationManager::begin()
{
	for (auto simulation : _attachedSimulations)
	{
		simulation->begin();
	}
}

void UGraphSimulationManager::_autoInstantiateSimulations()
{
	// nuke classes that might not be there anymore after code changes
	simulations.Remove(nullptr);

	for (TObjectIterator<UClass> it; it; ++it)
	{
		if (!it->IsChildOf(UObjectSimulation::StaticClass()) || *it == UObjectSimulation::StaticClass())
			continue;

		registerSimulation(*it);
	}
}






void UGraphSimulationManager::tick(float deltaT)
{
	_graph->tickCount++;

	// tick the engines
	for (auto engine : _attachedSimulations)
	{
		engine->tick(deltaT);
	}
}

void UGraphSimulationManager::tick_paused(float deltaT)
{
	// tick the engines
	for (auto engine : _attachedSimulations)
	{
		if( engine->simulationManager )
			engine->tick_paused(deltaT);
	}
}

void UGraphSimulationManager::didPause()
{
	for (auto engine : _attachedSimulations)
	{
		if (engine->simulationManager)
			engine->didPause();
	}
}

void UGraphSimulationManager::didResume()
{
	for (auto engine : _attachedSimulations)
	{
		if (engine->simulationManager)
			engine->didResume();
	}
}

UObjectSimulation * UGraphSimulationManager::registerSimulation(TSubclassOf<UObjectSimulation> simulationClass)
{
	UObjectSimulation * simulation = nullptr;

	for (auto s : simulations)
	{
		if (s->GetClass() == simulationClass)
		{
			simulation = s;
			break;
		}
	}

	// we need to create a simulation
	if (simulation == nullptr)
	{
		simulation = NewObject<UObjectSimulation>(this, *simulationClass);
		simulations.Add(simulation);
	}

	simulation->graph = _graph;
	simulation->actor = _actor;
	simulation->simulationManager = this;

	return simulation;
}

bool UGraphSimulationManager::isAttached(UObjectSimulation * simulation)
{
	return _attachedSimulations.Contains(simulation);
}

void UGraphSimulationManager::_attachSimulation(UObjectSimulation * simulation)
{
	checkf(simulations.Contains(simulation), TEXT("The simulation was not in the manager's simulations. The simulation must be registered first."));

	// only attach simulations that can run in the editor
	// or attach if we are in game mode
	if ((GetWorld()->WorldType != EWorldType::Editor || simulation->canRunInEditor()) &&
		!_attachedSimulations.Contains(simulation))
	{
		simulation->graph = _graph;
		simulation->actor = _actor;
		simulation->simulationManager = this;

		simulation->attach();

		_attachedSimulations.Add(simulation);
	}
}

void UGraphSimulationManager::_detachSimulation(UObjectSimulation * simulation)
{
	checkf(simulations.Contains(simulation), TEXT("The simulation was not in the manager's simulations. The simulation must be registered first."));

	if(	_attachedSimulations.Contains(simulation))
	{
		simulation->graph = _graph;
		simulation->actor = _actor;
		simulation->simulationManager = this;

		simulation->detach();

		_attachedSimulations.Remove(simulation);
	}
}
