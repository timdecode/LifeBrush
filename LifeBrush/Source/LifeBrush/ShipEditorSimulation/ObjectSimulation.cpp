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

void UGraphSimulationManager::attachSimulations()
{
	for (auto simulation : simulations)
	{
		if (GetWorld()->WorldType != EWorldType::Editor)
			simulation->attach();
		else if (simulation->canRunInEditor())
			simulation->attach();
	}
}

void UGraphSimulationManager::detachSimulations()
{
	for (auto simulation : simulations)
	{
		if (GetWorld()->WorldType != EWorldType::Editor)
			simulation->detach();
		else if (simulation->canRunInEditor())
			simulation->detach();
	}
}

void UGraphSimulationManager::begin()
{
	for (auto simulation : simulations)
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
	for (auto engine : simulations)
	{
		engine->tick(deltaT);
	}
}

void UGraphSimulationManager::tick_paused(float deltaT)
{
	// tick the engines
	for (auto engine : simulations)
	{
		engine->tick_paused(deltaT);
	}
}

void UGraphSimulationManager::didPause()
{
	for (auto engine : simulations)
	{
		engine->didPause();
	}
}

void UGraphSimulationManager::didResume()
{
	for (auto engine : simulations)
	{
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
