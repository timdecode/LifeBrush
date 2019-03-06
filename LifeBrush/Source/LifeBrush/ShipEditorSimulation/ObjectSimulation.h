// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "Graph.h"

#include "ObjectSimulation.generated.h"

class UGraphSimulationManager;

/**
*
*/
UCLASS( BlueprintType, EditInlineNew, DefaultToInstanced, abstract )
class LIFEBRUSH_API UObjectSimulation : public UObject, public ComponentListener
{
	GENERATED_BODY()

public:
	FGraph * graph = nullptr;
	UGraphSimulationManager * simulationManager = nullptr;
	AActor * actor = nullptr;

public:
	virtual bool canRunInEditor()
	{
		return false;
	}
	
	// Called after the graph, simulationManager and actor have been set.
	virtual void attach()
	{
	}

	virtual void detach() {}

	// This will be called after `attach`, but before the first `tick` call.
	virtual void begin() {}

	virtual void tick( float deltaT )
	{

	}

	virtual void tick_paused(float deltaT)
	{

	}

	// The simulation did enter a pause state, the next tick will be tick_paused.
	virtual void didPause() {}

	// The simulation did enter a playing state, the next tick will be a tick.
	virtual void didResume() {}

	virtual void snapshotToActor(AActor * actor) 
	{

	}

	// Called by FGraph::clear
	virtual void clear()
	{

	}

	virtual void componentAdded( FGraphNodeHandle node, ComponentType type )
	{

	}

	virtual void componentRemoved( FGraphNodeHandle node, ComponentType type )
	{

	}

	virtual void componentUpdated( FGraphNodeHandle node, ComponentType type )
	{

	}
};

UCLASS(BlueprintType, DefaultToInstanced)
class LIFEBRUSH_API UGraphSimulationManager : public UObject
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, Instanced, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<class UObjectSimulation*> simulations;

	UCameraComponent * camera = nullptr;

protected:
	FGraph * _graph = nullptr;
	AActor * _actor = nullptr;

public:

	// init must be called before init_simulations or registerSimulation.
	void init(FGraph& graph_in, AActor& actor, bool autoInstantiateSimulations = false);
	// init_simulations must be called before tick or tick_paused.
	void attachSimulations();

	void detachSimulations();

	void begin();

	void tick(float deltaT);
	void tick_paused(float deltaT);

	void didPause();
	void didResume();

	FGraph * graph() { return _graph; }
	AActor * actor() { return _actor; }

	template<class TObjectSimulation>
	TObjectSimulation * simulation()
	{
		for (auto s : simulations)
		{
			if (s->GetClass() == TObjectSimulation::StaticClass())
			{
				return Cast<TObjectSimulation>(s);
			}
		}

		return nullptr;
	}


	//---------------------------------------------------------------------------
	// Simulation Registration
	//---------------------------------------------------------------------------
	template<class TObjectSimulation>
	TObjectSimulation* registerSimulation()
	{
		UClass * simulationClass = TObjectSimulation::StaticClass();

		return static_cast<TObjectSimulation*>(registerSimulation(simulationClass));
	}

	UObjectSimulation * registerSimulation(TSubclassOf<UObjectSimulation> simulationClass);

protected:

	void _autoInstantiateSimulations();


	//UGraphSimulationManager& operator=(const UGraphSimulationManager& other)
	//{
	//	simulations.Empty();
	//	for (UObjectSimulation * simulation : other.simulations)
	//	{
	//		UObjectSimulation * duplicateSimulation = DuplicateObject(simulation, nullptr);
	//		simulations.Add(duplicateSimulation);
	//	}

	//	return *this;
	//}
};



