//
//  DiscreteElementEditorComponent.h
//  LifeBrush
//
//  Created by Timothy Davison on 2018-12-21.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "Algorithm/SynthesisContext.h"
#include "ShipEditorSimulation/ObjectSimulation.h"
#include "Simulation/FlexElements.h"


#include "DiscreteElementEditorComponent.generated.h"

class UElementGenerator;
class UFlexSimulationComponent;
struct OutputDomainExport;
struct FGraphSnapshot;
class UChunkedVolumeComponent;
struct FFlexSimulation;

// Manages UElementGenerators that interface with a UFlexSimulationComponent.
// 
// This component requires a UFlexSimulationComponent sibling.
UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UDiscreteElementEditorComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, Instanced, BlueprintReadWrite, Category = "LifeBrush")
	TArray<UElementGenerator*> generators;

	// The host actor for the UFlexSimulationComponent. If this is not set, we'll look at our own
	// parent to find one. The VRSketchyPawn will set this if it has a reference to this component's actor.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	AActor * simulationActor = nullptr;

	// This will be set by the VRSketchyPawn.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	AActor * exemplarActor = nullptr;

protected:
	UElementGenerator * _generator;

	bool _didInit = false;

	bool _paused = false;


	UPROPERTY(EditDefaultsOnly, Category="ComponentTick")
	struct FActorComponentTickFunction PrePhysicsTick;

	UPROPERTY(EditDefaultsOnly, Category = "ComponentTick")
	struct FActorComponentTickFunction PostPhysicsTick;



public:
	UDiscreteElementEditorComponent();

	virtual void InitializeComponent() override;
	virtual void BeginPlay() override;

	virtual void TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction) override;

	void preTickPhysics(float DeltaTime);
	void postTickPhysics(float DeltaTime);

	void setCurrentGenerator(UElementGenerator * generator);
	UElementGenerator* currentGenerator();

	template<typename GeneratorType>
	GeneratorType* generator();

	UFlexSimulationComponent * flexSimulationComponent();
	FFlexSimulation * flexSimulation();

	SynthesisContext * context();

	void stop();
	void start();

	OutputDomainExport exportOutputDomain();

	void loadElementDomain(FGraphSnapshot& toLoad);

	void init() { _init();  }

protected:
	void _init();

	void _autoInstantiateGenerators();

	void _initDependencies(UElementGenerator* generator,
		TSet<UElementGenerator*>& initialized, 
		TMap<UClass*, UElementGenerator*>& classToGenerator);

	UElementGenerator * _registerGenerator(TSubclassOf<UElementGenerator> generatorClass);

	virtual void RegisterComponentTickFunctions(bool bRegister) override;

protected:



};

template<typename GeneratorType>
GeneratorType* UDiscreteElementEditorComponent::generator()
{
	UClass * generatorClass = GeneratorType::StaticClass();

	for (auto& e : generators)
	{
		if (e && e->GetClass() == generatorClass)
			return static_cast<GeneratorType*>(e);
	}

	return static_cast<GeneratorType*>(_registerGenerator(generatorClass));
}
