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
struct OutputDomainExport;
struct FGraphSnapshot;
class UChunkedVolumeComponent;

UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UDiscreteElementEditorComponent : public UGraphComponent
{
	GENERATED_BODY()

public:
	std::unique_ptr<SynthesisContext> context;

	UCameraComponent * camera;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FBox limits = FBox(EForceInit::ForceInitToZero);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool drawLimits = false;

	UPROPERTY(EditAnywhere, Instanced, BlueprintReadWrite, Category = "ShipEditor")
	TArray<UElementGenerator*> generators;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FNvFlexParameters flexParams;

protected:
	UElementGenerator * _generator;

	bool _didInit = false;

	bool _paused = false;

	std::unique_ptr<FFlexSimulation> _flexSimulation;

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

	AActor * exemplarActor();

	void stop();
	void start();

	OutputDomainExport exportOutputDomain();

	void loadElementDomain(FGraphSnapshot& toLoad);


protected:
	void _initContextAndGraph();

	void _init();

	virtual void RegisterComponentTickFunctions(bool bRegister) override;

protected:
	// Limits Drawing
	// --------------
	void _updateLimits();
	void _updateDrawLimits();
	void _drawLimits();
	void _hideLimits();

	UPROPERTY(EditAnywhere, BlueprintReadOnly)
	UBoxComponent * _limitsBoxComponent = nullptr;
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

	return nullptr;
}
