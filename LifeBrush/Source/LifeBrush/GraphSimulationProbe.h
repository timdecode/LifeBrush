// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"


#include "GraphSimulationProbe.generated.h"	



// The static mesh should be a sphere.
UCLASS(ClassGroup = (Custom), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UGraphSimulationProbeComponent : public UStaticMeshComponent
{
	GENERATED_BODY()
public:
	UFlexSimulationComponent * flexSimulation;

	// Graph components with the same name will have their results merged when under the same AGraphSimulationProbe actor.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	FName groupName;

public:
	TMap<UScriptStruct*, uint32> generateCounts();

	// In world-space.
	FBox getProbeBox();
};

UCLASS()
class LIFEBRUSH_API AGraphSimulationProbe : public AStaticMeshActor
{
	GENERATED_BODY()

public:
	AGraphSimulationProbe();

	virtual void BeginPlay() override;
	virtual void Tick(float DeltaSeconds) override;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	AActor * flexSimulationActor;

	// Frequency of writes in seconds.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	float frequency = 1.0f;

	//// If this is empty, it will be populated with all types
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	TArray<FString> monitoredTypes;



	//// Filter monitored types to not include these structs (when monitoredTypes is auto-filled)
	//UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Agent Library")
	//TArray<UScriptStruct*> filteredTypes;

protected:
	TSet<UScriptStruct*> _runtimeMonitoredTypes;

	UFlexSimulationComponent * flexSimulation;

	FString baseFileName;

	float timeLeft = 0.0f;

	bool firstTick = true;

	FString separator = ", ";

	TMap < FName, TArray<UGraphSimulationProbeComponent*> > probeGroups;

protected:
	void _writeHeader(FString groupName);
	void _appendToFile(FString filePath, FString theString);
	void _writeCounts(TMap<UScriptStruct*, uint32>& counts, FString groupName);

	void _writeChannelStates();
	void _writeOneChannel();

	void _writeChannelStatesHeader();
	void _writeOneChannelHeader();

	FString _filePath(FString groupName);

	void _prepareFilePath();
};