// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "GraphSimulationProbe.h"
#include "Synapse/Synapse.h"

TMap<UScriptStruct*, uint32> UGraphSimulationProbeComponent::generateCounts()
{
	TMap<UScriptStruct*, uint32> result;

	if (!flexSimulation)
		return result;

	TMap<StructTypeManager<FGraphObject>::TypeId, uint32> typeCounts;
	
	FBox probe = getProbeBox();




	FGraph& graph = flexSimulation->graph;



	TArray<FGraphNodeHandle> toProbe;

	for (FGraphNode& node : graph.allNodes )
	{
		if( !node.isValid() ) continue;

		if (!probe.IsInside(node.position)) continue; 

		for (ComponentType theType : node.components)
		{
			typeCounts.FindOrAdd(theType)++;
		}
	}

	// convert
	for (auto& pair : typeCounts)
	{
		UScriptStruct * theStruct = FGraphObject::componentStruct(pair.Key);

		result.Add(theStruct, pair.Value);
	}

	return result;
}

FBox UGraphSimulationProbeComponent::getProbeBox()
{
	FTransform componentToWorld = this->GetComponentTransform();

	FBoxSphereBounds probe = this->CalcBounds(componentToWorld);

	// and then to simulation space
	FTransform simulationTransform = flexSimulation->GetOwner()->GetActorTransform().Inverse();
	probe.TransformBy(simulationTransform);

	return probe.GetBox();
}

AGraphSimulationProbe::AGraphSimulationProbe()
{
	PrimaryActorTick.bCanEverTick = true;
}

void AGraphSimulationProbe::BeginPlay()
{
	Super::BeginPlay();

	_prepareFilePath();

	// get our flex components
	if( flexSimulationActor )
		flexSimulation = flexSimulationActor->FindComponentByClass<UFlexSimulationComponent>();

	// get our probes
	USceneComponent * sceneComponent = this->GetRootComponent();

	TArray<USceneComponent*> children;

	sceneComponent->GetChildrenComponents(true, children);

	for (USceneComponent * child : children)
	{
		UGraphSimulationProbeComponent * probe = Cast<UGraphSimulationProbeComponent>(child);

		if (!probe) continue;

		probe->flexSimulation = flexSimulation;

		probeGroups.FindOrAdd(probe->groupName).Add(probe);
	}

	if (flexSimulation)
	{
		if (monitoredTypes.Num() == 0)
		{
			FGraph& graph = flexSimulation->graph;

			for (auto& storage : graph._componentStorage)
			{
				UScriptStruct * theStruct = storage.componentClass;

				_runtimeMonitoredTypes.Add(theStruct);
			}
		}
		else
		{
			for (FString& theName : monitoredTypes)
			{
				UScriptStruct * theStruct = FindObject<UScriptStruct>(ANY_PACKAGE, *theName);


				_runtimeMonitoredTypes.Add(theStruct);

			}
		}
	}

}


void AGraphSimulationProbe::Tick(float IgnoredDeltaT)
{
	Super::Tick(IgnoredDeltaT);

	if (!flexSimulation)
		return;

	if (!flexSimulation->flexSimulation()->isPlaying())
		return;

	if (firstTick)
	{
		for (auto& pair : probeGroups)
		{
			_writeHeader(pair.Key.ToString());
		}

		_writeChannelStatesHeader();
		_writeOneChannelHeader();

		firstTick = false;
	}

	float deltaT = flexSimulation->timeStep;

	timeLeft -= deltaT;


	if (timeLeft < 0.0f)
	{
		timeLeft = frequency;

		// get our data
		for (auto& pair : probeGroups)
		{
			TMap<UScriptStruct*, uint32> counts;

			for (UGraphSimulationProbeComponent * probe : pair.Value)
			{
				auto probeCounts = probe->generateCounts();

				for (auto& count : probeCounts)
					counts.FindOrAdd(count.Key) += count.Value;
			}

			FString groupName = pair.Key.ToString();

			_writeCounts(counts, groupName);
		}

		// write our channel states
		_writeChannelStates();
		_writeOneChannel();
	}
}

void AGraphSimulationProbe::_writeHeader(FString groupName)
{
	FGraph& graph = flexSimulation->graph;

	FString line = "";

	bool firstStorage = true;
	for (auto& storage : graph._componentStorage)
	{
		UScriptStruct * theStruct = storage.componentClass;

		if(!_runtimeMonitoredTypes.Contains(theStruct)) continue;

		if (firstStorage)
			firstStorage = false;
		else
			line += separator;

		line.Append(theStruct->GetName());
	}

	line += LINE_TERMINATOR;

	FString filePath = _filePath(groupName);

	_appendToFile(filePath, line);
}

void AGraphSimulationProbe::_appendToFile(FString filePath, FString theString)
{
	FFileHelper::SaveStringToFile(
		theString,
		*filePath,
		FFileHelper::EEncodingOptions::AutoDetect,
		&IFileManager::Get(),
		EFileWrite::FILEWRITE_Append);

}

void AGraphSimulationProbe::_writeCounts(TMap<UScriptStruct*, uint32>& counts, FString groupName)
{
	FString line = "";

	FGraph& graph = flexSimulation->graph;

	bool firstStorage = true;
	for (auto& storage : graph._componentStorage)
	{
		UScriptStruct * theStruct = storage.componentClass;

		if (!_runtimeMonitoredTypes.Contains(theStruct)) continue;

		if (firstStorage)
			firstStorage = false;
		else
			line += separator;

		auto count = counts.Find(theStruct);

		if (count)
			line.AppendInt(*count);
		else
			line.AppendInt(0);
	}

	line += LINE_TERMINATOR;

	FString filePath = _filePath(groupName);

	_appendToFile(filePath, line);
}

void AGraphSimulationProbe::_writeChannelStates()
{
	FGraph& graph = flexSimulation->graph;

	auto& voltageGatedSodiums = graph.componentStorage<FVoltageGatedIonChannel>();

	struct ChannelCounts
	{
		float averageVoltage = 0.0f;
		float averageDeltaVoltage = 0.0f;

		int numOpen = 0;
		int numClosed = 0;
		int numInactived = 0;
	} counts;

	int count = 0;
	// track one channel individually (we'll break out of this loop)
	for (FVoltageGatedIonChannel& gate : voltageGatedSodiums)
	{
		if (!gate.isValid()) continue;

		counts.averageVoltage += gate.averageValues.average();
		counts.averageDeltaVoltage += gate.deltaVoltage;

		counts.numOpen += gate.state == EVoltageGatedIonChannelState::Open ? 1 : 0;
		counts.numClosed += gate.state == EVoltageGatedIonChannelState::Closed ? 1 : 0;
		counts.numInactived += gate.state == EVoltageGatedIonChannelState::Inactivated ? 1 : 0;

		count++;
	}

	// compute averages
	if (count > 0)
	{
		counts.averageVoltage /= float(count);
		counts.averageDeltaVoltage /= float(count);
	}

	FString line = FString::SanitizeFloat(counts.averageVoltage) + separator +
		FString::SanitizeFloat(counts.averageDeltaVoltage) + separator +
		FString::FromInt(counts.numOpen) + separator +
		FString::FromInt(counts.numClosed) + separator +
		FString::FromInt(counts.numInactived) + LINE_TERMINATOR;

	FString fileName = "Channels";
	FString filePath = _filePath(fileName);

	_appendToFile(filePath, line);
}

void AGraphSimulationProbe::_writeChannelStatesHeader()
{
	FString fileName = "Channels";
	FString filePath = _filePath(fileName);

	FString line = "averageVoltate, averageDeltaVoltage, numOpen, numClosed, numInactivated";
	line += LINE_TERMINATOR;

	_appendToFile(filePath, line);
}

void AGraphSimulationProbe::_writeOneChannel()
{
	FGraph& graph = flexSimulation->graph;

	auto& voltageGatedSodiums = graph.componentStorage<FVoltageGatedIonChannel>();

	// track one channel individually (we'll break out of this loop)
	for (FVoltageGatedIonChannel& gate : voltageGatedSodiums)
	{
		if (!gate.isValid()) continue;

		FString stateAsString = "";
		{
			const UEnum * AlgorithmEnum = FindObject<UEnum>(ANY_PACKAGE, TEXT("EVoltageGatedIonChannelState"));
			int32 int32_algorithm = (int32)gate.state;

			stateAsString = AlgorithmEnum->GetDisplayNameTextByValue(int32_algorithm).ToString();
		}

		FString line = FString::SanitizeFloat(gate.averageValues.average()) + separator + // average voltage
			FString::SanitizeFloat(gate.deltaVoltage) + separator +// deltaVoltage
			FString::FromInt(gate.upCount) + separator +// upcount
			FString::FromInt(gate.downCount) + separator +// downcount
			stateAsString + LINE_TERMINATOR; // state (open/closed)

		FString fileName = "IndividualChannel";
		FString filePath = _filePath(fileName);

		_appendToFile(filePath, line);

		// done
		return;
	}
}



void AGraphSimulationProbe::_writeOneChannelHeader()
{
	FString fileName = "IndividualChannel";
	FString filePath = _filePath(fileName);

	FString line = "averageVoltage, deltaVoltage, upCount, downCount, state";
	line += LINE_TERMINATOR;

	_appendToFile(filePath, line);
}

FString AGraphSimulationProbe::_filePath(FString groupName)
{
	FString filePath = baseFileName + " " + groupName + ".txt";

	return filePath;
}

void AGraphSimulationProbe::_prepareFilePath()
{
	FString probeDir = FPaths::ProjectSavedDir() + "/GraphProbes";

	IPlatformFile& platformFile = FPlatformFileManager::Get().GetPlatformFile();

	platformFile.CreateDirectoryTree(*probeDir);

	FDateTime now = FDateTime::Now();

	baseFileName = probeDir + "/" + now.ToString();
}
