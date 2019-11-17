// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "MeshFilamentSimulation.h"
#include "Visualization/EdgeFactory.h"

void UMeshFilamentSimulation::attach()
{
	if (!_lineFactory && actor)
	{
		_lineFactory = NewObject<UColoredLineFactory>(actor, TEXT("meshFilamentLineFactory"));

		if (defaultLineMaterial)
			_lineFactory->material = defaultLineMaterial;

		_lineFactory->RegisterComponent();
	}

	_dirty = true;

	graph->addEdgeObjectListener<FFilamentConnection>((EdgeObjectListener*)this);
}

void UMeshFilamentSimulation::detach()
{
	graph->removeEdgeObjectListener<FFilamentConnection>((EdgeObjectListener*)this);

	if (_lineFactory)
	{
		_lineFactory->DestroyComponent();
		_lineFactory = nullptr;
	}

	_dirty = true;
}

void UMeshFilamentSimulation::begin()
{
	_updateFilaments();
}

void UMeshFilamentSimulation::tick(float deltaT)
{
	if (_areFilamentsFrozen)
		_updateFrozenComponents();

	if (!_areFilamentsFrozen || _dirty)
		_updateFilaments();

	_updateDesaturated();
}

void UMeshFilamentSimulation::tick_paused(float deltaT)
{
	if (_areFilamentsFrozen)
		_updateFrozenComponents();
	
	if( !_areFilamentsFrozen || _dirty )
		_updateFilaments();

	_updateDesaturated();
}


void UMeshFilamentSimulation::_updateFilaments()
{
	// indexed by section
	TMap<int32, FColoredLineBuilder> buildersBySection;

	// preallocate
	for (auto pair : groupToSection)
	{
		auto& builder = buildersBySection.Emplace(pair.Value);
		builder.numCircleComponents = 8;
	}


	auto n = graph->privateEdges.Num();

	if( _dirty )
	{
		auto& storage = graph->rawEdgeStorage<FFilamentConnection>();

		// sort by group then segmentID
		storage.sort([&](FFilamentConnection& a, FFilamentConnection& b) {
			if (a.group == b.group)
				return a.segmentID < b.segmentID;
			else
				return a.group < b.group;
		});
	}

	auto filaments = graph->edgeView<FFilamentConnection>();


	FFilamentConnection * lastConnection = nullptr;
	FGraphEdge * lastEdge = nullptr;

	FColoredLineBuilder * curBuilder = nullptr;

	FGraphNodeHandle lastPoint;

	{
		TypedEdgeStorage<FFilamentConnection>& filamentStorage = graph->rawEdgeStorage<FFilamentConnection>();

		auto num = filamentStorage.totalSize();

		auto findNext = [&](int i) -> FFilamentConnection* {
			while (i < num && !filamentStorage.isValid(i) )
				++i;

			if (i < num)
				return &filamentStorage[i];
			else
				return nullptr;
		};

		FGraphNodeHandle lastHandle;

		for (int i = 0; i < num; ++i)
		{
			if( !filamentStorage.isValid(i) ) continue;
			
			FFilamentConnection * filament = &filamentStorage[i];
			FFilamentConnection * nextFilament = findNext(i + 1);

			FGraphEdge& edge = graph->edge(filamentStorage.edge(i));

			FGraphNodeHandle ha = FGraphNodeHandle(edge.a);
			FGraphNodeHandle hb = FGraphNodeHandle(edge.b);

			FGraphNode& a = graph->node(edge.a);
			FGraphNode& b = graph->node(edge.b);

			FVector dir = (b.position - a.position).GetSafeNormal();

			if ((lastConnection == nullptr || lastConnection->group != filament->group) &&
				filament->visible )
			{
				auto section = groupToSection[filament->group];
				curBuilder = &buildersBySection[section];

				FVector aPosition = a.position - dir * filament->aExtension;

				curBuilder->begin(aPosition, filament->radius);
				curBuilder->addPoint(b.position, filament->radius);
			}
			
			if ((!nextFilament || nextFilament->group != filament->group) &&
				filament->visible)
			{
				FVector bPosition = b.position + dir * filament->bExtension;

				curBuilder->addPoint(a.position, filament->radius);
				curBuilder->end(bPosition, filament->radius);
			}
			else if( filament->visible )
				curBuilder->addPoint(a.position, filament->radius);

			lastConnection = filament;
		}
	}

	//filaments.each([&](FFilamentConnection& filament, FGraphEdge& edge) {
	//	bool newGroup = !lastConnection || lastConnection->group != filament.group;

	//	if (newGroup && lastConnection)
	//	{
	//		FVector p = graph->node(lastEdge->b).position;

	//		if(curBuilder)
	//			curBuilder->end(p, lastConnection->radius);
	//	}

	//	FVector a = graph->node(edge.a).position;

	//	if (newGroup)
	//	{
	//		auto section = groupToSection[filament.group];
	//		curBuilder = &builders[section];

	//		curBuilder->begin(a, filament.radius);
	//	}
	//	else
	//		curBuilder->addPoint(a, filament.radius);

	//	FVector b = graph->node(edge.b).position;

	//	curBuilder->addPoint(b, filament.radius);


	//	lastConnection = &filament;
	//	lastEdge = &edge;
	//});


	if (graph->numEdges() > 0)
	{
		bool topologyChanged = graph->numEdges() != lastCount || _dirty;

		TMap<int32, UMaterialInterface*> sectionToMaterial;
		for (auto& pair : materialToSection)
			sectionToMaterial.Add(pair.Value, pair.Key);

		for (auto& pair : buildersBySection)
		{
			auto section = pair.Key;
			auto material = sectionToMaterial[section];

			_lineFactory->commitWithFastPathOption(pair.Value, section, material, topologyChanged);
		}
	}

	_dirty = false;

	lastCount = graph->numEdges();
}

void UMeshFilamentSimulation::_updateFrozenComponents()
{
	auto& forzens = graph->componentStorage<FFrozenFilament>();

	for (FFrozenFilament& frozen : forzens)
	{
		if (!frozen.isValid()) continue;

		FGraphNode& node = graph->node(frozen.nodeHandle());

		node.position = frozen.position;
		node.orientation = frozen.orientation;
	}
}

int32 UMeshFilamentSimulation::addSection(UMaterialInterface * material)
{
	if (materialToSection.Contains(material))
		return materialToSection[material];

	int32 section = 0;

	for (auto& pair : materialToSection)
	{
		if (pair.Value >= section)
			section = pair.Value + 1;
	}

	materialToSection.Add(material, section);

	// we added a section, things are so dirty now
	_dirty = true;

	return section;
}


uint32 UMeshFilamentSimulation::nextGroup(UMaterialInterface * material)
{
	if (!material) material = defaultLineMaterial;

	auto theGroup = _nextGroup;
	_nextGroup++;

	int32 section = addSection(material);
	groupToSection.Add(theGroup, section);

	return theGroup;
}

void UMeshFilamentSimulation::edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type)
{
	const static auto FilmentConnectionType = edgeType<FFilamentConnection>();

	if (FilmentConnectionType != type)
		return;

	_dirty = true;
}

void UMeshFilamentSimulation::snapshotToActor(AActor * actor)
{
	URuntimeMeshComponent * rmc = _lineFactory->_runtimeMeshComponent;

	if (rmc)
	{
		URuntimeMeshComponent * newRMC = NewObject<URuntimeMeshComponent>(actor);

		newRMC->SetShouldSerializeMeshData(true);

		int32 newRMCSection = 0;

		auto sectionIds = rmc->GetRuntimeMesh()->GetSectionIds();

		for (auto& pair : materialToSection)
		{
			int32 section = pair.Value;
			UMaterialInterface * material = pair.Key;

			if (!rmc->DoesSectionExist(section))
				continue;

			auto readonlySection = rmc->GetSectionReadonly(section);

			if (!readonlySection || readonlySection->NumIndices() == 0)
				continue;

			auto builder = MakeRuntimeMeshBuilder(*readonlySection);

			readonlySection->CopyTo(builder, true);

			newRMC->CreateMeshSection(newRMCSection, builder, false, EUpdateFrequency::Infrequent);
			newRMC->SetSectionMaterial(newRMCSection, material);
			newRMC->SetMaterial(newRMCSection, material);

			newRMCSection++;
		}

		if (isDesaturated)
		{
			newRMC->SetCustomDepthStencilValue(1);
			newRMC->SetRenderCustomDepth(true);
		}

		newRMC->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		actor->AddInstanceComponent(newRMC);

		newRMC->RegisterComponent();

	}
}

void UMeshFilamentSimulation::freezeFilaments()
{
	TypedEdgeStorage<FFilamentConnection>& filamentStorage = graph->rawEdgeStorage<FFilamentConnection>();

	int num = filamentStorage.totalSize();

	for (int i = 0; i < num; ++i)
	{
		if (!filamentStorage.isValid(i)) continue;

		FFilamentConnection * filament = &filamentStorage[i];

		FGraphEdge& edge = graph->edge(filamentStorage.edge(i));

		FGraphNode& aNode = graph->node(edge.a);
		FGraphNode& bNode = graph->node(edge.b);

		FFrozenFilament& aFrozen = aNode.addComponent<FFrozenFilament>(*graph);
		aFrozen.position = aNode.position;
		aFrozen.orientation = aNode.orientation;

		FFrozenFilament& bFrozen = aNode.addComponent<FFrozenFilament>(*graph);
		bFrozen.position = bNode.position;
		bFrozen.orientation = bNode.orientation;
	}

	_areFilamentsFrozen = true;
}

void UMeshFilamentSimulation::unfreezeFilaments()
{
	// work on a copy, because we are mutating
	auto copiedFrozens = graph->componentStorage<FFrozenFilament>();

	for (FFrozenFilament& frozen : copiedFrozens)
	{
		if (!frozen.isValid()) continue;

		FGraphNode& node = graph->node(frozen.nodeHandle());

		node.removeComponent<FFrozenFilament>(*graph);
	}
}

bool UMeshFilamentSimulation::areFilamentsFrozen()
{
	return _areFilamentsFrozen;
}

void UMeshFilamentSimulation::setDesaturated(bool desaturated)
{
	isDesaturated = desaturated;


	_updateDesaturated();
}


void UMeshFilamentSimulation::_updateDesaturated()
{
	if (_lineFactory && _lineFactory->_runtimeMeshComponent)
	{
		if (isDesaturated)
		{
			_lineFactory->_runtimeMeshComponent->SetRenderCustomDepth(true);
			_lineFactory->_runtimeMeshComponent->SetCustomDepthStencilValue(1); // 1 is the magic value for the desaturation post process volume
		}
		else
		{
			_lineFactory->_runtimeMeshComponent->SetRenderCustomDepth(false);
			_lineFactory->_runtimeMeshComponent->SetCustomDepthStencilValue(0); // 1 is the magic value for the desaturation post process volume
		}
	}
}