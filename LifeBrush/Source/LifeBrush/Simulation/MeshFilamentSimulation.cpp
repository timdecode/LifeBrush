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

	_needsSort = true;

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

	_needsSort = true;
}

void UMeshFilamentSimulation::begin()
{
	_updateFilaments();
}

void UMeshFilamentSimulation::tick(float deltaT)
{
	_updateFilaments();
}

void UMeshFilamentSimulation::_updateFilaments()
{
	FColoredLineBuilder builder;

	builder.numCircleComponents = 8;

	auto n = graph->privateEdges.Num();

	if( _needsSort )
	{
		auto& storage = graph->rawEdgeStorage<FFilamentConnection>();

		// sort by group then segmentID
		storage.sort([&](FFilamentConnection& a, FFilamentConnection& b) {
			if (a.group == b.group)
				return a.segmentID < b.segmentID;
			else
				return a.group < b.group;
		});

		_needsSort = false;
	}

	auto filaments = graph->edgeView<FFilamentConnection>();


	FFilamentConnection * lastConnection = nullptr;
	FGraphEdge * lastEdge = nullptr;

	filaments.each([&](FFilamentConnection& filament, FGraphEdge& edge) {
		bool newGroup = !lastConnection || lastConnection->group != filament.group;

		if (newGroup && lastConnection)
		{
			FVector p =graph->node(lastEdge->b).position;

			builder.end(p, lastConnection->radius);
		}

		FVector a = graph->node(edge.a).position;

		if (newGroup)
		{
			builder.begin(a, filament.radius);
		}
		else
			builder.addPoint(a, filament.radius);

		FVector b = graph->node(edge.b).position;

		builder.addPoint(b, filament.radius);


		lastConnection = &filament;
		lastEdge = &edge;
	});


	if (graph->numEdges() > 0)
	{
		bool topologyChanged = graph->numEdges() != lastCount;

		_lineFactory->commitWithFastPathOption(builder, 0, defaultLineMaterial, topologyChanged);
	}

	lastCount = graph->numEdges();
}

void UMeshFilamentSimulation::edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type)
{
	const static auto FilmentConnectionType = edgeType<FFilamentConnection>();

	if (FilmentConnectionType != type)
		return;

	_needsSort = true;
}
