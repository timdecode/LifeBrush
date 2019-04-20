//
//  Created by Timothy Davison on 2018-12-28.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ShipEditorSimulation/MeshSimulation.h"

#include "Simulation/FlexElements.h"

#include "StringGenerator.h"

void UStringGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager_hack)
{
	Super::attach(context, simulationManager_hack);

	_context = context;

	_simulationManager = simulationManager_hack;

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	_simulationManager->attachSimulations();

	_initPath();
}

void UStringGenerator::detach()
{
	Super::detach();

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	meshSim->detach();

	_initPath();
}

void UStringGenerator::tick(float deltaT)
{
	if (!_context)
		return;

	_tick(deltaT);
}

void UStringGenerator::_tick(float deltaT)
{
	size_t n = _brushPoints.size();

	if (n < 2)
		return;

	FGraph& graph = _context->graph();

	graph.beginTransaction();

	float scale = radius * scaleFactor;

	for (int i = 0; i + 1 < n; ++i)
	{
		auto& cur = _brushPoints[i];
		auto& next = _brushPoints[i + 1];

		FVector pCur = unreal(cur.position);
		FVector pNext = unreal(next.position);

		FVector dir = (pNext - pCur);
		float d = dir.Size();

		dir /= d;

		while (_t <= d)
		{
			FVector p = pCur + dir * _t;

			FGraphNodeHandle handle(graph.addNode(p, FQuat::Identity, scale));

			FGraphNode& node = _context->graph().node(handle);

			FGraphMesh& mesh = node.addComponent<FGraphMesh>(graph);

			mesh.material = material;
			mesh.staticMesh = staticMesh;

			FElementObject& element = node.addComponent<FElementObject>(graph);

			element.radius = radius;
			element.type = -1;

			FFlexParticleObject& particle = node.addComponent<FFlexParticleObject>(graph);

			FVelocityGraphObject& velocity = node.addComponent<FVelocityGraphObject>(graph);


			if (_lastStringElement)
			{
				_context->graph().connectNodes<FFlexConnection>(handle, _lastStringElement, 2.0f * radius, 0.5f);
			}

			_lastStringElement = handle;

			_t += 2.0f * radius;
		}

		_t -= d;
	}

	_context->graph().endTransaction();

	// nuke the brush points we visited, but keep the last one
	_brushPoints.erase(_brushPoints.begin(), _brushPoints.begin() + n - 1);
}

void UStringGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex)
{
	_initPath();
}

void UStringGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	Eigen::Vector3f eigen_point = eigen(point);

	if (!_brushPoints.empty())
	{
		Eigen::Vector3f last = _brushPoints.back().position;
		
		float reject = std::pow(radius * 0.25f, 2.0f);

		if ((last - eigen_point).squaredNorm() < reject)
			return;
	}

	_brushPoints.emplace_back(eigen_point, radius, surfaceIndex);
}

void UStringGenerator::endBrushPath()
{

}

void UStringGenerator::_initPath()
{
	_t = 0.0f;
	_brushPoints.clear();
	_lastStringElement = FGraphNodeHandle::null;
}


