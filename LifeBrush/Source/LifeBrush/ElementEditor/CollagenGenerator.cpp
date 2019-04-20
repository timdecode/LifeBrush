//
//  Created by Timothy Davison on 2018-12-29.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ShipEditorSimulation/MeshSimulation.h"

#include "Simulation/FlexElements.h"

#include "CollagenGenerator.h"

void UCollagenGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::attach(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	meshSim->attach();

	rand.GenerateNewSeed();

	_loadFrontier();
}

void UCollagenGenerator::detach()
{
	Super::detach();

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	meshSim->detach();
}

void UCollagenGenerator::tick(float deltaT)
{
	if (!_context)
		return;

	_context->graph().pushTransactionContext();
	_spaceColonization();
	_context->graph().popTransactionContext();
}



void UCollagenGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex)
{
	_brushPoints.pts.clear();
	_index = std::make_unique< PositionFacePairCloudAdaptor >(3, _brushPoints, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));

	_brushBVH.removeAll();

	start = point;

	if (!_context)
		return;

	flexSimulation->addTickWork([&, point]() {
		_context->graph().beginTransaction();

		_conditionallySeedNC1Head(point);

		_context->graph().endTransaction();
	});
}

void UCollagenGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	Eigen::Vector3f eigenPoint = eigen(point);

	// short-circuit if we are near a previous point?
	if (!_brushPoints.pts.empty())
	{
		auto last = _brushPoints.pts.back().position;

		float reject = std::pow(radius * .25f, 2.0f);

		if ((eigenPoint - last).squaredNorm() < reject)
			return;
	}

	size_t i = _brushPoints.pts.size();
	_brushPoints.pts.emplace_back(eigenPoint, radius, surfaceIndex);

	_index->addPoints(i, i);
	_brushBVH.insertParticle(i, &eigenPoint.x(), radius);

	flexSimulation->addTickWork([&]() {
		_context->graph().beginTransaction();
		_conditionallySeedNC1Head(point);
		_context->graph().endTransaction();
	});
}

void UCollagenGenerator::endBrushPath()
{
	_brushPoints.pts.clear();
	_index = nullptr;
	_brushBVH.removeAll();
}

FQuat UCollagenGenerator::_randomQuat()
{
	FVector axis = rand.GetUnitVector();
	float angle = rand.GetFraction() * M_PI * 2.0f;

	return FQuat(axis, angle);
}

FGraphNode& UCollagenGenerator::_createNC1Head(FVector position, FQuat orientation)
{
	FGraph& graph = _context->graph();

	float scale = radius_segment * scaleFactor_segment;

	FGraphNodeHandle handle = _context->domain.insert(position, orientation, scale);
	FGraphNode& node = handle(graph);

	FGraphMesh& mesh = node.addComponent<FGraphMesh>(graph);
	mesh.material = material_NC1;
	mesh.staticMesh = staticMesh_NC1;

	FElementObject& element = node.component<FElementObject>(graph);
	element.radius = radius_NC1;

	FCollagen_NC1Head& nc1 = node.addComponent<FCollagen_NC1Head>(graph);

	FFlexParticleObject& particle = node.addComponent<FFlexParticleObject>(graph);

	FVelocityGraphObject& velocity = node.addComponent<FVelocityGraphObject>(graph);
	
	_spaceBVH.insertParticle(node.id, &node.position.X, element.radius);

	return node;
}

FGraphNode& UCollagenGenerator::_create7STail(FVector position, FQuat orientation)
{
	FGraph& graph = _context->graph();

	float scale = radius_segment * scaleFactor_segment;

	FGraphNodeHandle handle = _context->domain.insert(position, orientation, scale);
	FGraphNode& node = handle(graph);

	FGraphMesh& mesh = node.addComponent<FGraphMesh>(graph);
	mesh.material = material_7S;
	mesh.staticMesh = staticMesh_7S;

	FElementObject& element = node.component<FElementObject>(graph);
	element.radius = radius_7S;

	FCollagen_7STail& _7s = node.addComponent<FCollagen_7STail>(graph);

	FFlexParticleObject& particle = node.addComponent<FFlexParticleObject>(graph);

	FVelocityGraphObject& velocity = node.addComponent<FVelocityGraphObject>(graph);

	_spaceBVH.insertParticle(node.id, &node.position.X, element.radius);

	return node;
}

// Segments don't go in the bvh
FGraphNode& UCollagenGenerator::_createSegment(FVector position)
{
	FGraph& graph = _context->graph();

	float scale = radius_segment * scaleFactor_segment;

	FGraphNodeHandle handle = _context->domain.insert(position, FQuat::Identity, scale);
	FGraphNode& node = handle(graph);

	FCollagen_Segment& segment = node.addComponent<FCollagen_Segment>(graph);

	FGraphMesh& mesh = node.addComponent<FGraphMesh>(graph);
	mesh.material = material_segment;
	mesh.staticMesh = staticMesh_segment;

	FElementObject& element = node.component<FElementObject>(graph);
	element.radius = radius_segment;

	FFlexParticleObject& particle = node.addComponent<FFlexParticleObject>(graph);

	FVelocityGraphObject& velocity = node.addComponent<FVelocityGraphObject>(graph);

	return node;
}

FGraphNodeHandle UCollagenGenerator::_conditionallySeedNC1Head(FVector position)
{
	if (_hit(position, chainLength))
		return FGraphNodeHandle::null;

	// hack in a starting node
	FQuat orientation = _randomQuat();
	FGraphNode& nc1 = _createNC1Head(position, orientation);

	_frontier.insert(nc1.handle());

	return nc1.handle();
}

void UCollagenGenerator::_loadFrontier()
{
	if (!_context)
		return;

	_context->domain.rebalance();

	_spaceBVH.removeAll();

	FGraph& graph = _context->graph();
	for (auto& node : graph.allNodes)
	{
		if( !node.isValid() )
			continue;

		if( !node.hasComponent<FElementObject>() )
			continue;

		FElementObject& element = graph.component<FElementObject>(node.handle());

		_spaceBVH.insertParticle(node.id, &node.position.X, element.radius);
	}
}

void UCollagenGenerator::_link(FGraphNodeHandle start, FGraphNodeHandle end)
{
	FGraph& graph = _context->graph();
	Domain& domain = _context->domain;

	FGraphNode& startNode = graph.node(start);
	FGraphNode& endNode = graph.node(end);

	FVector dir = (endNode.position - startNode.position);

	const float length = dir.Size();

	dir = dir.GetSafeNormal();

	if (length < radius_segment)
	{
		return;
	}

	FGraphNodeHandle last = start;
	for (float t = radius_segment; t < chainLength && t < length; t += 2.0f * radius_segment)
	{
		FVector p = startNode.position + dir * t;

		FGraphNodeHandle segment = _createSegment(p);

		graph.connectNodes(last, segment);

		last = segment;
	}

	if (last != start)
	{
		graph.connectNodes(last, end);
	}
}

FQuat UCollagenGenerator::_rotationAroundY(float rad)
{
	return FQuat(FVector::RightVector, rad);
}

bool UCollagenGenerator::_inBounds(FVector p)
{
	if (_brushPoints.pts.empty())
		return false;

	const float r = 0.1f;

	BVH_t::AABB_t aabb(&p.X, r);

	auto overlaps = _brushBVH.query(aabb);

	return !overlaps.empty();
}


bool UCollagenGenerator::_hit(FVector occupied, float radius)
{
	decltype(_spaceBVH)::AABB_t query(&occupied.X, radius);

	bool hasOverlap = false;
	_spaceBVH.query(query, [&hasOverlap](unsigned int particleIndex) {
		hasOverlap = true;

		return false; // short-circuit on the first overlap
	});

	return hasOverlap;
}

void UCollagenGenerator::_spaceColonization()
{
	if (!_context)
		return;

	// Space colonization: 
	// - add points around horizon
	// - remove points from horizon that are not satisfiable or have been satisfied

	FGraph& graph = _context->graph();

	// work on a local version of the frontier, so, add/remove operations will affect the next round
	std::vector<FGraphNodeHandle> local;
	local.reserve(_frontier.size());

	for (FGraphNodeHandle handle : _frontier)
	{
		local.push_back(handle);
	}

	for (FGraphNodeHandle handle : local)
	{
		FGraphNode node = handle(graph);

		if( !_inBounds(node.position) )
			continue;

		if (_isSatisfied(handle))
		{
			_frontier.erase(handle);
			continue;
		}

		if (!_isSatisfiable(handle))
		{
			_frontier.erase(handle);
			continue;
		}


		if (node.hasComponent<FCollagen_NC1Head>())
			_process_nc1Head(node);

		if (node.hasComponent<FCollagen_7STail>())
			_process_7sTail(node);
	}
}

bool UCollagenGenerator::_isSatisfied(FGraphNodeHandle handle)
{
	FGraph& graph = _context->graph();

	FGraphNode& node = handle(graph);

	if (node.hasComponent<FCollagen_7STail>())
	{
		FCollagen_7STail& _7s = graph.component<FCollagen_7STail>(node.handle());

		return _7s.bound(graph);
	}

	if (node.hasComponent<FCollagen_NC1Head>())
	{
		FCollagen_NC1Head& _nc1 = graph.component<FCollagen_NC1Head>(node.handle());

		return _chainsOut(node, _nc1).size() == 4;
	}

	return false;
}

bool UCollagenGenerator::_isSatisfiable(FGraphNodeHandle handle)
{
	FGraph& graph = _context->graph();

	FGraphNode& node = handle(graph);

	if (node.hasComponent<FCollagen_NC1Head>())
	{
		FCollagen_NC1Head& _nc1 = graph.component<FCollagen_NC1Head>(node.handle());

		return _nc1.satisfiable;
	}

	if (node.hasComponent<FCollagen_7STail>())
	{
		FCollagen_7STail& _7s = graph.component<FCollagen_7STail>(node.handle());

		auto tails = _7s.collagenTails(graph);

		if (!tails.first)
			return false;

		FGraphNode& first = graph.node(tails.first);

		FVector dir = (first.position - node.position).GetUnsafeNormal();

		FVector p = -dir * radius_segment + node.position;

		auto nearest = _context->domain.nearest(eigen(p));

		if (nearest.distanceSquared < std::pow(radius_segment*.9f, 2.0f) && nearest.element != handle)
			return false;
	}

	return true;
}

void UCollagenGenerator::_process_nc1Head(FGraphNode node)
{
	FGraph& graph = _context->graph();
	Domain& domain = _context->domain;
	
	FCollagen_NC1Head& nc1 = graph.component<FCollagen_NC1Head>(node.handle());

	if (!nc1.satisfiable)
		return;

	auto chains = _chainsOut(node, nc1);

	FVector dir;

	int start = 0;

	if (chains.size() == 0)
	{
		dir = node.orientation.RotateVector(FVector::ForwardVector);
		start = 0;
	}
	else if (chains.size() > 1)
	{
		nc1.satisfiable = false;
		return;
	}
	else
	{
		FGraphNode& n0 = graph.node(chains[0]);

		dir = (n0.position - node.position).GetSafeNormal();

		start = 1;
	}


	const float halfChain = 0.5f * chainLength;

	const unsigned int maxAttempts = 10;

	for (int i = start; i < 4; ++i)
	{
		int attempt = 0;
		for (attempt = 0; attempt < maxAttempts; ++attempt)
		{
			const float alpha_lower = (2.0f * float(i - 1) * M_PI) / 4.0f;
			const float alpha_upper = (2.0f * float(i + 1) * M_PI) / 4.0f;

			const float alpha = rand.FRandRange(alpha_lower, alpha_upper);
			const FQuat rotation = node.orientation * _rotationAroundY(alpha);
			const FVector localDir = rotation * dir;

			FVector p_check = node.position + localDir * chainLength;

			auto neighbours = domain.nearestInRadius(eigen(p_check), halfChain);

			FGraphNodeHandle foundTail = FGraphNodeHandle::null;

			// Rule 1. Connect to neighbours
			for (auto h : neighbours)
			{
				FGraphNode& neighbour = h(graph);

				if (!neighbour.hasComponent<FCollagen_7STail>())
					continue;

				auto found = std::find(chains.begin(), chains.end(), h);

				if( found != chains.end() )
					continue;

				FCollagen_7STail& tail = neighbour.component<FCollagen_7STail>(graph);

				if (!tail.bound(graph))
				{
					_frontier.erase(tail.nodeHandle());

					graph.connectNodes(node.handle(), tail.nodeHandle());

					_link(node.handle(), tail.nodeHandle());

					foundTail = tail.nodeHandle();

					break;
				}
			}

			// Rule 2. Create a new tail
			if (foundTail)
				break;
			else
			{
				if( _hit(p_check, chainLength * .2f) )
					continue;

				FVector newDir = rand.GetUnitVector();

				if (FVector::DotProduct(newDir, localDir) > 0.0f)
					newDir *= -1.0f;

				float newAngle = rand.RandRange(0.0f, 2.0f * M_PI);
				FQuat newOrientation(newDir, newAngle);

				FGraphNodeHandle newTail = _create7STail(p_check, newOrientation);

				graph.connectNodes(node.handle(), newTail);

				_link(node.handle(), newTail);

				_frontier.insert(newTail);

				chains.push_back(newTail);

				break;
			}
		}

		if (attempt == maxAttempts)
		{
			nc1.satisfiable = false;
			_frontier.erase(node.handle());
		}
	}

	_frontier.erase(node.handle());
}

void UCollagenGenerator::_process_7sTail(FGraphNode node)
{
	FGraph& graph = _context->graph();
	Domain& domain = _context->domain;

	FCollagen_7STail& tail = graph.component<FCollagen_7STail>(node.handle());

	if (tail.bound(graph))
		return;

	auto collagens = tail.collagenTails(graph);

	if (!collagens.first)
	{
		return;
	}

	FGraphNode& first = graph.node(collagens.first);
	
	FVector otherDir = (first.position - node.position).GetSafeNormal();

	const int maxAttempts = 10;
	int attempt = 0;

	for (attempt = 0; attempt < maxAttempts; ++attempt)
	{
		FVector dir;
		do 
		{
			dir = rand.GetUnitVector();
		} while (FVector::DotProduct(dir, otherDir) >= 0.0f );

		FVector p = node.position + dir * chainLength;

		if (_hit(p, chainLength * .2f))
			continue;

		// create a NC1 node
		FGraphNodeHandle nc1 = _createNC1Head(p, _randomQuat());

		graph.connectNodes(nc1, node.handle());

		_frontier.insert(nc1);

		_link(node.handle(), nc1);

		break;
	}

	_frontier.erase(node.handle());
}

std::vector<FGraphNodeHandle> UCollagenGenerator::_chainsOut(FGraphNode& node, FCollagen_NC1Head& nc1)
{
	std::vector<FGraphNodeHandle> result;

	FGraph& graph = _context->graph();

	for (auto e : node.edges)
	{
		FGraphEdge& edge = graph.edge(FGraphEdgeHandle(e));

		auto other = edge.other(node.handle());

		if (other(graph).hasComponent< FCollagen_Segment>())
			result.emplace_back(other);
	}

	return result;
}

auto FCollagen_7STail::bound(FGraph& graph) -> bool
{
	FGraphNode& node = graph.node(nodeHandle());

	size_t n = 0;
	for (auto e : node.edges)
	{
		FGraphNodeHandle other = graph.edge(FGraphEdgeHandle(e)).other(nodeHandle());

		if (other(graph).hasComponent<FCollagen_Segment>())
			n++;
	}

	return n >= 2;
}

auto FCollagen_7STail::collagenTails(FGraph& graph) -> std::pair<FGraphNodeHandle, FGraphNodeHandle>
{
	FGraphNode& node = graph.node(nodeHandle());

	std::pair<FGraphNodeHandle, FGraphNodeHandle> pair;

	size_t n = 0;
	for (auto e : node.edges)
	{
		FGraphNodeHandle other = graph.edge(FGraphEdgeHandle(e)).other(nodeHandle());

		if (other(graph).hasComponent<FCollagen_Segment>())
		{
			if (n == 0)
				pair.first = other;
			else if (n == 1)
				pair.second = other;

			n++;
		}
	}

	return pair;
}
