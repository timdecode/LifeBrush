// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"


#include "SwarmGenerator.h"
#include "Simulation/MeshFilamentSimulation.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "Simulation/FlexElements.h"

void USwarmGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::attach(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;
	_graph = &context->graph();

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();
	meshSim->attach();

	USwarmSimulation * swarmSim = _simulationManager->registerSimulation<USwarmSimulation>();
	swarmSim->minNC1Separation = minNC1Separation;
	swarmSim->radius_segment = radius_segment;
	swarmSim->scaleFactor_segment = scaleFactor_segment;
	swarmSim->material_segment = material_segment;
	swarmSim->staticMesh_segment = staticMesh_segment;

	swarmSim->radius_7S = radius_7S;
	swarmSim->scaleFactor_7S = scaleFactor_7S;
	swarmSim->material_7S = material_7S;
	swarmSim->staticMesh_7S = staticMesh_7S;

	UMeshFilamentSimulation * filamentSim = _simulationManager->registerSimulation<UMeshFilamentSimulation>();
	filamentSim->attach();

	rand.GenerateNewSeed();

	_reloadBVH();
}

void USwarmGenerator::detach()
{
	Super::detach();

	UMeshSimulation * meshSim = _simulationManager->registerSimulation<UMeshSimulation>();

	meshSim->detach();
}

void USwarmGenerator::flexTick(
	float deltaT, 
	NvFlexVector<int>& neighbourIndices, 
	NvFlexVector<int>& neighbourCounts, 
	NvFlexVector<int>& apiToInternal, 
	NvFlexVector<int>& internalToAPI, 
	int maxParticles)
{
}

FQuat USwarmGenerator::_randomQuat()
{
	FVector axis = rand.GetUnitVector();
	float angle = rand.GetFraction() * M_PI * 2.0f;

	return FQuat(axis, angle);
}

void USwarmGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{

}

void USwarmGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	bool hit = false;

	decltype(_spaceBVH)::AABB_t query(&point.X, minNC1Separation * .3f);

	_spaceBVH.query(query, [&hit](unsigned int particleIndex) {
		hit = true;

		return false; // don't continue
	});

	if (hit)
		return;

	flexSimulation->addTickWork([&, point]() {
		_reloadBVH();

		_graph->beginTransaction();

		if (brushType == ESwarmGenerator_BrushType::Star)
			_createStar(point, _randomQuat());
		else if (brushType == ESwarmGenerator_BrushType::Anchor)
			_createAnchor(point, FQuat::Identity);

		_graph->endTransaction();
	});
}

void USwarmGenerator::endBrushPath()
{

}

FGraphNodeHandle USwarmGenerator::_createAnchor(FVector position, FQuat orientation)
{
	FGraphNodeHandle anchorHandle = _graph->node(_graph->addNode(position, orientation, radius_NC1 * scaleFactor_NC1));
	{
		// the star node is only guaranteed to point to valid memory in this block, where we haven't
		// added any nodes to the graph
		FGraphNode& anchorNode = anchorHandle(*_graph);

		_spaceBVH.insertParticle(anchorNode.id, &anchorNode.position.X, radius_NC1);

		FGraphMesh& mesh = anchorNode.addComponent<FGraphMesh>(*_graph);

		mesh.material = material_NC1;
		mesh.staticMesh = staticMesh_NC1;

		FBoid& boid = anchorNode.addComponent<FBoid>(*_graph);
		boid.radius = radius_NC1;

		anchorNode.addComponent<FStaticPositionObject>(*_graph).position = position;
		anchorNode.addComponent<FFlexParticleObject>(*_graph);
		anchorNode.addComponent<FBoidStar>(*_graph);
		anchorNode.addComponent<FBoidSeeker>(*_graph);
	}

	return anchorHandle;
}

FGraphNodeHandle USwarmGenerator::_createStar(FVector position, FQuat orientation)
{
	FGraphNodeHandle starHandle = _graph->node(_graph->addNode(position, orientation, radius_NC1 * scaleFactor_NC1));
	{
		// the star node is only guaranteed to point to valid memory in this block, where we haven't
		// added any nodes to the graph
		FGraphNode& starNode = starHandle(*_graph);

		_spaceBVH.insertParticle(starNode.id, &starNode.position.X, radius_NC1);

		FGraphMesh& mesh = starNode.addComponent<FGraphMesh>(*_graph);

		mesh.material = material_NC1;
		mesh.staticMesh = staticMesh_NC1;

		FBoid& boid = starNode.addComponent<FBoid>(*_graph);
		boid.radius = radius_NC1;

		starNode.addComponent<FFlexParticleObject>(*_graph);
		starNode.addComponent<FVelocityGraphObject>(*_graph);
		starNode.addComponent<FBoidStar>(*_graph);
	}

	const std::array<FQuat, 4> orientations = {
		FQuat(FVector::ForwardVector, M_PI * (0.0f / 2.0f)),
		FQuat(FVector::ForwardVector, M_PI * (1.0f / 2.0f)),
		FQuat(FVector::ForwardVector, M_PI * (2.0f / 2.0f)),
		FQuat(FVector::ForwardVector, M_PI * (3.0f / 2.0f))
	};

	const float offset = radius_segment + radius_NC1;

	for (FQuat localRotation : orientations)
	{
		FQuat newOrientation = orientation * localRotation;

		FVector dir = newOrientation * FVector::RightVector;

		FGraphNodeHandle wanderer = _createWanderer(position + dir * offset, newOrientation);

		FBoidGenerator& generator = wanderer(_graph).component<FBoidGenerator>(*_graph);
		{
			generator.last = starHandle;
			generator.starNode = starHandle;
		}


		float length = FVector::Dist(wanderer(_graph).position, position);
		_graph->connectNodes<FFlexConnection>(starHandle, wanderer, length, 0.9f);
	}

	return starHandle;
}

FGraphNodeHandle USwarmGenerator::_createWanderer(FVector position, FQuat orientation)
{
	FGraphNode& wandererNode = _graph->node(_graph->addNode(position, orientation, radius_segment * scaleFactor_segment));

	_spaceBVH.insertParticle(wandererNode.id, &wandererNode.position.X, radius_segment);

	FGraphMesh& mesh = wandererNode.addComponent<FGraphMesh>(*_graph);

	mesh.material = material_segment;
	mesh.staticMesh = staticMesh_segment;

	FBoid& boid = wandererNode.addComponent<FBoid>(*_graph);
	boid.radius = radius_segment;

	FFlexParticleObject& particle = wandererNode.addComponent<FFlexParticleObject>(*_graph);
	particle.group = 1;

	wandererNode.addComponent<FVelocityGraphObject>(*_graph);

	FBoidGenerator& boidGenerator = wandererNode.addComponent<FBoidGenerator>(*_graph);
	{
		boidGenerator.maxCount = 10;
		boidGenerator.filamentGroup = _simulationManager->simulation<UMeshFilamentSimulation>()->nextGroup();
	}


	FRandomWalkGraphObject& walker = wandererNode.addComponent<FRandomWalkGraphObject>(*_graph);

	return wandererNode.handle();
}

void USwarmGenerator::_reloadBVH()
{
	_spaceBVH.removeAll();

	auto& storage = _graph->componentStorage<FBoid>();

	for (FBoid& boid : storage)
	{
		if (!boid.isValid())
			continue;

		FGraphNode& node = _graph->node(boid.nodeHandle());

		_spaceBVH.insertParticle(boid.nodeIndex, &node.position.X, boid.radius);
	}
}


















void USwarmSimulation::begin()
{
	if (debugDrawBoids)
		_attachSegmentMeshes();
	else
		_detachSegmentMeshes();
}

void USwarmSimulation::detach()
{
	_detachSegmentMeshes();
}

void USwarmSimulation::_attachSegmentMeshes()
{
	auto& segmentStorage = graph->componentStorage<FBoidSegment>();
	auto& meshStorage = graph->componentStorage<FGraphMesh>();

	for (FBoidSegment& segment : segmentStorage)
	{
		FGraphMesh * mesh = meshStorage.componentPtrForNode(segment.nodeHandle());

		if( mesh )
			continue;

		mesh = &segment.nodeHandle()(graph).addComponent<FGraphMesh>(*graph);

		mesh->material = material_segment;
		mesh->staticMesh = staticMesh_segment;
	}
}

void USwarmSimulation::_detachSegmentMeshes()
{
	auto& segmentStorage = graph->componentStorage<FBoidSegment>();
	auto& meshStorage = graph->componentStorage<FGraphMesh>();

	for (FBoidSegment& segment : segmentStorage)
	{
		FGraphNode& node = segment.node(*graph);

		if( !node.hasComponent<FGraphMesh>() )
			continue;

		node.removeComponent<FGraphMesh>(*graph);
	}
}

void USwarmSimulation::tick(float deltaT)
{

}

void USwarmSimulation::flexTick(float deltaT, NvFlexVector<int>& neighbourIndices, NvFlexVector<int>& neighbourCounts, NvFlexVector<int>& apiToInternal, NvFlexVector<int>& internalToAPI, int maxParticles)
{
	graph->beginTransaction();

	_tickBoidStars(deltaT);

	_tickBoidGenerators();

	_tickBoidSeekers(deltaT);
	//_tickBoidSeekersFlexNeighbours(maxParticles, apiToInternal, neighbourCounts, internalToAPI, neighbourIndices);


	graph->endTransaction();
}

void USwarmSimulation::_tickBoidStars(float deltaT)
{
	_loadStarBVH();

	auto& storage = graph->componentStorage< FBoidStar>();

	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	for (FBoidStar& star : storage)
	{
		if( !star.isValid() )
			continue;

		FGraphNode& node = graph->node(star.nodeHandle());

		FVelocityGraphObject * velocity = velocities.componentPtrForNode(star.nodeHandle());

		if( !velocity )
			continue;

		decltype(_starBVH)::AABB_t query(&node.position.X, minNC1Separation);

		_starBVH.query(query, [&](unsigned int particleIndex) {
			FGraphNodeHandle handle(particleIndex);

			if (handle == star.nodeHandle())
				return true;

			FGraphNode& particleNode = graph->node(handle);

			FVector dir = particleNode.position - node.position;
			float length = dir.Size();
			dir = dir.GetSafeNormal();

			if (length < minNC1Separation)
			{
				// approximate a spring
				float hookLength = minNC1Separation - length;
			
				const float c = 10.0f;

				velocity->linearVelocity += -dir * hookLength * c * deltaT;
			}

			return true; // continue
		});
	}
}


void USwarmSimulation::_tickBoidSeekersFlexNeighbours(int maxParticles, NvFlexVector<int>& apiToInternal, NvFlexVector<int>& neighbourCounts, NvFlexVector<int>& internalToAPI, NvFlexVector<int>& neighbourIndices)
{
	const int stride = maxParticles;

	float bindingDistSqrd = std::pow(radius_segment * 2.0f, 2.0f);

	auto& seekerStorage = graph->componentStorage<FBoidSeeker>();
	auto& velocityStorage = graph->componentStorage<FVelocityGraphObject>();

	for (FBoidSeeker& seeker : seekerStorage)
	{
		if (!seeker.isValid())
			continue;

		FGraphNode& seekerNode = graph->node(seeker.nodeHandle());

		int flexIndex = apiToInternal[seeker.nodeIndex];
		int neighborCount = neighbourCounts[flexIndex];

		// find the nearest
		float minDistanceSqrd = std::numeric_limits<float>::max();
		FGraphNodeHandle minNodeHandle = FGraphNodeHandle::null;
		{
			FVector p = seekerNode.position;

			for (int ni = 0; ni < neighborCount; ++ni)
			{
				FGraphNodeHandle neighborNodeHandle = FGraphNodeHandle(internalToAPI[neighbourIndices[ni*stride + flexIndex]]);

				if (neighborNodeHandle == seeker.nodeHandle())
					continue;

				FBoidSeeker * neighbourSeeker = seekerStorage.componentPtrForNode(neighborNodeHandle);

				if (!neighbourSeeker || !neighbourSeeker->isValid())
					continue;

				// don't connect seekers in the same star
				if( neighbourSeeker->starNode == seeker.starNode )
					continue;

				FGraphNode& node = graph->node(neighborNodeHandle);

				float distSqrd = (node.position - p).SizeSquared();

				if (distSqrd < minDistanceSqrd)
				{
					minDistanceSqrd = distSqrd;
					minNodeHandle = neighborNodeHandle;
				}
			}
		}

		if (minNodeHandle)
		{
			FGraphNode& neighbourNode = graph->node(minNodeHandle);

			FBoidSeeker * neighbourSeeker = seekerStorage.componentPtrForNode(minNodeHandle);

			// attract to the neighbor
			FVelocityGraphObject * velocity = velocityStorage.componentPtrForNode(seeker.nodeHandle());

			if (!velocity)
				continue;

			FVector dir = (neighbourNode.position - seekerNode.position);
			float distSqrd = dir.SizeSquared();
			dir = dir.GetSafeNormal();

			// set velocity (not n-body for now)
			velocity->linearVelocity = dir * 5.0f;

			// bind
			if (distSqrd < bindingDistSqrd)
			{
				graph->connectNodes<FFlexConnection>(seekerNode.handle(), neighbourNode.handle(), radius_segment * 2.0f, 0.9f);

				seekerNode.removeComponent<FBoidSeeker>(*graph);
				neighbourNode.removeComponent<FBoidSeeker>(*graph);

				if( seekerNode.hasComponent<FRandomWalkGraphObject>() )
					seekerNode.removeComponent<FRandomWalkGraphObject>(*graph);
				if (neighbourNode.hasComponent<FRandomWalkGraphObject>())
					neighbourNode.removeComponent<FRandomWalkGraphObject>(*graph);

				FGraphMesh& seekerMesh = seekerNode.component<FGraphMesh>(*graph);
				seekerMesh.material = material_segment;
				seekerMesh.markDirty();

				FGraphMesh& neighbourMesh = neighbourNode.component<FGraphMesh>(*graph);
				neighbourMesh.material = material_segment;
				neighbourMesh.markDirty();
			}
		}
	}
}

void USwarmSimulation::_tickBoidGenerators()
{
	auto& generatorStorage = graph->componentStorage<FBoidGenerator>();

	float minSeparationSqrd = std::pow(radius_segment * 2.0f, 2.0f);

	for (FBoidGenerator& boid : generatorStorage)
	{
		if (!boid.isValid() || !boid.last)
			continue;

		FGraphNode boidNode = graph->node(boid.nodeHandle());

		FGraphNode& lastNode = graph->node(boid.last);

		FVector dir = (lastNode.position - boidNode.position);
		float distanceSqrd = dir.SizeSquared();

		// are we far enough from the last guy to create a new segment?
		if (distanceSqrd < minSeparationSqrd)
			continue;

		dir = dir.GetSafeNormal();

		FVector segmentPosition = lastNode.position + dir * radius_segment * 2.0f;

		FGraphNodeHandle segment;

		// is this the last segment in the chain?
		if (boid.counter >= boid.maxCount - 1)
		{
			graph->removeNode(boidNode.id);

			segment = _createSeeker(segmentPosition, boidNode.orientation, boid.starNode);
		}
		else
		{
			// create a segment
			segment = _createSegment(segmentPosition, boidNode.orientation);
		}

		auto oneRingEdge = graph->connectNodes<FFlexConnection>(boid.last, segment, radius_segment * 2.0f, 0.5f);

		auto& filament = graph->addOrReplaceEdgeObject<FFilamentConnection>(oneRingEdge);
		{
			filament.radius = 0.5f;
			filament.group = boid.filamentGroup;
			filament.segmentID = boid.counter;
		}

		if (boid.last2ring)
			graph->connectNodes<FFlexConnection>(boid.last2ring, segment, radius_segment * 4.0f, 0.5f);

		boid.last2ring = boid.last;
		boid.last = segment;
		boid.counter++;

		// don't touch boidNode anymore, it could be invalid
	}
}

void USwarmSimulation::_tickBoidSeekers(float deltaT)
{
	BVH_t seekerBVH;

	const float bindingDistSqrd = std::pow(radius_segment * 2.0f, 2.0f);


	auto& seekerStorage = graph->componentStorage<FBoidSeeker>();
	auto& velocityStorage = graph->componentStorage<FVelocityGraphObject>();

	for (FBoidSeeker& seeker : seekerStorage)
	{
		if (!seeker.isValid())
			continue;

		FGraphNode& node = graph->node(seeker.nodeHandle());

		if (!node.isValid())
			continue;

		seekerBVH.insertParticle(node.id, &node.position.X, radius_segment);
	}

	for (FBoidSeeker& seeker : seekerStorage)
	{
		FGraphNode& seekerNode = graph->node(seeker.nodeHandle());

		decltype(_starBVH)::AABB_t query(&seekerNode.position.X, minNC1Separation);

		float minDistanceSqrd = std::numeric_limits<float>::max();
		FGraphNodeHandle minNodeHandle = FGraphNodeHandle::null;

		seekerBVH.query(query, [&](unsigned int particleIndex) {
			FGraphNodeHandle handle(particleIndex);

			if (handle == seeker.nodeHandle())
				return true;

			FGraphNode& particleNode = graph->node(handle);

			FBoidSeeker * particleSeeker = seekerStorage.componentPtrForNode(handle);

			if (particleSeeker->starNode == seeker.starNode)
				return true;

			float distSqrd = (particleNode.position - seekerNode.position).SizeSquared();

			if (distSqrd < minDistanceSqrd)
			{
				minDistanceSqrd = distSqrd;
				minNodeHandle = handle;
			}

			return true; // continue
		});

		if (minNodeHandle)
		{
			FGraphNode& neighbourNode = graph->node(minNodeHandle);

			FBoidSeeker * neighbourSeeker = seekerStorage.componentPtrForNode(minNodeHandle);

			// attract to the neighbor
			FVelocityGraphObject * velocity = velocityStorage.componentPtrForNode(seeker.nodeHandle());

			if (!velocity)
				continue;

			FVector dir = (neighbourNode.position - seekerNode.position);
			float distSqrd = dir.SizeSquared();
			dir = dir.GetSafeNormal();

			// set velocity (not n-body for now)
			velocity->linearVelocity = dir * 5.0f;

			// bind
			if (distSqrd < bindingDistSqrd)
			{
				graph->connectNodes<FFlexConnection>(seekerNode.handle(), neighbourNode.handle(), radius_segment * 2.0f, 0.9f);

				seekerNode.removeComponent<FBoidSeeker>(*graph);
				neighbourNode.removeComponent<FBoidSeeker>(*graph);

				if (seekerNode.hasComponent<FRandomWalkGraphObject>())
					seekerNode.removeComponent<FRandomWalkGraphObject>(*graph);
				if (neighbourNode.hasComponent<FRandomWalkGraphObject>())
					neighbourNode.removeComponent<FRandomWalkGraphObject>(*graph);

				FGraphMesh& seekerMesh = seekerNode.component<FGraphMesh>(*graph);
				seekerMesh.material = material_segment;
				seekerMesh.markDirty();

				FGraphMesh& neighbourMesh = neighbourNode.component<FGraphMesh>(*graph);
				neighbourMesh.material = material_segment;
				neighbourMesh.markDirty();
			}
		}
	}
}

FGraphNodeHandle USwarmSimulation::_createSegment(FVector position, FQuat orientation)
{
	FGraphNode& segmentNode = graph->node(graph->addNode(position, orientation, radius_segment * scaleFactor_segment));

	if (debugDrawBoids)
	{
		FGraphMesh& mesh = segmentNode.addComponent<FGraphMesh>(*graph);

		mesh.material = material_segment;
		mesh.staticMesh = staticMesh_segment;
		mesh.visible = debugDrawBoids;
	}

	FBoid& boid = segmentNode.addComponent<FBoid>(*graph);
	boid.radius = radius_segment;

	segmentNode.addComponent<FFlexParticleObject>(*graph);
	segmentNode.addComponent<FVelocityGraphObject>(*graph);
	segmentNode.addComponent<FBoidSegment>(*graph);

	return segmentNode.handle();
}

FGraphNodeHandle USwarmSimulation::_createSeeker(FVector position, FQuat orientation, FGraphNodeHandle starNode)
{
	FGraphNode& segmentNode = graph->node(graph->addNode(position, orientation, radius_7S * scaleFactor_7S));

	FGraphMesh& mesh = segmentNode.addComponent<FGraphMesh>(*graph);

	mesh.material = material_7S;
	mesh.staticMesh = staticMesh_7S;

	FBoid& boid = segmentNode.addComponent<FBoid>(*graph);
	boid.radius = radius_7S;

	segmentNode.addComponent<FFlexParticleObject>(*graph);
	segmentNode.addComponent<FVelocityGraphObject>(*graph);
	segmentNode.addComponent<FRandomWalkGraphObject>(*graph);
	FBoidSeeker& seeker = segmentNode.addComponent<FBoidSeeker>(*graph);
	seeker.starNode = starNode;

	return segmentNode.handle();
}

void USwarmSimulation::_loadStarBVH()
{
	_starBVH.removeAll();

	auto& storage = graph->componentStorage<FBoidStar>();

	for (FBoidStar& star : storage)
	{
		if (!star.isValid())
			continue;

		FGraphNode& node = graph->node(star.nodeHandle() );

		if( !node.isValid() )
			continue;

		_starBVH.insertParticle(star.nodeIndex, &node.position.X, radius_segment);
	}
}

