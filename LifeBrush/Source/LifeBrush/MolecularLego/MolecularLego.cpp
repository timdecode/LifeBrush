// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include <algorithm>

#include "Simulation/FlexElements.h"
#include "Simulation/Aggregates.h"

#include "PhysicsUtilities.h"

#include "MolecularLego.h"

const FMLSpeciesID FMLSpeciesID::null(-1);

void UMLParticleSimulation::attach()
{
	const int nSpecies = 4;
	_attractionMatrix.resize(nSpecies, nSpecies);

	// negative: repulsion
	// positive: attraction
	_attractionMatrix <<  0.0f,  10.0f,  0.0f,  0.0f,
		                 10.0f,   0.0f,  0.0f,  0.0f,
		                  0.0f,   0.0f,  5.0f,  0.0f,
		                  0.0f,   0.0f,  0.0f,  5.0f;

	_bindingRadiusMatrix.resize(nSpecies, nSpecies);

	// negative: no binding
	// positive: the minimum radius to create a bond 
	_bindingRadiusMatrix << -1.0f,   5.0, -1.0f, -1.0f,
		                     5.0,  -1.0f, -1.0f, -1.0f,
		                    -1.0f,  -1.0f,  4.0f, -1.0f,
		                    -1.0f,  -1.0f, -1.0f,  4.0f;
}

void UMLParticleSimulation::detach()
{

}

void UMLParticleSimulation::tick(float deltaT)
{
	//recentlyBonded.Empty();

	//_updateParticleBVH();

	//tick_mlParticles(deltaT);

	//tick_bindOrientedParticles_sortedInteractions(deltaT);

	//tick_brownianMotion(deltaT);
}

TArray<FGraphNodeHandle> UMLParticleSimulation::_createStar(
	FVector position,
	FQuat orientation, 
	float scale, 
	float separation, 
	FGraphNodeHandle rigidBody, 
	UStaticMesh * optionalMesh/*= nullptr*/, 
	UMaterialInterface * optionalMaterial /*= nullptr*/)
{
	TArray<FGraphNodeHandle> star;

	FVector up = FVector::UpVector;

	for (int i = 0; i < 4; ++i)
	{
		FVector p = position + orientation.RotateVector(up) * separation;

		auto newNode = FGraphNodeHandle(graph->addNode(p, orientation, scale));

		auto& particle = newNode(*graph).addComponent<FFlexParticleObject>(*graph);
		particle.group = rigidBody.index;
		particle.selfCollide = false;

		newNode(*graph).addComponent<FVelocityGraphObject>(*graph);

		if (optionalMesh && optionalMaterial)
		{
			FGraphMesh& mesh = newNode(*graph).addComponent<FGraphMesh>(*graph);

			mesh.staticMesh = optionalMesh;
			mesh.material = optionalMaterial;
		}

		star.Add(newNode);

		up = FQuat(FVector::ForwardVector, FMath::DegreesToRadians(90.0f)).RotateVector(up);
	}

	return star;
}

void UMLParticleSimulation::tick_mlParticles(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& particles = graph->componentStorage<FMLParticle>();
	auto& aggregateIds = graph->componentStorage<FMLAggregateNO_id>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIds.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};

	graph->each_node_object<FMLParticle>([&](FGraphNode& node, FMLParticle& particle)
	{
		FVelocityGraphObject * velocityObject = velocities.componentPtrForNode(node.handle());

		if (!velocityObject)
			return;

		FGraphNodeHandle aggregateID = getAggregateID(node.handle());

		FVector& linearVelocity = velocityObject->linearVelocity;

		unrealAABB::AABB query(node.position, hackNeighbourhoodRadius);

		_particleBVH.query(query, [&](unsigned int particleIndex) {
			FMLParticle& otherParticle = particles[particleIndex];
			FGraphNode& otherNode = graph->node(otherParticle.nodeHandle());

			// don't interact with things in the same aggregate
			FGraphNodeHandle otherAggregateID = getAggregateID(otherParticle.nodeHandle());

			if (otherAggregateID == aggregateID)
				return true;

			FVector dir = otherNode.position - node.position;
			float r = dir.Size();

			// Binding radius
			// --------------

			float bindingRadius = _bindingRadiusMatrix(particle.speciesID, otherParticle.speciesID);

			// bind if we are close enough
			if (r < bindingRadius
				&& particle.numActiveBonds < hackMaxNumBonds
				&& otherParticle.numActiveBonds < hackMaxNumBonds)
			{
				// if they share the same aggregate, don't bind
				auto particleAggregate = FMLAggregateNO::getAggregate(*graph, particle.nodeHandle());
				auto otherAggregate = FMLAggregateNO::getAggregate(*graph, otherParticle.nodeHandle());

				// bind!
				if (particleAggregate != otherAggregate)
				{
					graph->connectNodes<FFlexConnection>(particle.nodeHandle(), otherParticle.nodeHandle(),
						0.5f * bindingRadius, // rest length
						0.5f // spring coefficient
						);

					particle.numActiveBonds++;
					otherParticle.numActiveBonds++;
				}
			}

			// Attraction
			// ----------

			float c = _attractionMatrix(particle.speciesID, otherParticle.speciesID);

			// normalize
			dir = r > 0.001 ? dir / r : FVector::ZeroVector;

			r = std::max(r, _hackRadius);

			if (r < hackNeighbourhoodRadius)
				linearVelocity += dir * (c / std::pow(r, 2.0f)) * deltaT;

			return true; // continue
		});
	});
}

bool UMLParticleSimulation::_facing(FGraphNode& a, FGraphNode& b)
{
	// check if they are pointing at each other
	const FVector normal = a.orientation.RotateVector(FVector::ForwardVector);
	const FVector otherNormal = b.orientation.RotateVector(FVector::ForwardVector);
	const FVector dir = b.position - a.position;

	const float normalDotProduct = FVector::DotProduct(normal, otherNormal);
	const float dirDotProduct = FVector::DotProduct(normal, dir);

	const bool facing = normalDotProduct < 0.0f && dirDotProduct > 0.0f;

	const float minAngle = FMath::DegreesToRadians(120.0f);

	return facing && FMath::Abs(FMath::Acos(-normalDotProduct)) <= minAngle;
}



void UMLParticleSimulation::tick_bindOrientedParticles(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();
	auto& aggregateIDs = graph->componentStorage<FMLAggregateNO_id>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIDs.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};



	// find and cache our helper particles
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		if (particle._didInit)
			return;

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;
		rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
			float distSqrd = (helperHandle(*graph).position - node.position).SizeSquared();

			helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
		});

		std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
			return a.second < b.second;
		});

		size_t n = particle.subs.size();

		int i;
		for ( i = 0; i < helperHandles.size() && i < n; ++i)
		{
			particle.subs[i] = helperHandles[i].first;
		}

		// null the rest
		for (; i < n; ++i)
			particle.subs[i] = FGraphNodeHandle::null;

		particle._didInit = true;
	});

	// update our interaction partners
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node, FMLOrientedParticle& particle)
	{
		if (particle.numActiveBonds >= hackMaxNumBonds)
			return;

		if (particle.interactionMultiplier <= 0.0f)
			return;

		

		//if (particle._interactionPartner)
		//{
		//	FMLOrientedParticle * parntner = orientedParticles.componentPtrForNode(particle._interactionPartner);
		//	float partnerScaledNeighbourhoodRadius = _hackNeighbourhoodRadius * parntner->radiusMultiplier;

		//	FGraphNode& partnerNode = graph->node(particle._interactionPartner);

		//	float scaledNeighbourhoodRadius = _hackNeighbourhoodRadius * particle.radiusMultiplier;

		//	float distSqrd = (partnerNode.position - node.position).SizeSquared();


		//	// it's still valid
		//	if ((distSqrd < std::pow(scaledNeighbourhoodRadius, 2.0f) || distSqrd < std::pow(partnerScaledNeighbourhoodRadius, 2.0f)) &&
		//		_facing(node, partnerNode))
		//		return;
		//}

		FGraphNodeHandle aggregateID = getAggregateID(node.handle());

		unrealAABB::AABB query(node.position, hackNeighbourhoodRadius);

		FGraphNodeHandle closest;
		float closestDistSqrd = std::numeric_limits<float>::max();

		float scaledNeighbourhoodRadius = hackNeighbourhoodRadius * particle.radiusMultiplier;

		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& otherParticle = orientedParticles.at(particleIndex);

			float partnerScaledNeighbourhoodRadius = hackNeighbourhoodRadius * otherParticle.radiusMultiplier;


			if (otherParticle.numActiveBonds >= hackMaxNumBonds)
				return true;

			if (otherParticle._bondedPartner)
				return true;

			const float c = _attractionMatrix(particle.speciesID, otherParticle.speciesID);

			if (c <= 0.0f)
				return true;

			FGraphNodeHandle otherAggregateID = getAggregateID(otherParticle.nodeHandle());

			// don't interact with yourself
			if (otherParticle.nodeHandle() == particle.nodeHandle() || aggregateID == otherAggregateID)
				return true;

			FGraphNode& otherNode = graph->node(otherParticle.nodeHandle());
			float distSqrd = (otherNode.position - node.position).SizeSquared();

			if (distSqrd < closestDistSqrd && _facing(node, otherNode) &&
				(distSqrd < std::pow(scaledNeighbourhoodRadius, 2.0f) || distSqrd < std::pow(partnerScaledNeighbourhoodRadius, 2.0f)))
			{
				closestDistSqrd = distSqrd;
				closest = otherParticle.nodeHandle();
			}

			return true;
		});

		if (closest)
		{
			FMLOrientedParticle * parntner = orientedParticles.componentPtrForNode(particle._bondedPartner);

			if (parntner)
				parntner->_bondedPartner = FGraphNodeHandle::null;

			parntner = nullptr;

			particle._bondedPartner = closest;
			orientedParticles.componentForNode(closest)._bondedPartner = particle.nodeHandle();
		}
		else if (particle._bondedPartner)
		{
			FMLOrientedParticle * parntner = orientedParticles.componentPtrForNode(particle._bondedPartner);

			if (parntner)
				parntner->_bondedPartner = FGraphNodeHandle::null;

			particle._bondedPartner = FGraphNodeHandle::null;
		}
	});

	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		// already bound
		if (particle.numActiveBonds >= hackMaxNumBonds || particle.interactionMultiplier <= 0.0f )
			return;

		if (!particle._bondedPartner)
			return;

		FGraphNodeHandle aggregateID = getAggregateID(node.handle());

		FVector position = node.position;

		float scaledNeighbourhoodRadius = hackNeighbourhoodRadius * particle.radiusMultiplier;

		unrealAABB::AABB query(node.position, scaledNeighbourhoodRadius);

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		// count the particles and find the center of mass
		int32 rigidParticleCount = 0;
		FVector com = FVector::ZeroVector;
		centerOfMass(rigidHandle, rigidParticleCount, com);
		if (rigidParticleCount == 0)
			return;

		FMLOrientedParticle& otherParticle = orientedParticles.componentForNode(particle._bondedPartner);

		if (otherParticle.interactionMultiplier <= 0.0f)
			return;

		FGraphNode& otherNode = graph->node(otherParticle.nodeHandle());

		const float c = _attractionMatrix(particle.speciesID, otherParticle.speciesID);

		// short circuits
		{
			// already bound
			if (otherParticle.numActiveBonds >= hackMaxNumBonds)
				return;

			// they should be different aggregates
			if (aggregateID == getAggregateID(otherNode.handle()))
				return;

			if (c <= 0.0f)
				return;


		}

		const float halfBindingOffset = bindingOffset / 2.0f;

		bool positionallyAligned = false;

		const FVector framePosition = node.position + node.orientation.RotateVector(FVector::ForwardVector * halfBindingOffset);
		const FVector targetFramePosition = otherNode.position + otherNode.orientation.RotateVector(FVector::ForwardVector * halfBindingOffset);
		const FQuat targetRotation = FRotationMatrix::MakeFromXZ(-otherNode.orientation.GetAxisX(), otherNode.orientation.GetAxisZ()).ToQuat();

		FQuat rotationDifference = targetRotation * node.orientation.Inverse();

		// calculate positional alignment forces
		if (false) {
			FVector dir = targetFramePosition - position;

			const float distSqrd = dir.SizeSquared();
			const float dist = dir.Size();

			if (dist < 0.1f)
				positionallyAligned = true;

			dir = dir.GetSafeNormal();

			const float acc = std::min(16.0f, c * (1 / dist));

			FVector deltaV = dir * acc * deltaT;

			rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
					velocity->linearVelocity += deltaV;
			});
		}

		// critically damped positional spring
		// see: http://mathproofs.blogspot.com/2013/07/critically-damped-spring-smoothing.html
		{
			const float w = 20.0f * particle.interactionMultiplier;

			const FVector oldVelocity = node.component<FVelocityGraphObject>(*graph).linearVelocity;

			const float distSqrd = (framePosition - targetFramePosition).SizeSquared();

			const float bindingThreshold = 1.5f;

			positionallyAligned = distSqrd < std::pow(bindingThreshold, 2.0f);

			const FVector newVelocity = (oldVelocity - std::pow(w, 2.0f) * deltaT * (framePosition - targetFramePosition)) / std::pow(1 + w * deltaT, 2.0f);

			// clamp the deltaV
			const float maxDeltaV = 20.0f * deltaT;
			const FVector deltaV = (newVelocity - oldVelocity).GetClampedToMaxSize(maxDeltaV);

			rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
					velocity->linearVelocity += deltaV;
			});
		}

		bool rotationallyAligned = false;

		// calculate rotational alignment forces
		if (true) {

			FVector torqueAxis;
			float torqueMagnitude;

			rotationDifference.ToAxisAndAngle(torqueAxis, torqueMagnitude);

			torqueMagnitude = std::sin(torqueMagnitude);

			rotationallyAligned = std::abs(torqueMagnitude) < 0.5f;

			FVector hackyTorque = torqueAxis * torqueMagnitude * 100.0f;
			FVector hackyTorque_i = hackyTorque / float(rigidParticleCount);

			rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
				{
					const FVector r = (nodeHandle(*graph).position - com);

					velocity->linearVelocity += FVector::CrossProduct(hackyTorque_i, r) * deltaT;
				}
			});
		}

		if (positionallyAligned && rotationallyAligned)
		{
			graph->beginTransaction();

			FMatrix aboutRotation = FRotationAboutPointMatrix::Make(rotationDifference, node.position);
			FMatrix translation = FTranslationMatrix::Make(targetFramePosition - node.position);

			FMatrix toOther = translation * aboutRotation;

			for (FGraphNodeHandle a : particle.subs)
			{
				if( !a ) continue;

				FVector p_a = toOther.TransformPosition(a(*graph).position);

				for (FGraphNodeHandle b : otherParticle.subs)
				{
					if( !b ) continue;

					FVector p_b = b(*graph).position;

					float dist = (p_a - p_b).Size();
					float strength = a == b ? 0.5f : 0.1f;

					graph->connectNodes<FFlexConnection>(a, b, dist, strength);
				}
			}

			graph->endTransaction();

			// we've mutated the graph, make sure we updated our particle by handle
			graph->componentPtr<FMLOrientedParticle>(particle.nodeHandle())->numActiveBonds++;
			graph->componentPtr<FMLOrientedParticle>(otherParticle.nodeHandle())->numActiveBonds++;

			recentlyBonded.Add(particle.nodeHandle());
			recentlyBonded.Add(otherParticle.nodeHandle());
		}

		return;
	});
}







FQuat UMLParticleSimulation::mirror(FQuat rotation)
{
	return FRotationMatrix::MakeFromXZ(-rotation.GetAxisX(), rotation.GetAxisZ()).ToQuat();
}

bool UMLParticleSimulation::_aligned(FGraphNode& node_a, FMLOrientedParticle& particle_a, FGraphNode& node_b, FMLOrientedParticle& particle_b)
{
	if (!particle_a.alignedToBind && !particle_b.alignedToBind)
		return true;

	FQuat q_a = node_a.orientation;
	FQuat q_b = mirror(node_b.orientation);

	float distance = q_a.AngularDistance(q_b);

	return distance < 0.5f;
}

void UMLParticleSimulation::tick_brownianMotion(float deltaT)
{
	FRandomStream rand;
	rand.GenerateNewSeed();



	auto& browns = graph->componentStorage<FMLBrownianMotion>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();

	for (FMLBrownianMotion& brown : browns)
	{
		brown.time -= deltaT;

		if( brown.time > 0.0f )
			continue;

		if (!brown._rigidBody)
		{
			brown._rigidBody = FFlexRigidBodyObject::getRigidBodyHandle(*graph, brown.nodeHandle());
		}

		brown.time = rand.FRandRange(minTime, maxTime);

		FVector dv = rand.GetUnitVector() * rand.FRandRange(minSpeed, maxSpeed);
		FVector dtau = rand.GetUnitVector() * rand.FRandRange(minTau, maxTau);

		if (!brown._rigidBody)
		{
			if (auto velocity = velocities.componentPtrForNode(brown.nodeHandle()))
			{
				velocity->linearVelocity = (velocity->linearVelocity * 0.6f + dv).GetClampedToMaxSize(15.0f);
			}
		}
		else
		{
			// count the particles and find the center of mass
			int32 rigidParticleCount = 0;
			FVector com = FVector::ZeroVector;
			centerOfMass(brown._rigidBody, rigidParticleCount, com);
			if (rigidParticleCount == 0)
				return;

			brown._rigidBody(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
				{
					velocity->linearVelocity = (velocity->linearVelocity * 0.6f + dv).GetClampedToMaxSize(15.0f);

					const FVector r = (nodeHandle(*graph).position - com);
					velocity->linearVelocity += FVector::CrossProduct(dtau, r) * deltaT;
				}
			});
		}


	}
}

void UMLParticleSimulation::tick_bindOrientedParticles2(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();
	auto& aggregateIDs = graph->componentStorage<FMLAggregateNO_id>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIDs.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};




	// find and cache our helper particles
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		if (particle._didInit)
			return;

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;
		rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
			float distSqrd = (helperHandle(*graph).position - node.position).SizeSquared();

			helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
		});

		std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
			return a.second < b.second;
		});

		size_t n = particle.subs.size();

		int i;
		for (i = 0; i < helperHandles.size() && i < n; ++i)
		{
			particle.subs[i] = helperHandles[i].first;
		}

		// null the rest
		for (; i < n; ++i)
			particle.subs[i] = FGraphNodeHandle::null;

		particle._didInit = true;
	});

	TSet< InteractionPair > interactions;

	// update our interaction partners
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node, FMLOrientedParticle& particle)
	{
		if (particle.numActiveBonds >= hackMaxNumBonds)
			return;

		if (particle.interactionMultiplier <= 0.0f)
			return;

		FGraphNodeHandle aggregateID = getAggregateID(node.handle());

		unrealAABB::AABB query(node.position, hackNeighbourhoodRadius);

		FGraphNodeHandle closest;
		float closestDistSqrd = std::numeric_limits<float>::max();

		float scaledNeighbourhoodRadius = hackNeighbourhoodRadius * particle.radiusMultiplier;

		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& otherParticle = orientedParticles.at(particleIndex);

			float partnerScaledNeighbourhoodRadius = hackNeighbourhoodRadius * otherParticle.radiusMultiplier;


			if (otherParticle.numActiveBonds >= hackMaxNumBonds)
				return true;

			if (otherParticle._bondedPartner)
				return true;

			const float c = _attractionMatrix(particle.speciesID, otherParticle.speciesID);

			if (c <= 0.0f)
				return true;

			FGraphNodeHandle otherAggregateID = getAggregateID(otherParticle.nodeHandle());

			// don't interact with yourself
			if (otherParticle.nodeHandle() == particle.nodeHandle() || aggregateID == otherAggregateID)
				return true;

			FGraphNode& otherNode = graph->node(otherParticle.nodeHandle());
			float distSqrd = (otherNode.position - node.position).SizeSquared();

			if (_facing(node, otherNode) &&
				(distSqrd < std::pow(scaledNeighbourhoodRadius, 2.0f) || distSqrd < std::pow(partnerScaledNeighbourhoodRadius, 2.0f)))
			{
				auto pair = InteractionPair{ particle.nodeHandle(), otherParticle.nodeHandle() };
				interactions.Emplace(pair);
			}

			return true;
		});
	});

	// Sort, so that we can visit some of the L1 and L2 cache in order
	interactions.Sort([&](const InteractionPair& a, const InteractionPair& b) {
		return a.a < b.a;
	});

	// Calculate our forces for each interaction pair
	for (auto& pair : interactions)
	{
		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		FMLOrientedParticle& particle_a = orientedParticles.componentForNode(pair.a);
		FMLOrientedParticle& particle_b = orientedParticles.componentForNode(pair.b);


		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());
		auto rigidHandle_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_b.handle());

		// count the particles and find the center of mass
		int32 rigidCount_a = 0;
		FVector com_a = FVector::ZeroVector;
		centerOfMass(rigidHandle_a, rigidCount_a, com_a);
		if (rigidCount_a == 0)
			return;

		int32 rigidCount_b = 0;
		FVector com_b = FVector::ZeroVector;
		centerOfMass(rigidHandle_b, rigidCount_b, com_b);
		if (rigidCount_b == 0)
			return;

		const float halfBindingOffset = bindingOffset / 2.0f;

		bool positionallyAligned = false;

		bool rotationallyAligned = false;

		{
			const FQuat q_a = node_a.orientation;
			const FQuat q_b = node_b.orientation;
			const FQuat q_b_mirror = mirror(q_b);

			const FQuat q_m = FQuat::Slerp(q_a, q_b_mirror, 0.5f);

			const FVector p_a = node_a.position;
			const FVector p_b = node_b.position;

			const FVector p_m = (p_b - p_a) * 0.5f + p_a;

			const FVector n_m = q_m.RotateVector(FVector::ForwardVector);

			const FVector l_a = FMath::ClosestPointOnInfiniteLine(p_m, p_m + n_m, p_a);
			const FVector l_b = FMath::ClosestPointOnInfiniteLine(p_m, p_m + n_m, p_b);

			auto velocityObject_a = velocities.componentPtrForNode(node_a.handle());
			auto velocityObject_b = velocities.componentPtrForNode(node_b.handle());

			if (!velocityObject_a || !velocityObject_b) return;

			const FVector v_a = velocityObject_a->linearVelocity;
			const FVector v_b = velocityObject_b->linearVelocity;

			const float w_a = 20.0f * particle_a.interactionMultiplier;
			const float w_b = 20.0f * particle_b.interactionMultiplier;

			// spring to mid-line
			//const FVector vToLine_a = PhysicsUtilities::criticallyDampedSpringVelocity(v_a, p_a, l_a, w_a, deltaT);
			//const FVector vToLine_b = PhysicsUtilities::criticallyDampedSpringVelocity(v_b, p_b, l_b, w_b, deltaT);

			//FVector dv_a = (vToLine_a - v_a) * deltaT;
			//FVector dv_b = (vToLine_b - v_b) * deltaT;

			FVector dv_a = FVector::ZeroVector;
			FVector dv_b = FVector::ZeroVector;

			// spring to binding offset
			//const FVector lRest_a = l_b + (l_a - l_b).GetSafeNormal() * bindingOffset;
			//const FVector lRest_b = l_a + (l_b - l_a).GetSafeNormal() * bindingOffset;
			//const FVector vToRest_a = PhysicsUtilities::criticallyDampedSpringVelocity(v_a, l_a, lRest_a, w_a, deltaT);
			//const FVector vToRest_b = PhysicsUtilities::criticallyDampedSpringVelocity(v_b, l_b, lRest_b, w_b, deltaT);

			const FVector lRest_a = p_a + (p_b - p_a).GetSafeNormal() * bindingOffset;
			const FVector lRest_b = p_b + (p_a - p_b).GetSafeNormal() * bindingOffset;

			const FVector vToRest_a = PhysicsUtilities::criticallyDampedSpringVelocity(v_a, p_a, lRest_a, w_a, deltaT);
			const FVector vToRest_b = PhysicsUtilities::criticallyDampedSpringVelocity(v_b, p_b, lRest_b, w_b, deltaT);

			dv_a += (vToRest_a - v_a) * deltaT;
			dv_b += (vToRest_b - v_b) * deltaT;

			// apply velocities to the rigid bodies
			const float maxDv = 20.0f * deltaT;

			dv_a.GetClampedToMaxSize(maxDv);
			dv_b.GetClampedToMaxSize(maxDv);

			rigidHandle_a(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
					velocity->linearVelocity += dv_a;
			});

			rigidHandle_b(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
					velocity->linearVelocity += dv_b;
			});

			// apply rotational forces
			const FQuat dq_a = q_m * q_a.Inverse();
			const FQuat dq_b = q_m * q_b_mirror.Inverse();

			float torqueMagnitude_a;
			float torqueMagnitude_b;

			FVector torque_a = _torqueSpring(dq_a, torqueMagnitude_a) * 100.0f * (1 / float(rigidCount_a));
			FVector torque_b = _torqueSpring(dq_b, torqueMagnitude_b) * 100.0f * (1 / float(rigidCount_b));

			rigidHandle_a(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
				{
					const FVector r = (nodeHandle(*graph).position - com_b);

					velocity->linearVelocity += FVector::CrossProduct(torque_a, r) * deltaT;
				}
			});

			rigidHandle_b(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle nodeHandle, FFlexRigidBodyConnection& rigidConnection) {
				if (auto velocity = velocities.componentPtrForNode(nodeHandle))
				{
					const FVector r = (nodeHandle(*graph).position - com_b);

					velocity->linearVelocity += FVector::CrossProduct(torque_b, r) * deltaT;
				}
			});

			// binding
			const float bindingAngleThreshold = 0.5f;
			const float bindingThresholdSqrd = std::pow(1.5f, 2.0f);

			const FVector target_a = p_b + q_b.RotateVector(FVector::ForwardVector) * bindingOffset;
			const FVector target_b = p_a + q_a.RotateVector(FVector::ForwardVector) * bindingOffset; 

			rotationallyAligned = std::abs(torqueMagnitude_a) < 1.5f && std::abs(torqueMagnitude_b) < bindingAngleThreshold;
			positionallyAligned =  FVector::DistSquared(p_a, target_a) < bindingThresholdSqrd
				&& FVector::DistSquared(p_b, target_b) < bindingThresholdSqrd;

			if (positionallyAligned && rotationallyAligned)
			{
				graph->beginTransaction();

				FMatrix trans = FTranslationMatrix::Make(-p_b)*FRotationMatrix::Make(q_b.Inverse())* FRotationMatrix::Make(mirror(q_a)) * FTranslationMatrix::Make(target_b);

				for (FGraphNodeHandle b_ : particle_b.subs)
				{
					if (!b_) continue;

					FVector aligned_b = trans.TransformPosition(b_(*graph).position);

					for (FGraphNodeHandle a_ : particle_a.subs)
					{
						if (!a_) continue;

						float dist = (a_(*graph).position - aligned_b).Size();
						float strength = a_ == b_ ? 0.5f : 0.1f;

						graph->connectNodes<FFlexConnection>(a_, b_, dist, strength);
					}
				}

				graph->endTransaction();

				// we've mutated the graph, make sure we updated our particle by handle
				graph->componentPtr<FMLOrientedParticle>(particle_a.nodeHandle())->numActiveBonds++;
				graph->componentPtr<FMLOrientedParticle>(particle_b.nodeHandle())->numActiveBonds++;

				recentlyBonded.Add(particle_a.nodeHandle());
				recentlyBonded.Add(particle_b.nodeHandle());
			}
		}


	}
}

void UMLParticleSimulation::tick_bindOrientedParticles_greedyFlexConnections(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();
	auto& aggregateIDs = graph->componentStorage<FMLAggregateNO_id>();
	auto& meshes = graph->componentStorage<FGraphMesh>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIDs.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};




	// find and cache our helper particles
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		if (particle._didInit)
			return;

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		particle._rigidBody = rigidHandle;

		std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;
		rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
			float distSqrd = (helperHandle(*graph).position - node.position).SizeSquared();

			helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
		});

		std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
			return a.second < b.second;
		});

		size_t n = particle.subs.size();

		int i;
		for (i = 0; i < helperHandles.size() && i < n; ++i)
		{
			particle.subs[i] = helperHandles[i].first;
		}

		// null the rest
		for (; i < n; ++i)
			particle.subs[i] = FGraphNodeHandle::null;

		particle._didInit = true;
	});

	TSet< InteractionPair > interactions;

	// update our interaction partners
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node_a, FMLOrientedParticle& particle_a)
	{
		if (particle_a.numActiveBonds >= hackMaxNumBonds)
			return;

		FGraphNodeHandle aggregateID_a = getAggregateID(node_a.handle());

		unrealAABB::AABB query(node_a.position, hackNeighbourhoodRadius);

		float radius_a = hackNeighbourhoodRadius * particle_a.radiusMultiplier;

		FGraphNodeHandle nearest_b;
		float nearestDistSqrd = std::numeric_limits<float>::max();

		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& particle_b = orientedParticles.at(particleIndex);

			const float radius_b = hackNeighbourhoodRadius * particle_b.radiusMultiplier;
			const float c = _attractionMatrix(particle_a.speciesID, particle_b.speciesID);
			const FGraphNodeHandle aggregateID_b = getAggregateID(particle_b.nodeHandle());

			FGraphNode& node_b = graph->node(particle_b.nodeHandle());
			float distSqrd = (node_b.position - node_a.position).SizeSquared();

			if (particle_b.numActiveBonds >= hackMaxNumBonds
				|| c <= 0.0f
				|| particle_b.nodeHandle() == particle_a.nodeHandle()
				|| aggregateID_a == aggregateID_b
				|| !_facing(node_a, node_b)
				|| !_aligned(node_a, particle_a, node_b, particle_b)
				|| distSqrd >= std::pow(std::max(radius_b, radius_a), 2.0f)
				|| distSqrd >= nearestDistSqrd )
			{
				return true;
			}

			// we already have a bond to this guy, no double bonds
			InteractionPair bondPair(particle_a._rigidBody, particle_b._rigidBody);
			if (_numBonds.FindOrAdd(bondPair) > 0)
				return true;

			nearest_b = particle_b.nodeHandle();
			nearestDistSqrd = distSqrd;

			return true;
		});

		if (nearest_b)
		{
			// we've mutated the graph, make sure we updated our particle by handle
			graph->componentPtr<FMLOrientedParticle>(particle_a.nodeHandle())->numActiveBonds++;
			graph->componentPtr<FMLOrientedParticle>(nearest_b)->numActiveBonds++;

			auto interactionPair = InteractionPair{ particle_a.nodeHandle(), nearest_b };
			interactions.Emplace(interactionPair);

			_numBonds.FindOrAdd(interactionPair)++;
		}

		return;
	});

	// Sort, so that we can visit some of the L1 and L2 cache in order
	interactions.Sort([&](const InteractionPair& a, const InteractionPair& b) {
		return a.a < b.a;
	});

	// Calculate our forces for each interaction pair
	for (auto& pair : interactions)
	{
		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		const FQuat q_a = node_a.orientation;
		const FQuat q_b = node_b.orientation;

		const FVector p_a = node_a.position;
		const FVector p_b = node_b.position;

		FMLOrientedParticle& particle_a = orientedParticles.componentForNode(pair.a);
		FMLOrientedParticle& particle_b = orientedParticles.componentForNode(pair.b);

		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());
		auto rigidHandle_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_b.handle());

		// count the particles and find the center of mass
		int32 rigidCount_a = 0;
		FVector com_a = FVector::ZeroVector;
		centerOfMass(rigidHandle_a, rigidCount_a, com_a);
		if (rigidCount_a == 0)
			return;

		int32 rigidCount_b = 0;
		FVector com_b = FVector::ZeroVector;
		centerOfMass(rigidHandle_b, rigidCount_b, com_b);
		if (rigidCount_b == 0)
			return;

		const FVector target_b = p_a + q_a.RotateVector(FVector::ForwardVector) * bindingOffset;

		node_a.scale = debugBindingSiteScale;
		node_b.scale = debugBindingSiteScale;

		graph->beginTransaction();

		FMatrix trans = FTranslationMatrix::Make(-p_b)*FRotationMatrix::Make(q_b.Inverse())* FRotationMatrix::Make(mirror(q_a)) * FTranslationMatrix::Make(target_b);

		for (FGraphNodeHandle b_ : particle_b.subs)
		{
			if (!b_) continue;

			FVector aligned_b = trans.TransformPosition(b_(*graph).position);

			for (FGraphNodeHandle a_ : particle_a.subs)
			{
				if (!a_) continue;

				float dist = (a_(*graph).position - aligned_b).Size();
				float strength = 0.1f;

				auto edgeHandle = graph->connectNodes<FFlexConnection>(a_, b_, dist, strength);

				auto& orientedConnection = graph->addOrReplaceEdgeObject<FMLOrientedConnection>(edgeHandle);
				orientedConnection.targetDistance = dist;
				orientedConnection.breakingDistance = dist * 3.0f;
			}
		}

		graph->endTransaction();


		recentlyBonded.Add(particle_a.nodeHandle());
		recentlyBonded.Add(particle_b.nodeHandle());
	}

	// Update connections
	auto flexConnections = graph->edgeStorage<FFlexConnection>();

	graph->beginTransaction();
	graph->edgeView<FMLOrientedConnection>().eachHandle([&](FMLOrientedConnection& connection, FGraphEdgeHandle handle) {
		FGraphEdge& edge = graph->edge(handle);

		FFlexConnection * flexConnection = flexConnections.objectPtr(handle);

		float distSqrd = FVector::DistSquared(graph->node(edge.a).position, graph->node(edge.b).position);
		float dist = std::sqrt(distSqrd);

		float length = (dist - connection.targetDistance) / 2 + connection.targetDistance;

		flexConnection->length = length;

		//if (dist > connection.breakingDistance)
		//{
		//	graph->removeConnection(handle);
		//}
	});
	graph->endTransaction();
}

void UMLParticleSimulation::tick_bindOrientedParticles4(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();
	auto& aggregateIDs = graph->componentStorage<FMLAggregateNO_id>();
	auto& meshes = graph->componentStorage<FGraphMesh>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIDs.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};

	// find and cache our helper particles
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		if (particle._didInit)
			return;

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		particle._rigidBody = rigidHandle;

		std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;
		rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
			float distSqrd = (helperHandle(*graph).position - node.position).SizeSquared();

			helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
		});

		std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
			return a.second < b.second;
		});

		size_t n = particle.subs.size();

		int i;
		for (i = 0; i < helperHandles.size() && i < n; ++i)
		{
			particle.subs[i] = helperHandles[i].first;
		}

		// null the rest
		for (; i < n; ++i)
			particle.subs[i] = FGraphNodeHandle::null;

		particle._didInit = true;
	});

	TSet< InteractionPair > interactions;

	// update our interaction partners
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node_a, FMLOrientedParticle& particle_a)
	{
		if (particle_a.numActiveBonds >= hackMaxNumBonds)
			return;

		FGraphNodeHandle aggregateID_a = getAggregateID(node_a.handle());

		unrealAABB::AABB query(node_a.position, hackNeighbourhoodRadius);

		float radius_a = hackNeighbourhoodRadius * particle_a.radiusMultiplier;

		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& particle_b = orientedParticles.at(particleIndex);

			const float radius_b = hackNeighbourhoodRadius * particle_b.radiusMultiplier;
			const float c = _attractionMatrix(particle_a.speciesID, particle_b.speciesID);
			const FGraphNodeHandle aggregateID_b = getAggregateID(particle_b.nodeHandle());

			FGraphNode& node_b = graph->node(particle_b.nodeHandle());
			float distSqrd = (node_b.position - node_a.position).SizeSquared();

			if (particle_b.numActiveBonds >= hackMaxNumBonds
				|| c <= 0.0f
				|| particle_b.nodeHandle() == particle_a.nodeHandle()
				|| aggregateID_a == aggregateID_b
				|| distSqrd >= std::pow(std::max(radius_b, radius_a), 2.0f)
				|| !_facing(node_a, node_b)
				)
			{
				return true;
			}

			// we already have a bond to this guy, no double bonds
			InteractionPair bondPair(particle_a._rigidBody, particle_b._rigidBody);
			if (_numBonds.FindOrAdd(bondPair) > 0)
				return true;

			auto interactionPair = InteractionPair{ particle_a.nodeHandle(), particle_b.nodeHandle() };
			interactions.Emplace(interactionPair);

			return true;
		});

		return;
	});

	// Sort, so that we can visit some of the L1 and L2 cache in order
	interactions.Sort([&](const InteractionPair& a, const InteractionPair& b) {
		return a.a < b.a;
	});

	// Calculate our forces for each interaction pair
	for (auto& pair : interactions)
	{
		FMLOrientedParticle& particle_a = orientedParticles.componentForNode(pair.a);
		FMLOrientedParticle& particle_b = orientedParticles.componentForNode(pair.b);

		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		const FQuat q_a = node_a.orientation;
		const FQuat q_b = node_b.orientation;

		const FVector p_a = node_a.position;
		const FVector p_b = node_b.position;

		const FVector target_b = p_a + q_a.RotateVector(FVector::ForwardVector) * bindingOffset;

		FMatrix trans = FTranslationMatrix::Make(-p_b)*FRotationMatrix::Make(q_b.Inverse())* FRotationMatrix::Make(mirror(q_a)) * FTranslationMatrix::Make(target_b);


		for (int i = 0; i < 4; ++i)
		{
			FGraphNodeHandle a_ = particle_a.subs[i];
			FGraphNodeHandle b_ = particle_b.subs[i];

			FGraphNode& node_a_ = graph->node(a_);
			FGraphNode& node_b_ = graph->node(b_);

			const FVector p = node_b_.position;
			const FVector p_b_ = trans.TransformPosition(node_b_.position);


			FVector dir = p_b_ - p;

			const float sizeSqrd = dir.SizeSquared();
			dir = dir.GetSafeNormal();

			const float c = 1000.0f;

			const float acceleration = std::min(16.0f, c * (1 / sizeSqrd));

			FVelocityGraphObject& velocity_a = velocities.componentForNode(pair.a);
			FVelocityGraphObject& velocity_b = velocities.componentForNode(pair.b);

			velocity_a.linearVelocity -= dir * acceleration * deltaT;
			velocity_b.linearVelocity += dir * acceleration * deltaT;
		}
	}

	// update our interaction partners
	TSet<InteractionPair> bonds;
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node_a, FMLOrientedParticle& particle_a)
	{
		if (particle_a.numActiveBonds >= hackMaxNumBonds)
			return;

		FGraphNodeHandle aggregateID_a = getAggregateID(node_a.handle());

		unrealAABB::AABB query(node_a.position, hackNeighbourhoodRadius);

		float radius_a = hackNeighbourhoodRadius * particle_a.radiusMultiplier;

		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& particle_b = orientedParticles.at(particleIndex);

			const float radius_b = hackNeighbourhoodRadius * particle_b.radiusMultiplier;
			const float c = _attractionMatrix(particle_a.speciesID, particle_b.speciesID);
			const FGraphNodeHandle aggregateID_b = getAggregateID(particle_b.nodeHandle());

			FGraphNode& node_b = graph->node(particle_b.nodeHandle());
			float distSqrd = (node_b.position - node_a.position).SizeSquared();

			if (particle_b.numActiveBonds >= hackMaxNumBonds
				|| c <= 0.0f
				|| particle_b.nodeHandle() == particle_a.nodeHandle()
				|| aggregateID_a == aggregateID_b
				|| distSqrd >= std::pow(std::max(radius_b, radius_a), 2.0f)
				|| distSqrd >= std::pow(flexBindingRadius, 2.0f)
				|| !_facing(node_a, node_b)
				|| !_aligned(node_a, particle_a, node_b, particle_b)
				)
			{
				return true;
			}

			// we already have a bond to this guy, no double bonds
			InteractionPair bondPair(particle_a._rigidBody, particle_b._rigidBody);
			if (_numBonds.FindOrAdd(bondPair) > 0)
				return true;

			// we've mutated the graph, make sure we updated our particle by handle
			graph->componentPtr<FMLOrientedParticle>(particle_a.nodeHandle())->numActiveBonds++;
			graph->componentPtr<FMLOrientedParticle>(particle_b.nodeHandle())->numActiveBonds++;

			auto interactionPair = InteractionPair{ particle_a.nodeHandle(), particle_b.nodeHandle() };
			bonds.Emplace(interactionPair);

			_numBonds.FindOrAdd(interactionPair)++;

			return true;
		});

		return;
	});

	// Calculate our forces for each interaction pair
	for (auto& pair : bonds)
	{
		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		const FQuat q_a = node_a.orientation;
		const FQuat q_b = node_b.orientation;

		const FVector p_a = node_a.position;
		const FVector p_b = node_b.position;

		FMLOrientedParticle& particle_a = orientedParticles.componentForNode(pair.a);
		FMLOrientedParticle& particle_b = orientedParticles.componentForNode(pair.b);

		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());
		auto rigidHandle_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_b.handle());

		// count the particles and find the center of mass
		int32 rigidCount_a = 0;
		FVector com_a = FVector::ZeroVector;
		centerOfMass(rigidHandle_a, rigidCount_a, com_a);
		if (rigidCount_a == 0)
			return;

		int32 rigidCount_b = 0;
		FVector com_b = FVector::ZeroVector;
		centerOfMass(rigidHandle_b, rigidCount_b, com_b);
		if (rigidCount_b == 0)
			return;


		const FVector target_b = p_a + q_a.RotateVector(FVector::ForwardVector) * bindingOffset;

		node_a.scale = debugBindingSiteScale;
		node_b.scale = debugBindingSiteScale;

		graph->beginTransaction();

		FMatrix trans = FTranslationMatrix::Make(-p_b)*FRotationMatrix::Make(q_b.Inverse())* FRotationMatrix::Make(mirror(q_a)) * FTranslationMatrix::Make(target_b);

		int b_i = 0;
		for (FGraphNodeHandle b_ : particle_b.subs)
		{
			if (!b_) continue;

			FVector aligned_b = trans.TransformPosition(b_(*graph).position);

			int a_i = 0;
			for (FGraphNodeHandle a_ : particle_a.subs)
			{
				if (!a_) continue;

				float dist = (a_(*graph).position - aligned_b).Size();

				auto edgeHandle = graph->connectNodes<FFlexConnection>(a_, b_, dist, a_i == b_i ? baseBindingStrength : baseBindingStrength * 0.5f);

				auto& orientedConnection = graph->addOrReplaceEdgeObject<FMLOrientedConnection>(edgeHandle);
				orientedConnection.targetDistance = dist;
				orientedConnection.breakingDistance = dist * 3.0f;

				a_i++;
			}

			b_i++;
		}

		graph->endTransaction();


		recentlyBonded.Add(particle_a.nodeHandle());
		recentlyBonded.Add(particle_b.nodeHandle());
	}

	// Update connections
	auto flexConnections = graph->edgeStorage<FFlexConnection>();

	//graph->beginTransaction();
	//graph->edgeView<FMLOrientedConnection>().eachHandle([&](FMLOrientedConnection& connection, FGraphEdgeHandle handle) {
	//	FGraphEdge& edge = graph->edge(handle);

	//	FFlexConnection * flexConnection = flexConnections.objectPtr(handle);

	//	float distSqrd = FVector::DistSquared(graph->node(edge.a).position, graph->node(edge.b).position);
	//	float dist = std::sqrt(distSqrd);

	//	float length = (dist - connection.targetDistance) * .95f + connection.targetDistance;

	//	flexConnection->length = length;

	//	//if (dist > connection.breakingDistance)
	//	//{
	//	//	graph->removeConnection(handle);
	//	//}
	//});
	//graph->endTransaction();
}

void UMLParticleSimulation::tick_bindOrientedParticles_sortedInteractions(float deltaT)
{
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();
	auto& aggregateIDs = graph->componentStorage<FMLAggregateNO_id>();
	auto& meshes = graph->componentStorage<FGraphMesh>();

	const float hackMaxNumBonds = 1;

	auto getAggregateID = [&](FGraphNodeHandle handle) -> FGraphNodeHandle {
		auto aggregate = aggregateIDs.componentPtrForNode(handle);

		return aggregate->aggregateNode ? aggregate->aggregateNode : FGraphNodeHandle::null;
	};

	// find and cache our helper particles
	_cacheBindingSubs();


	TSet< InteractionPair > interactions;

	// update our interaction partners
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode& node_a, FMLOrientedParticle& particle_a)
	{
		if (particle_a.numActiveBonds >= hackMaxNumBonds
			|| !particle_a.active)
			return;

		FGraphNodeHandle aggregateID_a = getAggregateID(node_a.handle());

		unrealAABB::AABB query(node_a.position, hackNeighbourhoodRadius);


		_orientedParticleBVH.query(query, [&](unsigned int particleIndex) {
			FMLOrientedParticle& particle_b = orientedParticles.at(particleIndex);

			const float b = _bindingRadiusMatrix(particle_a.speciesID, particle_b.speciesID) * std::max(particle_a.radiusMultiplier, particle_b.radiusMultiplier);

			const FGraphNodeHandle aggregateID_b = getAggregateID(particle_b.nodeHandle());

			FGraphNode& node_b = graph->node(particle_b.nodeHandle());
			float distSqrd = (node_b.position - node_a.position).SizeSquared();

			if (particle_b.numActiveBonds >= hackMaxNumBonds
				|| (!particle_b.active && !particle_b.isSlave)
				|| particle_b.nodeHandle() == particle_a.nodeHandle()
				|| aggregateID_a == aggregateID_b
				|| (!_facing(node_a, node_b) && (particle_a.speciesID == 0 || particle_a.speciesID == 1)) 
				//|| (!_aligned(node_a, particle_a, node_b, particle_b) && (particle_a.speciesID == )
				|| distSqrd >= std::pow(b, 2.0f)
				)
			{
				return true;
			}

			auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());
			auto rigidHandle_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_b.handle());

			auto& chain_a = _chainID.FindOrAdd(rigidHandle_a);
			auto& chain_b = _chainID.FindOrAdd(rigidHandle_b);

			// oh, yes, huge fucking hack
			// don't bind to the same chain
			//if( (particle_a.speciesID == 0 || particle_a.speciesID == 1)
			//	&& (chain_a > 0 && chain_a == chain_b) )
			//{
			//	return true;
			//}

			InteractionPair rigidInteraction(rigidHandle_a, rigidHandle_b);
			auto numBonds = _numBonds.Find(rigidInteraction);
			if (numBonds && *numBonds > 0)
				return true;

			auto interactionPair = InteractionPair{ particle_a.nodeHandle(), particle_b.nodeHandle() };
			interactions.Emplace(interactionPair);

			return true;
		});
	});

	// Sort, ascending. Smallest first.
	interactions.Sort([&](const InteractionPair& a, const InteractionPair& b) {
		float dsqrd_a = FVector::DistSquared(graph->node(a.a).position, graph->node(a.b).position);
		float dsqrd_b = FVector::DistSquared(graph->node(b.a).position, graph->node(b.b).position);

		return dsqrd_a < dsqrd_b;
	});

	TSet<FGraphNodeHandle> chosen;

	// Calculate our forces for each interaction pair
	for (auto& pair : interactions)
	{
		// only bind to the nearest
		if (chosen.Contains(pair.a) || chosen.Contains(pair.b))
			continue;

		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		const FQuat q_a = node_a.orientation;
		const FQuat q_b = node_b.orientation;

		const FVector p_a = node_a.position;
		const FVector p_b = node_b.position;

		FMLOrientedParticle& particle_a = orientedParticles.componentForNode(pair.a);
		FMLOrientedParticle& particle_b = orientedParticles.componentForNode(pair.b);

		const FVector target_b = p_a + q_a.RotateVector(FVector::ForwardVector) * bindingOffset;

		node_a.scale = debugBindingSiteScale;
		node_b.scale = debugBindingSiteScale;

		graph->beginTransaction();

		FMatrix trans = FTranslationMatrix::Make(-p_b)*FRotationMatrix::Make(q_b.Inverse())* FRotationMatrix::Make(mirror(q_a)) * FTranslationMatrix::Make(target_b);

		for (FGraphNodeHandle b_ : particle_b.subs)
		{
			if (!b_) continue;

			FVector aligned_b = trans.TransformPosition(b_(*graph).position);

			for (FGraphNodeHandle a_ : particle_a.subs)
			{
				if (!a_) continue;

				float dist = (a_(*graph).position - aligned_b).Size();

				auto edgeHandle = graph->connectNodes<FFlexConnection>(a_, b_, dist, baseBindingStrength);

				auto& orientedConnection = graph->addOrReplaceEdgeObject<FMLOrientedConnection>(edgeHandle);
				orientedConnection.targetDistance = dist;
				orientedConnection.breakingDistance = dist * 3.0f;
			}
		}

		graph->endTransaction();

		// we've mutated the graph, make sure we updated our particle by handle
		graph->componentPtr<FMLOrientedParticle>(particle_a.nodeHandle())->numActiveBonds++;
		graph->componentPtr<FMLOrientedParticle>(particle_b.nodeHandle())->numActiveBonds++;

		particle_a._bondedPartner = particle_b.nodeHandle();
		particle_b._bondedPartner = particle_a.nodeHandle();

		recentlyBonded.Add(particle_a.nodeHandle());
		recentlyBonded.Add(particle_b.nodeHandle());

		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());
		auto rigidHandle_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_b.handle());

		InteractionPair rigidInteraction(rigidHandle_a, rigidHandle_b);
		_numBonds.FindOrAdd(rigidInteraction)++;

		// update chains
		auto& chain_a = _chainID.FindOrAdd(rigidHandle_a);
		auto& chain_b = _chainID.FindOrAdd(rigidHandle_b);

		if (chain_a == 0 && chain_b == 0)
		{
			chain_a = _nextChain;
			chain_b = _nextChain;

			_nextChain++;
		}
		else
		{
			auto id = std::max(chain_a, chain_b);
			auto other = std::min(chain_a, chain_b);

			// assign the id to both chains
			for (auto& pair : _chainID)
			{
				if( pair.Value == other ) pair.Value = id;
			}
		}
	}

	// Update connections
	auto flexConnections = graph->edgeStorage<FFlexConnection>();

	graph->beginTransaction();
	graph->edgeView<FMLOrientedConnection>().eachHandle([&](FMLOrientedConnection& connection, FGraphEdgeHandle handle) {
		FGraphEdge& edge = graph->edge(handle);

		FFlexConnection * flexConnection = flexConnections.objectPtr(handle);

		float distSqrd = FVector::DistSquared(graph->node(edge.a).position, graph->node(edge.b).position);
		float dist = std::sqrt(distSqrd);

		float length = (dist - connection.targetDistance) / 2 + connection.targetDistance;

		flexConnection->length = length;

		//if (dist > connection.breakingDistance)
		//{
		//	graph->removeConnection(handle);
		//}
	});
	graph->endTransaction();
}


void UMLParticleSimulation::_cacheBindingSubs()
{
	graph->each_node_object<FMLOrientedParticle>([&](FGraphNode node, FMLOrientedParticle& particle)
	{
		if (particle._didInit)
			return;

		auto rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node.handle());

		particle._rigidBody = rigidHandle;

		std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;
		rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
			float distSqrd = (helperHandle(*graph).position - node.position).SizeSquared();

			helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
		});

		std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
			return a.second < b.second;
		});

		size_t n = particle.subs.size();

		int i;
		for (i = 0; i < helperHandles.size() && i < n; ++i)
		{
			particle.subs[i] = helperHandles[i].first;
		}

		// null the rest
		for (; i < n; ++i)
			particle.subs[i] = FGraphNodeHandle::null;

		particle._didInit = true;
	});
}

FVector UMLParticleSimulation::_torqueSpring(const FQuat &dq, float& torqueMagnitude_out)
{
	FVector torqueAxis;

	dq.ToAxisAndAngle(torqueAxis, torqueMagnitude_out);

	torqueMagnitude_out = std::sin(torqueMagnitude_out);

	return torqueAxis * torqueMagnitude_out;
}

void UMLParticleSimulation::centerOfMass(FGraphNodeHandle rigidHandle, int32 &rigidParticleCount, FVector &com)
{
	rigidHandle(*graph).each<FFlexRigidBodyConnection>(*graph, [&](
		FGraphNodeHandle handle,
		FFlexRigidBodyConnection& connection)
	{
		rigidParticleCount++;
		com += graph->node(handle).position;
	});

	if (rigidParticleCount == 0)
		return;

	com /= float(rigidParticleCount);
}

void UMLParticleSimulation::_updateParticleBVH()
{
	{
		_particleBVH.removeAll();

		auto& particles = graph->componentStorage<FMLParticle>();

		for (FMLParticle& particle : particles)
		{
			if (!particle.isValid())
				continue;

			FGraphNode& node = graph->node(particle.nodeHandle());

			if (!node.isValid())
				continue;

			auto componentIndex = particles.elementIndex(particle.nodeHandle());

			_particleBVH.insertParticle(componentIndex, node.position, _hackRadius);
		}
	}

	{
		_orientedParticleBVH.removeAll();

		auto& orientedParticles = graph->componentStorage<FMLOrientedParticle>();

		for (FMLOrientedParticle& particle : orientedParticles)
		{
			if (!particle.isValid())
				continue;

			FGraphNode& node = graph->node(particle.nodeHandle());

			if (!node.isValid())
				continue;

			auto componentIndex = orientedParticles.elementIndex(particle.nodeHandle());

			_orientedParticleBVH.insertParticle(componentIndex, node.position, _hackRadius);
		}
	}
}