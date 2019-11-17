// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "PhysicsUtilities.h"
#include "Hungarian2.h"

#include <ctime>

#include "MolecularLego_Relaxation.h"

void UMLElementSimulation::tick(float deltaT)
{
	_updateBVHs();

	if (disabled)
		return;

	// we must find our rotations before finding positions
	auto rotationalInteractions = _interactions();

	_rotationalOptimization(deltaT, rotationalInteractions);

	auto positionalnteractions = _interactions();

	_positionalOptimization(deltaT, positionalnteractions);

	_previousInteractions = std::move(positionalnteractions);

	// do tick work
	for (auto& work : _tickWork)
		work();

	_tickWork.clear();
}

void UMLElementSimulation::tick_paused(float deltaT)
{
	_updateBVHs();
}

void UMLElementSimulation::attach()
{
	ruleGraph.init();

	_initElements();
	_initOccluders();

	_initCachedRules();

	graph->addComponentListener<FMLOccluder>(this);
	graph->addComponentListener<FMLElement>(this);

}

void UMLElementSimulation::begin()
{

}

auto UMLElementSimulation::_initCachedRules() -> void
{
	auto& rules = ruleGraph.componentStorage<FMLElementRule>();
	auto& elements = ruleGraph.componentStorage<FMLElement>();


	// write our rules to a compact _cachedRules vector
	cachedRules.clear();
	cachedRules.reserve(rules.size());

	for (FMLElementRule& rule : rules)
	{
		FGraphNode& ruleNode = ruleGraph.node(rule.nodeHandle());

		// build the rule
		cachedRules.emplace_back();
		Rule& cachedRule = cachedRules.back();

		cachedRule.ruleHandle = rule.nodeHandle();

		auto generateRulePosition = [&](FGraphNodeHandle handle) {
			FGraphNode& node = ruleGraph.node(handle);
			FMLElement * element = elements.componentPtrForNode(handle);

			if (!element) return;

			RulePosition rulePosition;
			rulePosition.type = element->type;
			rulePosition.positionType = element->positionType;
			rulePosition.position = node.position;
			rulePosition.orientation = node.orientation;
			rulePosition.handle = node.handle();

			cachedRule.positions.push_back(rulePosition);
		};

		// generate a position for the rule node
		generateRulePosition(rule.nodeHandle());

		// and each of its connections
		ruleNode.each<FMLElementRuleConnection>(ruleGraph, [&](FGraphNodeHandle handle, FMLElementRuleConnection& ruleConnection) {
			generateRulePosition(handle);
		});
	}
}

void UMLElementSimulation::detach()
{
	graph->removeComponentListener<FMLOccluder>(this);
}

void AMLRuleActor::writeToGraph(FGraph& ruleGraph)
{
	FVector p = this->GetActorLocation();
	FQuat q = this->GetActorQuat();

	FQuat q_inverse = q.Inverse();

	// create a rule
	FGraphNodeHandle ruleHandle = FGraphNodeHandle(ruleGraph.addNode(p, q, 1.0f));
	{
		FMLElementRule& newRule = ruleHandle(ruleGraph).addComponent<FMLElementRule>(ruleGraph);
		newRule = rule;
		newRule.nodeIndex = ruleHandle.index;

		if (hasElement)
		{
			FMLElement& newElement = ruleHandle(ruleGraph).addComponent<FMLElement>(ruleGraph);
			newElement = element;
			newElement.nodeIndex = ruleHandle.index;
		}
	}



	TArray<USceneComponent*> childSceneComponents;
	GetRootComponent()->GetChildrenComponents(true, childSceneComponents);

	// connect the elements to the rule, from the child AMLElementActors
	for (USceneComponent * child : childSceneComponents)
	{
		AMLElementActor * elementActor = Cast<AMLElementActor>(child->GetOwner());

		if (!elementActor) continue;

		FVector newP = child->GetComponentLocation();
		FQuat newQ = child->GetComponentQuat();

		FGraphNodeHandle elementHandle = FGraphNodeHandle(ruleGraph.addNode(newP, newQ, 1.0f));

		FMLElement& newElement = elementHandle(ruleGraph).addComponent<FMLElement>(ruleGraph);
		newElement = elementActor->element;
		newElement.nodeIndex = elementHandle.index;

		ruleGraph.connectNodes<FMLElementRuleConnection>(ruleHandle, elementHandle);
	}
}

void UMLElementSimulation::componentAdded(FGraphNodeHandle handle, ComponentType type)
{
	static const ComponentType TOccluderType = componentType<FMLOccluder>();
	static const ComponentType TElementType = componentType<FMLElement>();

	if (TOccluderType == type)
	{
		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, handle);

		FMLOccluder& occluder = graph->componentStorage<FMLOccluder>().componentForNode(handle);

		occluder.transientHost = rigidHandle;

		FVector p = graph->node(handle).position;

		int32 particleIndex = handle.index;

		if (occluderBVH.containsParticle(particleIndex))
			occluderBVH.updateParticle(particleIndex, p, _hackRadius);
		else
			occluderBVH.insertParticle(particleIndex, p, _hackRadius);
	}
	else if (TElementType == type)
	{
		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, handle);

		FMLElement& element = graph->componentStorage<FMLElement>().componentForNode(handle);

		element.rigid = handle;

		FVector p = graph->node(handle).position;

		int32 particleIndex = handle.index;

		if (elementBVH.containsParticle(particleIndex))
			elementBVH.updateParticle(particleIndex, p, _hackRadius);
		else
			elementBVH.insertParticle(particleIndex, p, _hackRadius);
	}
}

auto UMLElementSimulation::addTickWork(std::function<void()> work) -> void
{
	_tickWork.push_back(work);
}

auto UMLElementSimulation::_initOccluders() -> void
{
	auto& occluders = graph->componentStorage<FMLOccluder>();

	for (auto& occluder : occluders)
	{
		if( !occluder.isValid() ) continue;

		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, occluder.nodeHandle());
		
		if( !rigidHandle ) continue;

		occluder.transientHost = rigidHandle;
	}
}

auto UMLElementSimulation::_initElements() -> void
{
	auto& elements = graph->componentStorage<FMLElement>();

	for (auto& element : elements)
	{
		if( !element.isValid() ) continue;

		FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, element.nodeHandle());

		if( !rigidHandle ) continue;

		element.rigid = rigidHandle;
	}
}


auto UMLElementSimulation::_greedyPairs(
	const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
	const Rule& rule,
	std::vector<std::pair<int16, int16>>& pairings_out) -> float
{
	using namespace std;

	float cost = 0.0f;

	int aSize = aPositions.size();
	int bSize = rule.positions.size();

	pairings_out.clear();
	pairings_out.reserve(aSize);

	float thresholdDistance = minAssignmentDistance;

	struct FastIndexSetIndex
	{
		int index;      // the index at this position
		int runStart;   // the start of a run of indices
	};
	std::vector<FastIndexSetIndex> indices(bSize + 1);

	for (int i = 0; i < bSize + 1; ++i)
		indices[i] = { i, i };

	for (int aIndex = 0; aIndex < aSize; ++aIndex)
	{
		const FVector& aPosition = aPositions[aIndex].position;
		const FQuat& aRotation = aPositions[aIndex].rotation;

		FGraphNodeHandle aHandle = aNeighbours[aIndex];
		FElementObject& aElement = aGraph.component<FElementObject>(aHandle);

		int bNearestIndex = -1;
		float minCost = numeric_limits<float>::max();


		int bIndex = indices[0].index;
		while (bIndex < bSize)
		{
			const RulePosition& rulePosition = rule.positions[bIndex];

			const FVector& bPosition = rulePosition.position;
			const FQuat& bRotation = rulePosition.orientation;

			const float d = (aPosition - bPosition).Size();

			float costSquared = d;

			costSquared += aRotation.AngularDistance(bRotation) * 5.0f;


			// add in the cost of unmatched types
			UE_LOG(LogTemp, Warning, TEXT("disabled costs"));
			//if (aElement.type != rulePosition.type)
			//	costSquared += typeCost * typeCost;


			if (costSquared < minCost && d < thresholdDistance)
			{
				minCost = d;
				bNearestIndex = bIndex;
			}

			bIndex = indices[bIndex + 1].index;
		}

		if (bNearestIndex >= 0)
		{
			// erase an index
			const FastIndexSetIndex nullIndex = { -1, -1 };

			auto& bPrevious = bNearestIndex - 1 >= 0 ? indices[bNearestIndex - 1] : nullIndex;
			auto& bNext = indices[bNearestIndex + 1];
			auto& bCur = indices[bNearestIndex];

			bCur.index = bNext.index;

			if (bPrevious.index == bNearestIndex)
			{
				indices[bCur.runStart].index = bNext.index;
				indices[bNext.index].runStart = bCur.runStart;
			}
			else
				indices[bCur.index].runStart = bNearestIndex;

			pairings_out.emplace_back(aIndex, bNearestIndex);

			cost += minCost;

			UE_LOG(LogTemp, Warning, TEXT("disabled costs"));

			//if (aElement.type != rule.positions[bNearestIndex].type)
			//	cost += typeCost;
		}
		else
			pairings_out.emplace_back(aIndex, -1);
	}

	return cost;
}

auto UMLElementSimulation::_density(const FVector p_a) -> float
{
	auto& elements = graph->componentStorage<FMLElement>();

	float density = 0.0f;

	// density
	unrealAABB::AABB query(p_a, hackNeighbourhoodRadius);

	const float inverseSupportRadius = 1.0 / minAssignmentDistance;

	elementBVH.query(query, [&](unsigned int particleIndex) {
		FGraphNodeHandle h_b(particleIndex);

		const FVector p_b = graph->node(h_b).position;

		const FVector dir = p_a - p_b;

		float w = std::exp(inverseSupportRadius * -1.0f * (dir).SizeSquared());

		density += w;

		return true; // keep querying
	});

	return density;
}

auto UMLElementSimulation::_optimalPairs(
	FGraphNode& aNode,
	const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
	FGraphNodeHandle ruleHandle, FGraphNodeHandle ruleBase,
	const RuleInteractions& previousInteractions,
	std::vector<std::pair<FGraphNodeHandle, FGraphNodeHandle>>& pairings_out) -> float
{
	auto& elements = graph->componentStorage<FMLElement>();

	FGraphNode& ruleNode = ruleGraph.node(ruleHandle);

	FGraphNode& baseNode = ruleGraph.node(ruleBase);
	FMLElement& baseElement = elements.componentForNode(ruleBase);

	FMLElement& aElement = elements.componentForNode(aNode.handle());

	std::vector<FGraphNodeHandle> ruleHandles;

	ruleNode.each<FMLElementRuleConnection>(ruleGraph, [&](FGraphNodeHandle other, FMLElementRuleConnection& connection) {
		if( other != ruleBase )
			ruleHandles.push_back(other);
	});

	const size_t n_rule = ruleHandles.size();
	const size_t n_neighbours = aNeighbours.size();

	if (n_rule == 0 || n_neighbours == 0)
		return 0.0f;

	const size_t n = std::max(n_neighbours, n_rule - 1);

	Eigen::MatrixXf cost(n, n);

	const float maxCost = 500.0f;

	cost.setConstant(maxCost);

	auto baseTransform = [&](const FVector& p) { return baseNode.orientation.RotateVector(p - baseNode.position); };

	auto baseRotate = [&](const FQuat& q) { return baseNode.orientation * q; };

	// build the cost matrix for Hungarian assignment
	int ri = 0;
	for( FGraphNodeHandle other : ruleHandles )
	{
		FGraphNode& otherNode = ruleGraph.node(other);

		const FVector p_rule = baseTransform(otherNode.position);
		const FQuat q_rule = baseRotate(otherNode.orientation);

		float angularCost = p_rule.Size() * 0.5f;

		for (int ai = 0; ai < n_neighbours; ++ai)
		{
			FVector pred = (p_rule - aPositions[ai].position);
			const float dist = (pred).Size();

			const int row = ai;
			const int col = ri;

			if (dist > minAssignmentDistance)
				cost(row, col) = maxCost;
			else
			{
				cost(row, col) = dist;
				cost(row, col) += std::sin(aPositions[ai].rotation.AngularDistance(q_rule)) * angularCost;

				bool isStrongInteraction = previousInteractions.interactionPairs.Contains(InteractionPair(aElement.nodeHandle(), aNeighbours[ai]));
				isStrongInteraction &= dist < minAssignmentDistance;

				// Strong-interaction rule:
				// The cost is reduced if two nodes have previously strongly interacted.
				if (isStrongInteraction)
					cost(row, col) *= 0.1f;
			}
		}

		ri++;
	};

	// solve the Hungarian assignment problem for the cost matrix
	auto assignment = Hungarian::solve(cost);

	// build our assignment vector
	float assignmentCost = 0.0f;

	for (int i = 0; i < aPositions.size(); ++i)
	{
		int rule_i = assignment.second[i];

		if (i >= assignment.second.size() || rule_i >= ruleHandles.size()) continue;

		const float theCost = cost(i, rule_i);

		assignmentCost += theCost;

		if( theCost >= maxCost )
			continue;

		pairings_out.emplace_back(std::make_pair(aNeighbours[i], ruleHandles[rule_i]));
	}

	return assignmentCost;
}

auto UMLElementSimulation::_optimalPairs(
	FGraphNode& aNode,
	const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
	const Rule& rule, size_t positionIndex, const RuleInteractions& previousInteractions,
	std::vector<std::pair<int16, int16>>& pairings_out) -> float
{
	const float maxCost = 500.0f;

	const size_t n_rule = rule.positions.size();
	const size_t n_neighbours = aNeighbours.size();

	auto& elements = graph->componentStorage<FMLElement>();

	FMLElement& aNodeElement = elements.componentForNode(aNode.handle());

	const RulePosition& ruleBase = rule.positions[positionIndex];

	if (aPositions.size() == 0 || rule.positions.size() == 0
		|| aNodeElement.type != ruleBase.type)
	{
		return maxCost * n_rule;
	}

	const size_t n = std::max(n_neighbours, n_rule - 1);

	Eigen::MatrixXf cost(n,n);


	cost.setConstant(maxCost);



	FQuat inverseRotation = ruleBase.orientation.Inverse();

	auto subjectTransform = [&](const FVector& p) { return inverseRotation.RotateVector(p - ruleBase.position); };

	auto subjectRotate = [&](const FQuat& q) { return inverseRotation * q; };


	for (int ri = 0; ri < n_rule; ++ri)
	{
		if (ri == positionIndex) continue;

		const FVector p_rule = subjectTransform(rule.positions[ri].position);
		const FQuat q_rule = subjectRotate(rule.positions[ri].orientation);

		const auto type_rule = rule.positions[ri].type;

		float angularCost = p_rule.Size() * 0.5f;

		for (int ai = 0; ai < n_neighbours; ++ai)
		{
			FMLElement& element = elements.componentForNode(aNeighbours[ai]);

			FVector pred = (p_rule - aPositions[ai].position);
			const float dist = (pred).Size();

			const int row = ai;
			const int col = ri < positionIndex ? ri : ri - 1;

			if (dist > minAssignmentDistance)
				cost(row, col) = maxCost;
			else if (type_rule != element.type)
				cost(row, col) = maxCost;
			else
			{
				cost(row, col) = dist;
				cost(row, col) += std::sin(aPositions[ai].rotation.AngularDistance(q_rule)) * angularCost;

				bool isStrongInteraction = previousInteractions.interactionPairs.Contains(InteractionPair(aNodeElement.nodeHandle(), aNeighbours[ai]));
				isStrongInteraction &= dist < minAssignmentDistance;

				// Strong-interaction rule:
				// The cost is reduced if two nodes have previously strongly interacted.
				if (isStrongInteraction)
					cost(row, col) *= 0.1f;
			}
		}
	}

	auto assignment = Hungarian::solve(cost);
	
	float assignmentCost = 0.0f;

	for (int i = 0; i < aPositions.size(); ++i)
	{
		int rule_i = assignment.second[i] < positionIndex ? assignment.second[i] : assignment.second[i] + 1;

		if( i >= assignment.second.size() || rule_i >= rule.positions.size() ) continue;

		assignmentCost += cost(i, assignment.second[i]);

		const FVector p_rule = subjectTransform(rule.positions[rule_i].position);
		const auto type_rule = rule.positions[rule_i].type;

		FMLElement& element = elements.componentForNode(aNeighbours[i]);

		if ((aPositions[i].position - p_rule).SizeSquared() < std::pow(minAssignmentDistance, 2.0f)
			&& element.type == type_rule )
		{
			pairings_out.emplace_back(std::make_pair(i, rule_i));
		}
	}

	return assignmentCost;
}

auto UMLElementSimulation::_ruleSatisfied(Rule& rule, std::vector<std::pair<int16, int16>> pairs) -> bool
{
	auto& positions = rule.positions;

	unsigned int ruleOccupied = 0;

	for (auto& rulePosition : positions)
	{
		ruleOccupied += (rulePosition.positionType == ERulePositionType::Occupied 
			|| rulePosition.positionType == ERulePositionType::OccoupiedBase) ? 1 : 0;
	}

	// the rule base position is always assigned, so our initial occupied count is one
	unsigned int countOccoupied = 1;

	for (auto& pair : pairs)
	{
		auto ri = pair.second;

		auto positionType = positions[ri].positionType;

		// Empty positions must not be filled
		if (positionType == ERulePositionType::Empty)
		{
			return false;
		}
		// Occupied positions must be filled
		else if (positionType == ERulePositionType::Occupied || positionType == ERulePositionType::OccoupiedBase)
		{
			countOccoupied++;
		}
		// Optional positions may be filled
		else if (positionType == ERulePositionType::Optional)
		{
		}
	}

	return countOccoupied == ruleOccupied;
}

void UMLElementSimulation::_updateBVHs()
{
	_ticksLeftBeforeRebuild--;

	//update
	_updateBVH<FMLElement>(elementBVH);
	_updateBVH<FMLOccluder>(occluderBVH);

	// rebuild
	if (_ticksLeftBeforeRebuild < 0)
	{
		std::clock_t startClock = std::clock();

		_ticksLeftBeforeRebuild = 27000;

		//if( elementBVH.nParticles() > 0 )
		//	elementBVH.rebuild();
		//if( occluderBVH.nParticles() > 0 )
		//	occluderBVH.rebuild();

		std::clock_t endClock = std::clock();

		double seconds = (endClock - startClock) / (double)CLOCKS_PER_SEC;

		UE_LOG(LogTemp, Warning, TEXT("rebuild time %f"), seconds);
	}
}




void UMLElementSimulation::_rotationalOptimization(float deltaT, std::vector<RuleInteractions>& interactions)
{
	using namespace Eigen;
	using namespace std;

	typedef SparseMatrix<float> Matrix;
	typedef Triplet<float> Triplet;

	auto& elements = graph->componentStorage<FMLElement>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();
	auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();

	TMap<FGraphNodeHandle, FGraphNodeHandle> rigidHandles;
	for (auto& element : elements)
	{
		if (!element.isValid()) continue;

		rigidHandles.Add(element.nodeHandle(), element.rigid);
	}

	std::map<FGraphNodeHandle, unsigned int> problem;
	std::vector<FGraphNodeHandle> indexToElement;

	// setup the problem map
	// we associate with each element an index in the A matrix (the problem matrix)
	for (RuleInteractions& ruleInteractions : interactions)
	{
		for (auto& keyValue : ruleInteractions.interactionPairs)
		{
			InteractionPair& interactionPair = keyValue.Key;

			auto addToProblem = [&](FGraphNodeHandle elementHandle) {
				FGraphNodeHandle* rigidHandle = rigidHandles.Find(elementHandle);

				if (!rigidHandle) return;

				auto found = problem.find(*rigidHandle);

				if (found == problem.end())
				{
					problem[*rigidHandle] = indexToElement.size();
					indexToElement.push_back(*rigidHandle);
				}
			};

			addToProblem(interactionPair.a);
			addToProblem(interactionPair.b);
		}
	}

	const auto n = indexToElement.size();

	std::vector<std::vector< PredictedPosition > > predictions(n);

	for (RuleInteractions& ruleInteractions : interactions)
	{
		_buildPredictions(ruleInteractions, problem, rigidHandles, predictions);
	}

	std::vector<FQuatAverager> predictedRotations(n);

	// build a variance for each updated element position
	// Compute the eigen vectors and values for the prediction clusters
	for (unsigned int i = 0; i < n; ++i)
	{
		auto& positions = predictions[i];

		auto& rotationAverage = predictedRotations[i];

		FGraphNode& node = indexToElement[i].node(*graph);

		rotationAverage.average(node.orientation, 0.1f);

		for( auto& prediction : positions )
		{
			float w = prediction.weight;

			rotationAverage.average(prediction.rotation, w);
		}
	}


	// update the rotations
	for (auto keyValue : problem)
	{
		unsigned int i = keyValue.second;

		FGraphNode& node_a = graph->node(keyValue.first);

		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());

		if (!rigidHandle_a) continue;

		FQuat q_new = predictedRotations[i].currentUnnormalized().GetNormalized();

		FQuat dq = q_new * node_a.orientation.Inverse();

		FVector dp = FVector::ZeroVector;

		rigids.componentPtrForNode(rigidHandle_a)->applyRotationTranslationInferVelocity(dq, dp, *graph);
	}
}

void UMLElementSimulation::_positionalOptimization(float deltaT, std::vector<RuleInteractions>& interactions)
{
	using namespace Eigen;
	using namespace std;

	typedef SparseMatrix<float> Matrix;
	typedef Triplet<float> Triplet;

	auto& elements = graph->componentStorage<FMLElement>();
	auto& velocities = graph->componentStorage<FVelocityGraphObject>();

	TMap<FGraphNodeHandle, FGraphNodeHandle> rigidHandles;
	for (auto& element : elements)
	{
		if( !element.isValid() ) continue;

		rigidHandles.Add(element.nodeHandle(), element.rigid);
	}

	std::map<FGraphNodeHandle, unsigned int> problem;
	std::vector<FGraphNodeHandle> indexToElement;

	// setup the problem map
	// we associate with each element an index in the A matrix (the problem matrix)
	for (RuleInteractions& ruleInteractions : interactions)
	{
		for (auto& keyValue : ruleInteractions.interactionPairs)
		{
			InteractionPair& interactionPair = keyValue.Key;

			auto addToProblem = [&](FGraphNodeHandle elementHandle) {
				FGraphNodeHandle* rigidHandle = rigidHandles.Find(elementHandle);

				if (!rigidHandle) return;

				auto found = problem.find(*rigidHandle);

				if (found == problem.end())
				{
					problem[*rigidHandle] = indexToElement.size();
					indexToElement.push_back(*rigidHandle);
				}
			};

			addToProblem(interactionPair.a);
			addToProblem(interactionPair.b);
		}
	}

	const auto n = indexToElement.size();

	std::vector<std::vector< PredictedPosition > > predictions(n);

	for (auto keyValue : problem)
	{
		// this first prediction is a dummy for the covariance calculation
		// for some reason we do this twice, is this a mistake?
		FGraphNodeHandle handle = keyValue.first;
		unsigned int i = keyValue.second;

		auto node = handle(*graph);

		float weight = 1.0f;

		predictions[i].emplace_back(i, eigen(node.position), node.orientation, weight);
		predictions[i].emplace_back(i, eigen(node.position), node.orientation, weight);
	}

	for (RuleInteractions& ruleInteractions : interactions)
	{
		_buildPredictions(ruleInteractions, problem, rigidHandles, predictions);
	}

	// allocate the linear system
	// allocate Ax = b
	VectorXf b[3];
	for (int c = 0; c < 3; ++c)
	{
		b[c] = VectorXf(n);
		b[c].setZero();
	}

	std::vector<Triplet> a_t[3];

	// add the initial positions
	for (auto keyValue : problem)
	{
		FGraphNode& e = graph->node(keyValue.first);
		const unsigned int i = keyValue.second;

		const float weight = 0.5;

		for (int c = 0; c < 3; ++c)
		{
			a_t[c].emplace_back(i, i, weight);
			b[c](i) += weight * e.position[c];
		}
	}

	// build a variance for each updated element position
	// Compute the eigen vectors and values for the prediction clusters
	for (unsigned int i = 0; i < n; ++i)
	{
		auto& positions = predictions[i];

		// ignore the first item, it's a dummy for the actual position used in the covariance calculation
		float normalizeTerm = 1.0f / (positions.size() - 1);
		for (int index = 1; index < positions.size(); ++index)
		{
			auto& prediction = positions[index];

			int j = prediction.elementIndex;
			FGraphNode& element = graph->node(indexToElement[j]);

			const Vector3f predictedPosition = prediction.position;

			const Vector3f predictedDirection = predictedPosition - eigen(element.position);

			float w = prediction.weight;

			if (w != w || w == std::numeric_limits<float>::infinity())
				w = 0.0f;

			for (int c = 0; c < 3; ++c)
			{
				a_t[c].emplace_back(i, i, w);
				a_t[c].emplace_back(j, j, w);
				a_t[c].emplace_back(i, j, -w);
				a_t[c].emplace_back(j, i, -w);

				b[c](j) -= w * predictedDirection(c);
				b[c](i) += w * predictedDirection(c);
			}
		}
	}

	// finally, we can build A
	Matrix a[3] = { Matrix(n,n), Matrix(n,n), Matrix(n,n) };
	for (int c = 0; c < 3; ++c)
	{
		a[c].setFromTriplets(a_t[c].begin(), a_t[c].end());
		a[c].finalize();
		a[c].makeCompressed();
	}

	VectorXf x[3] = { VectorXf(n), VectorXf(n), VectorXf(n) };

	for (int c = 0; c < 3; ++c)
	{
		Eigen::SimplicialCholesky<Matrix> cholesky(a[c]);
		x[c] = cholesky.solve(b[c]);
	}

	// update the positions
	for (auto keyValue : problem)
	{
		unsigned int i = keyValue.second;

		FGraphNode& node_a = graph->node(keyValue.first);

		auto rigidHandle_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, node_a.handle());

		if (!rigidHandle_a) continue;

		Vector3f ep;
		for (int c = 0; c < 3; ++c)
			ep(c) = x[c](i);

		FVector p = unreal(ep);

		FVector dp = p - node_a.position;

		rigidHandle_a(*graph).component<FFlexRigidBodyObject>(*graph).applyRotationTranslationInferVelocity(FQuat::Identity, dp, *graph);
	}
}


auto UMLElementSimulation::_ruleOccluded(
	int i_rule, 
	int i_base, 
	FGraphNode& node,
	std::vector<FGraphNodeHandle>& handles,
	const std::vector<std::pair<int16, int16>>& pairings
) -> bool
{
	auto& rule = cachedRules[i_rule];
	auto& occluders = graph->componentStorage<FMLOccluder>();

	const FVector p_base = rule.positions[i_base].position;
	const FQuat qInv_base = rule.positions[i_base].orientation.Inverse();

	for (auto& pair : pairings)
	{
		FGraphNodeHandle handle = handles[pair.first];

		const FVector pred = qInv_base.RotateVector(rule.positions[pair.second].position - p_base);
		const FVector p = node.orientation.RotateVector(pred) + node.position;

		unrealAABB::AABB query(p, 0.1f * minAssignmentDistance);

		bool hit = false;

		// check if the rule application overlaps with any occluders
		// - we don't overlap with the target of the rule
		// - we might need to exclude all of the handles in the rule
		occluderBVH.query(query, [&](unsigned int particleIndex) {
			FMLOccluder& occluder = occluders.at(particleIndex);

			if (occluder.transientHost == handle)
				return true;

			hit = true;
			return false;
		});

		if (hit)
		{
			return true;
		}
	}

	return false;
}

//std::vector<FGraphNodeHandle> UMLElementSimulation::_nonOverlappingElements()
//{
//	auto& elements = graph->componentStorage<FMLElement>();
//
//	std::vector<FGraphNodeHandle> filterdElements;
//	TSet<FGraphNodeHandle> removed; 
//
//	auto weight = [&](FGraphNodeHandle handle, RuleInteractions) {
//		if( )
//	};
//
//	for (FMLElement& e_a : elements)
//	{
//		if( removed.Contains(e_a.nodeHandle()) ) continue;
//
//		FGraphNode& node_a = graph->node(e_a.nodeHandle());
//
//		// find overlaps
//		unrealAABB::AABB query(&node_a.position.X, minAssignmentDistance);
//
//		float weight_a = _previousInteractions.
//
//		_elementBVH.query(query, [&](unsigned int particleIndex) {
//			FMLElement& e_b = elements.at(particleIndex);
//
//			if (&e_a == &e_b) return true;
//			if (removed.Contains(e_b.nodeHandle())) return true;
//			
//
//			return true; // keep querying
//		});
//	}
//}

std::vector<UMLElementSimulation::RuleInteractions> UMLElementSimulation::_interactions()
{
	auto& elements = graph->componentStorage<FMLElement>();

	const auto n = cachedRules.size();

	std::vector<RuleInteractions> interactions(n);

	std::vector<PositionQuat> aPositions;
	std::vector<FGraphNodeHandle> aHandles;

	// make sure we have previous, dummy, interactions
	if (_previousInteractions.size() < n)
		_previousInteractions.resize(n);

	for (FMLElement& e_a : elements)
	{
		FGraphNode& node_a = graph->node(e_a.nodeHandle());

		FQuat q_a = node_a.orientation;
		FVector p_a = node_a.position;
		auto i_a = elements.elementIndex(e_a.nodeHandle());

		const float density = _density(p_a);

		if (density >= maxAssignmentDensity)
			continue;

		aPositions.clear();
		aHandles.clear();

		FQuat inverseRotation = q_a.Inverse();

		// find the neighborhood 
		unrealAABB::AABB query(p_a, hackNeighbourhoodRadius);

		elementBVH.query(query, [&](unsigned int particleIndex) {
			FGraphNodeHandle h_b(particleIndex);

			FMLElement& e_b = elements.componentForNode(h_b);

			if (&e_a == &e_b) return true;

			FGraphNode& node_b = graph->node(h_b);

			aPositions.push_back({
				inverseRotation.RotateVector(node_b.position - p_a),
				inverseRotation * node_b.orientation });

			aHandles.push_back(node_b.handle());

			return true; // keep querying
		});

		// find prediction vectors
		int i_rule = 0;
		for (auto& rule : cachedRules)
		{
			if (rule.positions.empty()) continue;

			RuleInteractions& ruleInteractions = interactions[i_rule];
			RuleInteractions& previousInteractions = _previousInteractions[i_rule];

			// find the min pairings
			for (int positionIndex = 0; positionIndex < rule.positions.size(); ++ positionIndex)
			{
				auto& position = rule.positions[positionIndex];

				if( position.positionType != ERulePositionType::OccoupiedBase
					||  position.type != e_a.type)
					continue;

				std::vector<std::pair<int16, int16>> pairings;

				float theCost = _optimalPairs(node_a, aPositions, aHandles, *graph, rule, positionIndex, previousInteractions, pairings);

				bool satisified = _ruleSatisfied(rule, pairings);
				bool occluded = _ruleOccluded(i_rule, positionIndex, node_a, aHandles, pairings);

				if (satisified && !occluded)
				{
					ruleInteractions.rule = &rule;
					ruleInteractions.weights.Add(e_a.nodeHandle(), theCost);

					for (auto& pair : pairings)
					{
						FGraphNodeHandle neighbourHandle = aHandles[pair.first];

						InteractionPair interactionPair(e_a.nodeHandle(), neighbourHandle);

						InteractionPairRuleIndices& ruleIndices = ruleInteractions.interactionPairs.FindOrAdd(interactionPair);

						auto& weight = e_a.nodeHandle() < neighbourHandle ? ruleIndices.first_weight : ruleIndices.second_weight;

						// reassign indices if we have a better weight
						if (theCost < weight)
						{
							// select which index to assign to
							auto& ruleIndex = e_a.nodeHandle() < neighbourHandle ? ruleIndices.first_a : ruleIndices.second_a;
							auto& ruleBase = e_a.nodeHandle() < neighbourHandle ? ruleIndices.first_base : ruleIndices.second_base;

							ruleIndex = pair.second;
							ruleBase = positionIndex;
							weight = theCost;
						}
					}
				}
			}

			i_rule++;
		}
	}

	// filter to two-way interactions
	for (auto& interaction : interactions)
		_filterInteractions(interaction);

	_normalizeInteractions(interactions);

	return interactions;
}

void UMLElementSimulation::_filterInteractions(RuleInteractions& interactions)
{
	// filtered will contain the two-way interaction pairs
	TMap<InteractionPair, InteractionPairRuleIndices> filtered;

	Rule& rule = *interactions.rule;

	// find our strong interactions
	// and count n
	TSet<FGraphNodeHandle> isStrong;
	for (auto& keyValue : interactions.interactionPairs)
	{
		InteractionPair& pair = keyValue.Key;
		InteractionPairRuleIndices& ruleIndices = keyValue.Value;

		bool aHasRule = ruleIndices.first_a >= 0;
		bool bHasRule = ruleIndices.second_a >= 0;

		if (aHasRule && bHasRule )
			//// the rule indices must reciprocate
			//&& rule.positions[pairIndices.first_a].reciprocal == pairIndices.second_a)
		{
			filtered.Add(pair, ruleIndices);

			isStrong.Add(pair.a);
			isStrong.Add(pair.b);
		}
	}

	interactions.interactionPairs = std::move(filtered);
}

void UMLElementSimulation::_normalizeInteractions(std::vector<RuleInteractions>& interactions)
{
	// find the max weight
	float max = 0.0f;

	for (RuleInteractions& ruleInteractions : interactions)
	{
		for (auto& keyValue : ruleInteractions.weights)
		{
			max = std::max(keyValue.Value, max);
		}
	}

	// normalize and invert costs, so the least cost has the highest weight
	for (RuleInteractions& ruleInteractions : interactions)
	{
		for (auto& keyValue : ruleInteractions.weights)
		{
			keyValue.Value = 1.0f - keyValue.Value / max;
		}
	}
}

void UMLElementSimulation::_buildPredictions(
	RuleInteractions& interactions,
	std::map<FGraphNodeHandle, unsigned int>& problem,
	TMap<FGraphNodeHandle, FGraphNodeHandle>& rigidHandles,
	std::vector<std::vector<PredictedPosition>>& predictions_in_out)
{
	Rule& rule = *interactions.rule;

	// We want to enumerate by neighborhood. This will keep more of the problem in the caches.
	interactions.interactionPairs.KeySort([](const InteractionPair& a, const InteractionPair& b) {
		return a.a < b.a;
	});

	for (auto& keyValue : interactions.interactionPairs)
	{
		InteractionPair& pair = keyValue.Key;
		InteractionPairRuleIndices& rulePair = keyValue.Value;

		auto emplacePrediction = [&](FGraphNode& subject, FGraphNode& predictor, int8_t ruleBase, int8_t ruleIndex) {
			auto& basePosition = rule.positions[ruleBase];
			auto& indexPosition = rule.positions[ruleIndex];

			FVector r = basePosition.orientation.UnrotateVector(indexPosition.position - basePosition.position);

			r = predictor.orientation.RotateVector(r);

			FGraphNodeHandle subjectRigidHandle = rigidHandles.FindChecked(subject.handle());
			FGraphNodeHandle predictorRigidHandle = rigidHandles.FindChecked(predictor.handle());

			FGraphNode& predictorRigid = graph->node(predictorRigidHandle);
			FGraphNode& subjectRigid = graph->node(subjectRigidHandle);


			FQuat q_new = predictorRigid.orientation * (basePosition.orientation.Inverse() * indexPosition.orientation);
			FQuat dq = q_new * subject.orientation.Inverse();

			FQuat quat_a = dq * subjectRigid.orientation;

			FVector dir = subjectRigid.position - subject.position;
			FVector pred_a = predictor.position + r + dq.RotateVector(dir);

			unsigned int i_a = problem[subjectRigidHandle];
			unsigned int i_b = problem[predictorRigidHandle];

			predictions_in_out[i_a].emplace_back(
				i_b,
				eigen(pred_a),
				quat_a,
				interactions.weights[predictor.handle()]
			);
		};

		FGraphNode& node_a = graph->node(pair.a);
		FGraphNode& node_b = graph->node(pair.b);

		if(rulePair.second_a >= 0) emplacePrediction(node_a, node_b, rulePair.second_base, rulePair.second_a);
		if(rulePair.first_a >= 0) emplacePrediction(node_b, node_a, rulePair.first_base, rulePair.first_a);
	}
}
