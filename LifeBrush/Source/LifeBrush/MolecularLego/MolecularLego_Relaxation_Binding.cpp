// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "PhysicsUtilities.h"
#include "Hungarian2.h"

#include "MolecularLego_Relaxation.h"

auto UMLElementSimulation::_bindingHandles(FGraphNodeHandle rigidHandle, FVector target, int num) -> std::vector<FGraphNodeHandle>
{
	std::vector<std::pair<FGraphNodeHandle, float>> helperHandles;

	FGraphNode& rigidNode = graph->node(rigidHandle);

	// find the distance of each sub-node from the rigidNode
	rigidNode.each<FFlexRigidBodyConnection>(*graph, [&](FGraphNodeHandle helperHandle, FFlexRigidBodyConnection& connection) {
		float distSqrd = (helperHandle(*graph).position - rigidNode.position).SizeSquared();

		helperHandles.emplace_back(std::make_pair(helperHandle, distSqrd));
	});

	// sort descending
	std::sort(helperHandles.begin(), helperHandles.end(), [&](auto& a, auto& b) {
		return a.second > b.second;
	});

	// build our result, it's the num furthest nodes
	// This should be much smarter, like finding the combination of nodes that is close to the plane orthonormal to the rigidHandle and target
	// direction vector, but far from each other.
	std::vector<FGraphNodeHandle> result;
	result.reserve(num);

	for (size_t i = 0; i < num && i < helperHandles.size(); ++i)
	{
		result.push_back(helperHandles[i].first);
	}

	return result;
}


void UMLElementSimulation::_bind(FGraphNodeHandle handle_a, FGraphNodeHandle handle_b, FVector targetPosition_b, FQuat targetRotation_b)
{
	FGraphNode& node_a = graph->node(handle_a);
	FGraphNode& node_b = graph->node(handle_b);

	FGraphNodeHandle rigid_a = FFlexRigidBodyObject::getRigidBodyHandle(*graph, handle_a);
	FGraphNodeHandle rigid_b = FFlexRigidBodyObject::getRigidBodyHandle(*graph, handle_b);

	auto bindingNodes_a = _bindingHandles(rigid_a, targetPosition_b, 4);
	auto bindingNodes_b = _bindingHandles(rigid_b, node_a.position, 4);

	// cross link them
	const float baseBindingStrength = 0.1f;

	FMatrix trans = FTranslationMatrix::Make(-node_b.position) // to identity
		* FRotationMatrix::Make(node_b.orientation.Inverse())  // ...
		* FRotationMatrix::Make(targetRotation_b)              // to target
		* FTranslationMatrix::Make(targetPosition_b);          // ...

	graph->beginTransaction();

	for (FGraphNodeHandle b_ : bindingNodes_b)
	{
		if (!b_) continue;

		const FVector& p_b = b_(*graph).position;

		FVector aligned_b = trans.TransformPosition(p_b);

		for (FGraphNodeHandle a_ : bindingNodes_a)
		{
			if (!a_) continue;

			const FVector& p_a = a_(*graph).position;

			float dist = (p_a - aligned_b).Size();

			auto edgeHandle = graph->connectNodes<FFlexConnection>(a_, b_, dist, baseBindingStrength);
		}
	}

	graph->endTransaction();
}

void UMLElementSimulation::_createBonds(float deltaT, std::vector< std::vector<Neighbourhood> >& ruleNeighbourhoods)
{
	auto& elements = graph->componentStorage<FMLElement>();

	for (auto& predictionNeighbourhoods : ruleNeighbourhoods)
	{
		for (Neighbourhood& neighbourhood : predictionNeighbourhoods)
		{
			FMLElement& element = elements.componentForNode(neighbourhood.element);

			FGraphNode node_base = graph->node(element.nodeHandle());
			FGraphNode rule_base = ruleGraph.node(neighbourhood.ruleElement);

			if (!neighbourhood.satisfied) continue;
			if (element.bound) continue;

			for (auto& pair : neighbourhood.pairs)
			{
				FGraphNode& node = graph->node(pair.first);
				FGraphNode& rule = ruleGraph.node(pair.second);

				FVector offset = rule_base.orientation.UnrotateVector(rule.position - rule_base.position);

				FVector target = node_base.position + node_base.orientation.RotateVector(offset);
				FQuat targetRotation = node_base.orientation * (rule_base.orientation.Inverse() * rule.orientation);

				_bind(node_base.handle(), node.handle(), target, targetRotation);
			}

			element.bound = true;
		}
	}

}

