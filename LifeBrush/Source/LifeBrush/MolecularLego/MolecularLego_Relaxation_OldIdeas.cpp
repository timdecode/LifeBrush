// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "MolecularLego_Relaxation.h"

#include "Simulation/FlexElements.h"




void UMLElementSimulation::_repulsion(
	std::vector<Neighbourhood> predictionNeighbourhoods,
	std::map<FGraphNodeHandle, unsigned int>& problem,
	std::vector<FGraphNodeHandle>& indexToElement,
	std::vector<std::vector<PredictedPosition>>& predictions_in_out)
{
	auto& elements = graph->componentStorage<FMLElement>();

	auto& rules = ruleGraph.componentStorage<FMLElementRule>();

	std::vector<Neighbourhood> predictions;

	const float minDistSqrd = std::pow(minAssignmentDistance, 2.0f);

	for (FMLElement& e_a : elements)
	{
		FGraphNode& node_a = graph->node(e_a.nodeHandle());

		FQuat q_a = node_a.orientation;
		FVector p_a = node_a.position;

		auto found_i = problem.find(e_a.nodeHandle());

		size_t i = 0;

		if (found_i == problem.end())
		{
			i = indexToElement.size();
			problem[e_a.nodeHandle()] = i;
			indexToElement.push_back(e_a.nodeHandle());
			predictions_in_out.emplace_back();
		}
		else
			i = found_i->second;

		// find the neighborhood 
		unrealAABB::AABB query(p_a, minAssignmentDistance);

		std::vector<PositionQuat> aPositions;
		std::vector<FGraphNodeHandle> aHandles;


		FQuat inverseRotation = q_a.Inverse();

		elementBVH.query(query, [&](unsigned int particleIndex) {
			FGraphNodeHandle h_b(particleIndex);

			FMLElement& e_b = elements.componentForNode(h_b);

			if (&e_a == &e_b) return true;

			auto found_j = problem.find(h_b);

			size_t j = 0;

			if (found_j == problem.end())
			{
				j = indexToElement.size();
				problem[e_b.nodeHandle()] = j;
				indexToElement.push_back(e_b.nodeHandle());
				predictions_in_out.emplace_back();
			}
			else
				j = found_j->second;

			FGraphNode& node_b = graph->node(h_b);

			FVector dir = node_b.position - node_a.position;


			if (dir.SizeSquared() > minDistSqrd)
				return true;

			const float distSqrd = dir.Size();

			dir = dir.GetSafeNormal();

			FVector p = node_b.position + dir * (minAssignmentDistance - std::sqrt(distSqrd));

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			predictions_in_out[j].emplace_back(
				i,
				eigen(p),
				node_b.orientation
			);

			return true; // keep querying
		});
	}
}

