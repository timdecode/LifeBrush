// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "Graph.h"
#include "MirrorPoints.h"

/**
 * 
 */
struct LIFEBRUSH_API FGraphCopyContext
{
public:
	void clear();

public:
	// The duplicated nodes in the target graph.
	TArray<FGraphNodeHandle> duplicatedNodes; 

public:
	static FGraphCopyContext copySubgraph(
		FGraph& sourceGraph,
		FGraph& targetGraph,
		TArray<FGraphNodeHandle>& nodes,
		TArray<FGraphEdgeHandle>& edges,
		const FTransform localTransform);


	static FGraphNodeHandle copyAggregate(
		FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph,
		FGraph& targetGraph,
		const FVector position,
		const FQuat rotation_in
	);

	static void copyComponent(ComponentType componentType, FGraphNode& sourceNode, FGraphNode& targetNode, FGraph& sourceGraph, FGraph& targetGraph);


};
