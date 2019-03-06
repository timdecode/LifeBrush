// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/ObjectSimulation.h"

#include "Aggregates.generated.h"

// An aggregate is a meta-agent composed of FGraphNdoes connected together by FMLAggregateEO edges. Each each
// connects to a central node, with a FMLAggregateNO. Removing the FMLAggregateNO removes the FMLAggregateEO edges.
// To create an aggregate, use FMLAggregateNO::aggregateNodes or FMLAggregateNO::aggregateEdges.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLAggregateEO : public FGraphEdgeObject
{
	GENERATED_BODY()

public:

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLAggregateNO_id : public FGraphObject
{
	GENERATED_BODY()

public:
	FMLAggregateNO_id() {}
	FMLAggregateNO_id(FGraphNodeHandle aggregateNode) : aggregateNode(aggregateNode) {}

	UPROPERTY()
	FGraphNodeHandle aggregateNode = FGraphNodeHandle::null;
};

// An aggregate is a meta-agent composed of FGraphNdoes connected together by FMLAggregateEO edges. Each each
// connects to a central node, with a FMLAggregateNO. Removing the FMLAggregateNO removes the FMLAggregateEO edges.
//
// An aggregate can be removed as one with FMLAggregateNO::removeAggregate.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLAggregateNO : public FGraphObject
{
	GENERATED_BODY()

public:
	// Aggregates the nodes together with UMLAggregateEO connections. All connections will be connected to the hostNode. 
	// @param nodes Must not contain the hostNode. If it does, I'll crash you.
	// @param graph The graph containing the node handles.
	// @param hostNode All connections in the aggregate will stem from the hostNode.
	static TArray<FGraphEdgeHandle> aggregateNodes(const TArray<FGraphNodeHandle>& nodes, FGraph& graph, FGraphNodeHandle hostNode);

	// Adds a FMLAggregateEO to each edge an a FMLAggregateNO to the host node. Each each must already be connected 
	// to the hostNode. If not, I'll crash you.
	// @param edges The edge handles to which we will add FMLAggregateEO-s. Each edge must be connected to the hostNode already.
	// @param graph The graph containing the node handles.
	// @param hostNode All connections in the aggregate will stem from the hostNode. We'll add a FMLAggregateNO to the host node.
	static void aggregateEdges(const TArray<FGraphEdgeHandle>& edges, FGraph& graph, FGraphNodeHandle hostNode);

	// Assuming the hostNode has a FMLAggregateNO, this function will remove all of the nodes connected to it with a FMLAggregateEO.
	static void removeAggregate(FGraph& graph, FGraphNodeHandle hostNode);

	// Returns the aggregate node object for the provided subNodeHandle, if it is part of an aggregate.
	static FMLAggregateNO* getAggregate(FGraph& graph, FGraphNodeHandle subNodeHandle);

	// Query to find the nodes and edges connected to this aggregate.
	// @param nodes_out The nodes connected to the aggregate, including this node. This node will the best first entry in nodes_out.
	// @param edges_out The edges connecting the aggregate.
	void edgesAndNodesInAggregate(FGraph& graph, TArray<FGraphNodeHandle>& nodes_out, TArray<FGraphEdgeHandle>& edges_out);
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UMLAggregateSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	virtual void attach() override;
	virtual void detach() override;

	virtual void componentAdded(FGraphNodeHandle node, ComponentType type) override;

protected:
	void _updateIDs(FGraphNode& node);
};
