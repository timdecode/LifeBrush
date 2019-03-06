// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Aggregates.h"

TArray<FGraphEdgeHandle> FMLAggregateNO::aggregateNodes(const TArray<FGraphNodeHandle>& nodes, FGraph& graph, FGraphNodeHandle hostNodeHandle)
{
	// add the FMLAggregateNO
	FGraphNode& hostNode = graph.node(hostNodeHandle);

	checkfSlow(hostNode.isValid(), TEXT("The hostNodeHandle must be valid"));
	checkfSlow(!hostNode.hasComponent<FMLAggregateNO>(), TEXT("The hostNodeHandle already has a FMLAggregateNO"));

	hostNode.addComponent<FMLAggregateNO>(graph);
	hostNode.addComponent<FMLAggregateNO_id>(graph, hostNodeHandle);


	// connect the aggregate
	TArray<FGraphEdgeHandle> edges;
	for (const FGraphNodeHandle& node : nodes)
	{
		checkfSlow(node != hostNodeHandle, TEXT("The hostNode must not be in nodes."));

		graph.node(node).addComponent<FMLAggregateNO_id>(graph, hostNodeHandle);
		auto newEdge = graph.connectNodes<FMLAggregateEO>(node, hostNodeHandle);

		edges.Add(newEdge);
	}

	return edges;
}

void FMLAggregateNO::aggregateEdges(const TArray<FGraphEdgeHandle>& edges, FGraph& graph, FGraphNodeHandle hostNodeHandle)
{
	// add the FMLAggregateNO
	FGraphNode& hostNode = graph.node(hostNodeHandle);

	checkfSlow(hostNode.isValid(), TEXT("The hostNodeHandle must be valid"));
	checkfSlow(!hostNode.hasComponent<FMLAggregateNO>(), TEXT("The hostNodeHandle already has a FMLAggregateNO"));

	hostNode.addComponent<FMLAggregateNO>(graph);

	// connect the aggregate
	auto aggregateStorage = graph.edgeStorage<FMLAggregateEO>();

	for (const FGraphEdgeHandle edgeHandle : edges)
	{
		checkfSlow(!aggregateStorage.contains(edgeHandle), TEXT("We have already aggregated this edge."));
		
		FGraphEdge& edge = graph.edge(edgeHandle);

		checkfSlow(edge.a == hostNode || edge.b == hostNode, TEXT("The edges must be connected to the hostNode."));

		graph.addOrReplaceEdgeObject<FMLAggregateEO>(edgeHandle);
	}
}

void FMLAggregateNO::removeAggregate(FGraph& graph, FGraphNodeHandle hostNodeHandle)
{
	FGraphNode& hostNode = graph.node(hostNodeHandle);

	checkfSlow(hostNode.isValid(), TEXT("The hostNodeHandle must be valid"));
	checkfSlow(hostNode.hasComponent<FMLAggregateNO>(), TEXT("The hostNodeHandle does not have a FMLAggregateNO"));

	graph.beginTransaction();

	for (auto ei : hostNode.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		auto other = graph.edge(edgeHandle).other(hostNodeHandle);

		graph.removeNode(other);
	}

	graph.removeNode(hostNodeHandle);

	graph.endTransaction();
}

FMLAggregateNO* FMLAggregateNO::getAggregate(FGraph& graph, FGraphNodeHandle subNodeHandle)
{
	auto& node = graph.node(subNodeHandle);

	if (node.hasComponent<FMLAggregateNO>())
		return &node.component<FMLAggregateNO>(graph);

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if (auto aggregateEdge = graph.edgeObjectPtr<FMLAggregateEO>(edgeHandle))
		{
			FGraphEdge& edge = graph.edge(edgeHandle);

			return &edge.other(subNodeHandle).node(graph).component<FMLAggregateNO>(graph);
		}
	}

	return nullptr;
}

void FMLAggregateNO::edgesAndNodesInAggregate(FGraph& graph, TArray<FGraphNodeHandle>& nodes_out, TArray<FGraphEdgeHandle>& edges_out)
{
	auto aggregateStorage = graph.edgeStorage<FMLAggregateEO>();

	FGraphNode& node = graph.node(nodeHandle());

	nodes_out.Add(nodeHandle());

	for (auto ei : node.edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if( !aggregateStorage.objectPtr(edgeHandle) ) continue;

		edges_out.Add(edgeHandle);
	
		FGraphEdge& edge = graph.edge(edgeHandle);

		nodes_out.Add(edge.other(nodeHandle()));
	}
}

// ----------------------------------------------
// UMLAggregateSimulation
// ----------------------------------------------

void UMLAggregateSimulation::attach()
{
	graph->addComponentListener<FMLAggregateNO>(this);

	// update the FMLAggregateNO_id for each node in the aggregate
	graph->each_node_object<FMLAggregateNO>([&](FGraphNode& node, FMLAggregateNO& aggregateNode) {
		_updateIDs(node);
	});
}

void UMLAggregateSimulation::detach()
{
	graph->removeComponentListener<FMLAggregateNO>(this);
}

void UMLAggregateSimulation::_updateIDs(FGraphNode& node)
{
	node.addComponent<FMLAggregateNO_id>(*graph, node.handle());

	node.each<FMLAggregateEO>(*graph, [&](FGraphNodeHandle other, FMLAggregateEO& edgeObject) {
		FGraphNode& otherNode = graph->node(other);

		otherNode.addComponent<FMLAggregateNO_id>(*graph, node.handle());
	});
}

void UMLAggregateSimulation::componentAdded(FGraphNodeHandle handle, ComponentType type)
{
	static ComponentType aggregateNodeType = componentType<FMLAggregateNO>();

	if (type == aggregateNodeType)
	{
		_updateIDs(graph->node(handle));
	}
}


