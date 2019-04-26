// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include <vector>

#include "ShipEditorSimulation/GraphSnapshot.h"

void FGraphSnapshot::snapshot(FGraph& graph_in)
{
	graph.init();

	graph.clear();


	tickCount = graph_in.tickCount;

	// copy all nodes, even ones marked for recycling
	graph.allNodes = graph_in.allNodes;

	graph.privateEdges = graph_in.privateEdges;

	graph._componentStorage.Empty();
	for (FComponentStorage& storage : graph_in._componentStorage)
	{
		size_t n = storage.size();
		UScriptStruct * componentClass = storage.componentClass;

		graph._componentStorage.Add(storage);
		
		FGraphObject * objects = storage.at(0, componentClass);
		FGraphObject * newObjects = storage.at(0, componentClass);

		componentClass->CopyScriptStruct(newObjects, objects, n);
	}
}

void FGraphSnapshot::restore(FGraph& graph_in)
{
	graph_in.clear();

	graph_in.beginTransaction();

	graph_in.tickCount = tickCount;

	// load nodes
	std::vector<int32> oldToNewNodes(graph.allNodes.Num());

	for (FGraphNode& node : graph.allNodes)
	{
		if (!node.isValid())
		{
			oldToNewNodes[node.id] = -1;
			continue;
		}

		FGraphNodeHandle newNode = FGraphNodeHandle(graph_in.addNode(node.position, node.orientation, node.scale));

		oldToNewNodes[node.id] = newNode.index;
	}

	// load components
	for (FComponentStorage& storage : graph._componentStorage)
	{
		size_t n = storage.size();

		UScriptStruct * componentClass = storage.componentClass;
		ComponentType componentType = FGraphObject::componentType(componentClass);

		FComponentStorage& graphStorage = graph_in.componentStorage(componentType);

		if (componentClass->StructFlags & STRUCT_IsPlainOldData)
		{
			// directly copy the data over
			graphStorage = storage;
		}
		else
		{
			int32 structSize = componentClass->GetStructureSize();

			// reserve size
			graphStorage = storage;

			if (n > 0)
			{
				void * source = storage.pool.begin(structSize);
				void * dest = graphStorage.pool.begin(structSize);

				componentClass->CopyScriptStruct(dest, source, n);
			}
		}

		for (int i = 0; i < n; ++i)
		{
			FGraphObject * object = graphStorage.at(i, componentClass);

			object->nodeIndex = oldToNewNodes[object->nodeIndex];

			FGraphNode& node = graph_in.node(object->nodeIndex);
			node.components.Add(componentType);

			graph_in.componentAdded(node, componentType);
		}
	}

	// load edges
	for (FGraphEdge& edge : graph.privateEdges)
	{
		if (!edge.isValid())
			continue;

		// TODO
		UE_LOG(LogTemp, Warning, TEXT("FGraphSnapShot - connection restoration is not implemented"));


		//graph_in.connectNodes(
		//	FGraphNodeHandle(oldToNewNodes[edge.a]),
		//	FGraphNodeHandle(oldToNewNodes[edge.b]),
		//	edge.type, 
		//	edge.halfExtents
		//);
	}

	graph_in.endTransaction();
}