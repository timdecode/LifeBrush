// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"

#include "ShipEditorSimulation/MeshSimulation.h"
#include "GraphCopyContext.h"

void FGraphCopyContext::clear()
{
	duplicatedNodes.Empty();
}

FGraphCopyContext FGraphCopyContext::copySubgraph(
	FGraph& sourceGraph,
	FGraph& targetGraph,
	TArray<FGraphNodeHandle>& nodes,
	TArray<FGraphEdgeHandle>& edges,
	const FTransform localTransform)
{
	FGraphCopyContext context;

	std::unordered_map<int32, int32> toTargetNodes;
	std::unordered_map<FGraphEdgeHandle, FGraphEdgeHandle> toTargetEdges;

	float scale = FMath::Abs( localTransform.GetScale3D().GetMin() );

	targetGraph.beginTransaction();

	// copy the edges first (we don't care about edge.a or edge.b yet)
	for(FGraphEdgeHandle ei : edges)
	{
		// copy the edge, we don't want a reference inside of the array, because we could copy back into the array
		// and Unreal hates that (asserts on targetGraph.privateEdges.Add(edge) if targetGraph == sourceGraph. 
		FGraphEdge edge = sourceGraph.edge(ei);

		if(!edge.isValid())
			continue;

		auto targetIndex = FGraphEdgeHandle(targetGraph.privateEdges.Add( edge ));
		FGraphEdge& targetEdge = targetGraph.privateEdges[targetIndex.index];

		targetEdge.halfExtents *= scale;

		toTargetEdges[ei] = targetIndex;
	}

	// copy the graph edge objects
	for (FEdgeStorage& sourceStorage : sourceGraph._edgeStorage)
	{
		UScriptStruct * scriptStruct = sourceStorage.scriptStruct();

		auto type = FGraphEdgeObject::edgeType(scriptStruct);

		FEdgeStorage& targetStorage = targetGraph.rawEdgeStorage(type);

		for (FGraphEdgeHandle sourceHandle : edges)
		{
			FGraphEdgeHandle targetHandle = toTargetEdges[sourceHandle];

			FGraphEdgeObject * sourceObject = sourceStorage.at(sourceHandle);

			if( sourceObject )
				targetStorage.add(targetHandle, *sourceObject);
		}
	}

	// copy the nodes
	FTransform unflippedTransform = localTransform;
	unflippedTransform.SetScale3D( unflippedTransform.GetScale3D().GetAbs() );


	FVector centroid = FVector::ZeroVector;

	for(FGraphNodeHandle ni : nodes)
	{
		// again, don't use references here, because we might end up copying back into the same array
		FGraphNode sourceNode = sourceGraph.node(ni);

		if(!sourceNode.isValid())
			continue;

		int32 targetIndex = targetGraph.allNodes.Add( sourceNode );
		FGraphNode& targetNode = targetGraph.allNodes[targetIndex];

		targetNode.position = localTransform.TransformPosition( sourceNode.position );
		targetNode.orientation = localTransform.TransformRotation(sourceNode.orientation);
		targetNode.scale = scale * sourceNode.scale;

		targetNode.id = targetIndex;

		// this will be populated in the next step
		targetNode.edges.Empty();

		// copy the components
		// this will become a method on the components themselves
		targetNode.components.Empty();

		for(ComponentType type : sourceNode.components)
		{
			copyComponent(type, sourceNode, targetNode, sourceGraph, targetGraph);
		}

		// update the centroid
		centroid += targetNode.position;

		// update the mapping
		toTargetNodes[sourceNode.id] = targetIndex;

		context.duplicatedNodes.Add( FGraphNodeHandle(targetIndex) );
	}

	centroid /= float( sourceGraph.allNodes.Num() );

	// link the edges and nodes up
	for(auto pair : toTargetEdges)
	{
		FGraphEdge& edge = sourceGraph.privateEdges[pair.first.index];
		FGraphEdge& targetEdge = targetGraph.privateEdges[pair.second.index];

		targetEdge.a = toTargetNodes[edge.a];
		targetEdge.b = toTargetNodes[edge.b];

		FGraphNode& targetA = targetGraph.allNodes[targetEdge.a];
		FGraphNode& targetB = targetGraph.allNodes[targetEdge.b];

		targetA.edges.Add( pair.second.index );
		targetB.edges.Add( pair.second.index );
	}

	targetGraph.endTransaction();

	return context;
}

void FGraphCopyContext::copyComponent(ComponentType componentType, FGraphNode& sourceNode, FGraphNode& targetNode, FGraph& sourceGraph, FGraph& targetGraph)
{
	FGraphObject * component = sourceNode.component(sourceGraph, componentType);
	FGraphObject * targetComponent = targetNode.addComponent(targetGraph, componentType);

	UScriptStruct * scriptStruct = FGraphObject::componentStruct(componentType);

	// it's very important that the C++ constructor gets invoked
	scriptStruct->GetCppStructOps()->Construct(targetComponent);

	scriptStruct->CopyScriptStruct(targetComponent, component);

	targetComponent->nodeIndex = targetNode.id;
}
