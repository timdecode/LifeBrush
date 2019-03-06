//
//  Utility.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-07-20.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ShipEditorSimulation/Graph.h"
#include "TimStructBox.h"

#include "Utility.h"


TArray<FTimStructBox> Utility::boxComponents(struct FGraphNode& node, struct FGraph& graph)
{
	TArray<FTimStructBox> result;
	result.Reserve(node.components.Num());

	for (auto t : node.components)
	{
		result.Emplace();

		FTimStructBox& box = result.Last();

		box.scriptStruct = FGraphObject::componentStruct(t);
		box.initStruct();

		FGraphObject * sourceComponent = graph.component(node.handle(),t);
		box.scriptStruct->CopyScriptStruct(box.structMemory, sourceComponent);
	}

	return result;
}


void Utility::unboxComponents(TArray<FTimStructBox>& components, struct FGraphNode& destination, struct FGraph& graph)
{
	for (FTimStructBox& box : components)
	{
		if (!box.scriptStruct)
			continue;

		auto componentType = FGraphObject::componentType(box.scriptStruct);

		FGraphObject * destinationComponent = destination.addComponent(graph, componentType);

		if(box.structMemory)
			box.scriptStruct->CopyScriptStruct(destinationComponent, box.structMemory);
	}
}


