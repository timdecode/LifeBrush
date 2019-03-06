//
//  Domain.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-06-25.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "Domain.h"

void Domain::componentAdded(FGraphNodeHandle nodeHandle, ComponentType type)
{
	if (componentType<FElementObject>() != type)
		return;

	FGraphNode& node = nodeHandle(graph);

	_kdTree.add(node);
}

void Domain::componentRemoved(FGraphNodeHandle nodeHandle, ComponentType type)
{
	if (componentType<FElementObject>() != type)
		return;

	_kdTree.remove(nodeHandle);
}
