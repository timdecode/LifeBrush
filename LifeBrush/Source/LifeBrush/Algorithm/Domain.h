//
//  Domain.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-06-25.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include <numeric> 
#include <algorithm>
#include <iterator>

#include "Element.h"
#include "ElementKDTree.h"
#include "Entity.hpp"
#include "Utility.h"

#include "ShipEditorSimulation/GraphCopyContext.h"


struct Domain : public ComponentListener
{
protected:

	ElementKDTree _kdTree;

public:
	Domain(FGraph& graph_in) : graph(graph_in), _kdTree(graph_in)
	{
		//graph.addListener<FElementObject>(this);
	}

	virtual ~Domain() {}

	FGraph& graph;

	// For a node created in the graph directly, this notifies the domain of the new node and adds it to the internal
	// data structures. Don't call this twice with the same node. Don't call this with a node created using Domain::insert.
	void didInsert(FGraphNodeHandle handle)
	{
		_kdTree.add(handle(graph));
	}

	auto insert(FVector newPosition, FQuat orientation = FQuat::Identity, float scale = 1.0f) -> FGraphNodeHandle
	{
		FGraphNodeHandle newNodeHandle(graph.addNode(newPosition, orientation, scale));
		 
		newNodeHandle(graph).addComponent<FElementObject>(graph);

		_kdTree.add(newNodeHandle(graph));

		return newNodeHandle;
	}

	auto insert(FElementObject& elementObject, FGraph& sourceGraph, FVector newPosition) -> FGraphNodeHandle
	{
		FGraphNode& elementNode = elementObject.node(sourceGraph);

		FGraphNodeHandle newNodeHandle(graph.addNode(newPosition, elementNode.orientation, elementNode.scale));
		FGraphNode& newNode = newNodeHandle(graph);

		for (auto t : elementNode.components)
		{
			FGraphCopyContext::copyComponent(t, elementNode, newNode, sourceGraph, graph);
		}

		_kdTree.add(newNode);


		return newNodeHandle;
	}

	// copies the elements into the domain
	auto insert(FElementObject& elementObject, FGraph& sourceGraph) -> FGraphNodeHandle
	{
		FGraphNode& elementNode = elementObject.node(sourceGraph);

		FGraphNodeHandle newNodeHandle(graph.addNode(elementNode.position, elementNode.orientation, elementNode.scale));
		FGraphNode& newNode = newNodeHandle(graph);

		for (auto t : elementNode.components)
		{
			FGraphCopyContext::copyComponent(t, elementNode, newNode, sourceGraph, graph);
		}

		_kdTree.add(newNode);


		return newNodeHandle;
	}

	// copies the template element into the domain
	auto insert(FGraphNode& elementNode, FGraph& sourceGraph) -> FGraphNodeHandle
	{
		if (!elementNode.hasComponent<FElementObject>())
			return FGraphNodeHandle();

		FGraphNodeHandle newNodeHandle(graph.addNode(elementNode.position, elementNode.orientation, elementNode.scale));
		FGraphNode& newNode = newNodeHandle(graph);

		for (auto t : elementNode.components)
		{
			FGraphCopyContext::copyComponent(t, elementNode, newNode, sourceGraph, graph);
		}

		_kdTree.add(newNode);

		return newNodeHandle;
	}


    void erase( const std::vector<FGraphNodeHandle>& toRemove )
    {
		if (toRemove.size() == 0)
			return;

		graph.beginTransaction();

		for (auto e : toRemove)
		{
			graph.removeNode(e.index);
			_kdTree.remove(e);
		}

		graph.endTransaction();

		rebalance();
	}

	void emptyRecycleBin()
	{
	}
    
    void setElements(const std::vector<FElementObject>& elements, FGraph& sourceGraph, bool createEntities = false)
    {
        clear();

		//graph.removeListener<FElementObject>(this);

		for (auto e : elements)
		{
			FGraphNode& node = sourceGraph.node(e.nodeIndex);
			FGraphNodeHandle newNodeHandle(graph.addNode(node.position, node.orientation, node.scale));

			for (auto t : node.components)
			{
				FGraphCopyContext::copyComponent(t, node, newNodeHandle(graph), sourceGraph, graph);
			}
		}

		//graph.addListener<FElementObject>(this);

		rebalance();
	}
    
    size_t size()
    {
        return graph.numNodes();
    }

    
	FElementObject& elementAt(const unsigned int index)
    {
        return graph.node(index).component<FElementObject>(graph);
    }
    
    void clear()
    {
		graph.clear();

        _kdTree.clear();
    }

	void componentAdded(FGraphNodeHandle nodeHandle, ComponentType type);

	void componentRemoved(FGraphNodeHandle nodeHandle, ComponentType type);

    // queries
    ElementKDTree::NearestResult nearest(const Eigen::Vector3f& point)
    {
        return _kdTree.nearest(point);
    }
    
    // queries
    ElementKDTree::NearestResult nearest(const Eigen::Vector3f point, float radius)
    {
        auto result = _kdTree.nearest(point);

		if(result.distanceSquared < radius * radius)
			return result;
		else
			return { FGraphNodeHandle(), -1.0f };
    }
    
    std::vector<FGraphNodeHandle> nearestKNeighbours(const Eigen::Vector3f& point, size_t count, std::vector<float>& distances_out)
    {
        return _kdTree.nearestKNeighbours(point, count, distances_out);
    }
    
    std::vector<FGraphNodeHandle> nearestKNeighbours(const Eigen::Vector3f& point, size_t count)
    {
        return _kdTree.nearestKNeighbours(point, count);
    }
    
    std::vector<FGraphNodeHandle> nearestInRadius(const Eigen::Vector3f& point, float radius, std::vector<float>& distances_out, int count = -1, bool sort = true)
    {
        return _kdTree.nearestInRadius(point, radius, distances_out, count, sort);
    }
    
    std::vector<FGraphNodeHandle> nearestInRadius(const Eigen::Vector3f& point, float radius, int count = -1, bool sort = true)
    {
        return _kdTree.nearestInRadius(point, radius, count, sort);
    }

	// Hacked, we should really add a nearestInBox implementation to the k-d tree. This would be relatively trivial given that it
	// is a k-d tree.
	std::vector<FGraphNodeHandle> nearestInBox( const Eigen::Vector3f& point, const Eigen::Vector3f& halfExtents, int count = -1, bool sort = true )
	{
		// c^2 = a^2 + b^2
		// so, a = b: c^2 = 2a^2
		// hence, we have our radius
		float maxHalfExtents = std::max( halfExtents( 0 ), std::max( halfExtents( 1 ), halfExtents( 2 ) ) );
		float radius = sqrt( 2 ) * maxHalfExtents;

		Eigen::AlignedBox3f aabb = Eigen::AlignedBox3f( point - halfExtents, point + halfExtents );

		auto inRadius = _kdTree.nearestInRadius( point, radius, count, sort );

		std::vector<FGraphNodeHandle> inBox;

		for(FGraphNodeHandle e : inRadius)
		{
			auto position = e(graph).position;

			if( aabb.contains(eigen(position)) )
				inBox.emplace_back( e );
		}

		return inBox;
	}
    
    std::vector<FGraphNodeHandle> neighbours(FGraphNodeHandle& element, float radius, int count = -1, bool sort = true, std::vector<float>* distances_out = nullptr, bool surfaceVolumeInteraction = true)
    {
        return _kdTree.neighbours(element, radius, count, sort, distances_out, surfaceVolumeInteraction );
    }
    
    // range for compatability
    // let's switch to this idea:
    // http://stackoverflow.com/a/352162
    auto begin() -> FElementObject*
    {
		return graph.componentStorage<FElementObject>().begin();
    }
    
    auto end() -> FElementObject*
    {
		return graph.componentStorage<FElementObject>().end();
	}

	std::vector<FGraphNodeHandle> elements()
	{
		std::vector<FGraphNodeHandle> elements;
		
		auto& storage = graph.componentStorage<FElementObject>();
		
		for (auto& e : storage)
		{
			if( e.isValid() )
				elements.push_back(e.nodeHandle());
		}

		return elements;
	}
    
    void rebalance()
    {
		auto toSet = elements();
		_kdTree.setElements(toSet);
    }

	Eigen::AlignedBox3f computeAABB()
	{
		auto& storage = graph.componentStorage<FElementObject>();

		if(storage.size() == 0 )
			return Eigen::AlignedBox3f( Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero() );

		const auto firstPostion = storage[0].node(graph).position;

		Eigen::AlignedBox3f b(eigen(firstPostion), eigen(firstPostion));

		for (auto& e : storage)
		{
			const auto p = e.node(graph).position;

			b.extend(eigen(p));
		}

		return b;
	}
        

};
