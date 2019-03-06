//
//  ElementKDTree.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-12-01.
//  Copyright © 2015 Epic Games, Inc. All rights reserved.
//

#pragma once

#include <vector>
#include <iostream>
#include <memory>
#include <unordered_set>

#include "Element.h"
#include "PointCloud.hpp"
#include "Utility.h"



struct ElementKDTree
{
	typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
		nanoflann::L2_Simple_Adaptor<float, PointTypedCloud<float, Eigen::Vector3f> >,
		PointTypedCloud<float, Eigen::Vector3f>,
		3 /* dim */
	> IndexAdaptor;

public:
    FGraph& graph;

    ElementKDTree(FGraph& graph) : graph(graph)
    {
        _setupIndex();
    }
    
    
    // Add the element without claiming ownership over the pointer
    void add(FGraphNode& node)
    {
        size_t newIndex = _cloud.kdtree_get_point_count();

        auto p = node.position;

        _cloud.pts.emplace_back(eigen(p));
        
        _nodeHandles.push_back(node.handle());

		_indexOfHandle[node.handle()] = newIndex;

		_index->addPoints( newIndex, newIndex );
    }

	void remove(FGraphNodeHandle handle)
	{
		size_t index = _indexOfHandle[handle];

		_index->removePoint(index);

		_nodeHandles[index] = FGraphNodeHandle::null;
		_indexOfHandle.erase(handle);
	}
    
    size_t size()
    {
		// we use _indexOfHandle to track size, as we keep deleted nodes around in _nodeHandles (as nulls)
        return _indexOfHandle.size();
    }
    
    template<typename Container>
    void setElements(Container& elements)
    {
        clear();

		_setupIndex();

		Element::eachNode(elements, graph, [&](FGraphNode& node) {
			_cloud.pts.emplace_back(eigen(node.position));

			_indexOfHandle[node.handle()] = _nodeHandles.size();

			_nodeHandles.push_back(node.handle());
		});
        
		if( _cloud.pts.size() > 0 )
			_index->addPoints( 0, _cloud.pts.size() - 1 );
    }

    void clear()
    {
        _cloud.pts.clear();
        _nodeHandles.clear();
		_indexOfHandle.clear();
        
        _setupIndex();
    }
    
    struct NearestResult
    {
        FGraphNodeHandle element;
        float distanceSquared;
    };
    
    // queries
    NearestResult nearest(const Eigen::Vector3f point)
    {
        if( _nodeHandles.size() == 0 )
            return NearestResult();
        
        bool found = false;
        
        float distanceSquared;
        size_t index;
        
        _index->knnSearch(&point(0), 1, &index, &distanceSquared);
        
        return {_nodeHandles[index], distanceSquared};
    }
    
    std::vector<FGraphNodeHandle> nearestKNeighbours(const Eigen::Vector3f point, size_t count, std::vector<float>& distances_out)
    {
		std::vector<FGraphNodeHandle> result;

		if (count == 0)
			return result;
        
        if( count > _nodeHandles.size() )
            count = _nodeHandles.size();
        
        
        if( _nodeHandles.size() == 0 )
            return result;
        
        result.reserve(count);
        
        std::vector<size_t> indices(count);
        std::vector<float> distances(count);
        
        nanoflann::KNNResultSet<float,size_t> resultSet(count);
        resultSet.init(&indices[0], &distances[0]);
        _index->findNeighbors(resultSet, &point(0), nanoflann::SearchParams());
        
        for( int i = 0; i < resultSet.size(); ++i )
        {
            FGraphNodeHandle handle = _nodeHandles[indices[i]];
            
            result.emplace_back(handle);
            distances_out.push_back(distances[i]);
        }
        
        return result;
    }
    
	std::vector<FGraphNodeHandle> nearestKNeighbours(const Eigen::Vector3f point, size_t count)
    {
		std::vector<FGraphNodeHandle> result;

        if( count == 0 )
            return result;
        
        if( _nodeHandles.size() == 0 )
            return result;
        
        if( count > _nodeHandles.size() )
            count = _nodeHandles.size();
        
        std::vector<size_t> indices(count);
        std::vector<float> distances(count);
        
        nanoflann::KNNResultSet<float,size_t> resultSet(count);
        resultSet.init(&indices[0], &distances[0]);
        _index->findNeighbors(resultSet, &point(0), nanoflann::SearchParams());
        
		result.resize(resultSet.size());

        for( int i = 0; i < resultSet.size(); ++i )
        {
			result[i] = _nodeHandles[indices[i]];
        }
        
        return result;
    }
    
    std::vector<FGraphNodeHandle> nearestInRadius(const Eigen::Vector3f point, float radius, std::vector<float>& distances_out, int count = -1, bool sort = true)
    {
		std::vector<FGraphNodeHandle> result;

        if( _nodeHandles.size() == 0 )
			return result;

        if( radius <= 0.0f )
			return result;

        std::vector<std::pair<size_t, float> > dynoResult;
        
        nanoflann::SearchParams searchParams;
        
        searchParams.sorted = sort;
        
        _index->radiusSearch(&point(0), radius * radius, dynoResult, searchParams);

		result.resize(dynoResult.size());

        size_t i = 0;
        for( auto pair : dynoResult )
        {
			result[i] = _nodeHandles[pair.first];

            ++i;

            distances_out.push_back(pair.second);
        }
        
        return result;
    }
    
	std::vector<FGraphNodeHandle> nearestInRadius(const Eigen::Vector3f point, float radius, int count = -1, bool sort = true)
    {
		std::vector<FGraphNodeHandle> result;

		if (_nodeHandles.size() == 0)
			return result;

		if (radius <= 0.0f)
			return result;
        
        std::vector<std::pair<size_t, float> > dynoResult;
        
        nanoflann::SearchParams searchParams;
        
        searchParams.sorted = sort;
        
        _index->radiusSearch(&point(0), radius * radius, dynoResult, searchParams);
        
		result.resize(dynoResult.size());

        size_t i = 0;
        for( auto& pair : dynoResult )
        {
			result[i] = _nodeHandles[pair.first];
            
            ++i;
        }
        
        return result;
    }
    
    std::vector<FGraphNodeHandle> neighbours(FGraphNodeHandle& elementHandle, float radius, int count = -1, bool sort = true, std::vector<float>* distances_out = nullptr, bool surfaceVolumeInteraction  = true)
    {
		std::vector<FGraphNodeHandle> result;

		if (_nodeHandles.size() == 0)
			return result;

		if (radius <= 0.0f)
			return result;

		ElementTuple element(elementHandle, graph);
        
        std::vector<std::pair<size_t, float> > dynoResult;
        
        nanoflann::SearchParams searchParams;
        
        searchParams.sorted = sort;
        

		FVector& p = element.node.position;
        _index->radiusSearch(&p[0], radius * radius, dynoResult, searchParams);
        
        if( distances_out != nullptr )
            distances_out->clear();
        
		if(dynoResult.size() == 0)
			return result;

		result.reserve(dynoResult.size() - 1);
        
        for( auto& pair : dynoResult )
        {
            FGraphNodeHandle handle = _nodeHandles.at(pair.first);
			FElementObject& e = handle(graph).component<FElementObject>(graph);
            
            if( handle == elementHandle)
                continue;

			// don't allow surface elements to interact with volume elements
			if(!surfaceVolumeInteraction && element.element.surfaceIndex.isOnSurface() && !e.surfaceIndex.isOnSurface() )
				continue;

            result.push_back(handle);
            
            if( distances_out != nullptr )
                distances_out->push_back(pair.second);
        }
        
        return result;
    }
    
    void rebalance()
    {
        _updatePositions();

	
	}
    
    auto begin() -> decltype( std::vector<FGraphNodeHandle>().begin() )
    {
        return _nodeHandles.begin();
    }
    
    auto end() -> decltype( std::vector<FGraphNodeHandle>().begin() )
    {
        return _nodeHandles.end();
    }
    
private:
    void _setupIndex()
    {
		_index = std::make_unique<IndexAdaptor>(3, _cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
    }
    
    void _updatePositions()
    {
		std::vector<Eigen::Vector3f> points = std::move( _cloud.pts );
		_cloud.pts = std::vector<Eigen::Vector3f>();

		_setupIndex();

        for( int i = 0; i < _nodeHandles.size(); ++i )
        {
            auto& nodeHandle = _nodeHandles[i];
            const auto position = nodeHandle(graph).position;
            
			points[i] = { position[0], position[1], position[2] };
        }

		_cloud.pts = points;

		if( _cloud.pts.size() > 0 )
			_index->addPoints( 0, _cloud.pts.size() - 1 );
    }


    
	PointTypedCloud<float, Eigen::Vector3f> _cloud;
    std::unique_ptr<IndexAdaptor> _index = nullptr;

    std::vector<FGraphNodeHandle> _nodeHandles;
	std::unordered_map<FGraphNodeHandle, size_t> _indexOfHandle;
};