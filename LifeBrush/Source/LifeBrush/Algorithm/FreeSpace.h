// Copyright 2016, Timothy Davison. All rights reserved.

#pragma once

#include "Eigen/Dense"
#include <vector>
#include "Algorithm/SurfaceIndex.h"
#include "Element.h"

struct FreeSpaceCluster
{
public:
	// the centroid of the element positions in this cluster
	Eigen::Vector3f centroid = Eigen::Vector3f::Zero();

	// each column is a basis vector of the distribution contained in this cluster
	Eigen::Matrix3f basisVectors = Eigen::Matrix3f::Zero();

	float radius = 0.0f;

	// the elements that belong to this cluster
	// the offset vectors from the centroid of the elements
	// these vectors come from overlapping the k-similar neighbourhoods with the source of the cluster
	std::vector < std::pair < Eigen::Vector3f, FGraphNodeHandle> > offsetsAndElements;
};

struct OutputFreeSpaceCluster
{
	FreeSpaceCluster * exampleCluster = nullptr;

	Eigen::Vector3f position = Eigen::Vector3f::Zero();

	FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface;
};

//template<typename ElementType>
//class LIFEBRUSH_API SpatialIndex
//{
//public:
//	ElementType& emplace_back();
//	ElementType& back();
//
//	void remove( size_t index );
//	ElementType& at( size_t index );
//
//	bool overlap( Eigen::Vector3f point, float radius );
//
//protected:
//	std::vector<size_t> _recycledIndices;
//};



