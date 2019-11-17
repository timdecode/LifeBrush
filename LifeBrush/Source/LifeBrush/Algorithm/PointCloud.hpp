//
//  PointCloud.hpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-09-08.
//  Copyright © 2015 EpicGames. All rights reserved.
//

#pragma once

#include "Eigen/Dense"
#include "Algorithm/nanoflann.hpp"

struct EigenPointCloud
{
    std::vector<Eigen::Vector3f>  pts;
    
    // Must return the number of data points
    inline size_t kdtree_get_point_count() const { return pts.size(); }
    
    // Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
    inline float kdtree_distance(const float *p1, const size_t idx_p2,size_t /*size*/) const
    {
        const auto& p2_3f = pts[idx_p2];
        const auto& p1_3f = *reinterpret_cast<const Eigen::Vector3f*>(p1);
        
        return (p2_3f - p1_3f).squaredNorm();
    }
    
    // Returns the dim'th component of the idx'th point in the class:
    // Since this is inlined and the "dim" argument is typically an immediate value, the
    //  "if/else's" are actually solved at compile time.
    inline float kdtree_get_pt(const size_t idx, int dim) const
    {
        return pts[idx](dim);
    }
    
    // Optional bounding-box computation: return false to default to a standard bbox computation loop.
    //   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
    //   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};



// dynamic adapter
typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
	nanoflann::L2_Simple_Adaptor<float, EigenPointCloud>,
EigenPointCloud,
3 /* dim */
> EigenDynamicPointCloudIndexAdaptor;

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, EigenPointCloud>,
	EigenPointCloud,
	3 /* dim */
> EigenStaticPointCloudAdaptor;

template <typename T>
struct PointCloudPoint
{
    T  x,y,z;
};



template <typename Scalar, typename Point>
struct PointTypedCloud
{
	std::vector<Point>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline Scalar kdtree_get_pt(const size_t idx, int dim) const
	{
		return pts[idx][dim];
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

template<typename Scalar, typename Point>
struct Point_kdTree
{
	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<Scalar, PointTypedCloud<Scalar, Point> >,
		PointTypedCloud<Scalar, Point>,
		3 /* dim */
	> IndexAdaptor;

	PointTypedCloud<Scalar, Point> cloud;
	IndexAdaptor kdTree;

	Point_kdTree(const std::vector<Point>& points) : kdTree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */))
	{
		cloud.pts = points;
		kdTree.buildIndex();
	}
};

struct FVectorPointCloud
{
	std::vector<FVector>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline float kdtree_distance(const float *p1_, const size_t idx_p2, size_t /*size*/) const
	{
		const FVector& p1 = *reinterpret_cast<const FVector*>(p1_);
		const FVector& p2 = pts[idx_p2];

		return FVector::DistSquared(p1, p2);
	}

	// Returns the dim'th component of the idx'th point in the class:
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		return pts[idx].Component(dim);
	}

	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }

};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, FVectorPointCloud >,
	FVectorPointCloud,
	3 /* dim */
> FVectorPointCloudIndexAdaptor;

typedef Point_kdTree<float, FVector> FVector_kdTree;




struct SurfacePosition
{
	SurfacePosition() {}
	SurfacePosition(FVector position, float radius, FSurfaceIndex surfaceIndex) : position(position), radius(radius), surfaceIndex(surfaceIndex) {}

	FVector position;

	FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface;


	float& operator[](int32 i) { return (&position.X)[i]; }

	float radius = 0.0f;
};

struct SurfacePositionCloud
{
	std::vector<SurfacePosition>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t /*size*/) const
	{
		const auto& p2_3f = pts[idx_p2].position;
		const auto& p1_3f = *reinterpret_cast<const FVector*>(p1);

		return (p2_3f - p1_3f).SizeSquared();
	}

	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		return pts[idx].position[dim];
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
	nanoflann::L2_Simple_Adaptor<float, SurfacePositionCloud>,
	SurfacePositionCloud,
	3 /* dim */
> SurfacePositionDynamicCloud;