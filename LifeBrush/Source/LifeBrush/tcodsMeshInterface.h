//
//  tcodsMeshInterface.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-08-21.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#pragma once

#include <map>
#include <unordered_map>
#include <unordered_set>

#include "Algorithm/FastBVH/BVH.h"
#include "Algorithm/tcods/HalfEdge.h"
#include "Algorithm/tcods/MeshIO.h"

#include "Algorithm/PointCloud.hpp"
#include "Algorithm/nanoflann.hpp"

#include "Algorithm/SurfaceIndex.h"

#include "tim_kdtree_face.h"

#include "UniformGrid.h"
#include "ChunkedMarchingCubes.h"

#include "Utility.h"



struct SamplePoint
{
    float x;
    float y;
    float z;

	FSurfaceIndex surfaceIndex;
};

/**
 An interface for navigating the surface of a UStaticMesh instance within the 
 Region Growing Algorithm. The coordinate space of all operations inside of the
 mesh interface are in world space, as determined by the toWorld transform passed in
 tcodsMeshInterface::buildMesh.
 */
struct tcodsMeshInterfaceBase
{
public:
	struct VertexIndex
	{
		int32 sectionIndex;
		int32 vertexIndex;
	};

   
public:
    
	void rebuildBvh();

	/** \param uStaticMeshVertex A vertex directly from the UStaticMesh instance used to
     build the mesh interface.
     */
    auto indexOfVertex(const FVector& uStaticMeshVertex, int32& index_out) -> bool;
	auto indexOfVertex(const FVector& uStaticMeshVertex, VertexIndex& index_out) -> bool;

	bool hasSection(uint32_t section)
	{
		return _sections.find(section) != _sections.end();
	}

	tcods::Mesh& mesh(uint32_t section)
	{
		auto found = _sections.find(section);

		assert(found != _sections.end());

		return *(found->second).get();
	}
    
   
    struct SurfacePoint 
	{ 
		SurfacePoint() : point(0.0f, 0.0f, 0.0f) 
		{
			surfaceIndex.setOffSurface();  
		}

		SurfacePoint(FVector point, FSurfaceIndex surfaceIndex) : point(point), surfaceIndex(surfaceIndex) {}

		FVector point; 

		FSurfaceIndex surfaceIndex;
	};
    
    auto nearestPointOnMesh(const FVector& pointInWorld) -> SurfacePoint;

	auto getIntersection(const FVector& rayOrigin, const FVector& rayDir) -> std::pair<bool, FVector>
    {
        bvh::Ray ray;
        ray.origin = rayOrigin;
        ray.direction = rayDir;
        ray.direction.Normalize();
        
        bvh::IntersectionInfo intersection;
        if( !bvh.getIntersection(ray, intersection, false, true) )
            return std::pair<bool, FVector>(false, FVector::ZeroVector);
        
        FVector hitPoint = intersection.t * ray.direction + ray.origin;
        
        return std::pair<bool, FVector>(true, hitPoint);
    }

    auto getIntersectionAndFace(const FVector& rayOrigin, const FVector& rayDir)
    -> std::pair<bool, SurfacePoint> // hit, point, face index
    {
        bvh::Ray ray;
        ray.origin = rayOrigin;
        ray.direction = rayDir;
        ray.direction.Normalize();
        
        bvh::IntersectionInfo intersection;
		if (!bvh.getIntersection(ray, intersection, false, true))
			return std::make_pair(false, SurfacePoint());
        
        FVector hitPoint = intersection.t * ray.direction + ray.origin;
        
		return std::make_pair(true, SurfacePoint(hitPoint, intersection.triangle.surfaceIndex));
    }
    
    inline auto rotationAtSurfacePoint(const FVector& pointInWorld) -> FQuat
    {
        auto pair = rotationAndNormalAtSurfacePoint(pointInWorld);
        
        return pair.first;
    }
    
    inline auto rotationAndNormalAtSurfacePoint(const FVector& pointInWorld) -> std::pair<FQuat, FVector>
    {
        auto nearest = nearestPointOnMesh(pointInWorld);
        
        return rotationAndNormalAtIndex(nearest.surfaceIndex);
    }
    
    inline auto rotationAndNormalAtIndex(const FSurfaceIndex surfaceIndex) -> std::pair<FQuat, FVector>
    {
        auto face = mesh(surfaceIndex.sectionIndex).faces[surfaceIndex.faceIndex];
        
        using namespace tcods;
		using namespace DDG;
        
        double alpha = face.alpha;
        
        Vector e1, e2; face.frame(e1, e2);
        Vector n = cross(e2, e1);
        Vector u = e1 * cos(alpha) + e2 * sin(alpha);

		// the vector orthogonal to u in the plane of e1, e2
		const double pi_2 = Math::pid() * .5;
		Vector v = e1 * cos( alpha + pi_2 ) + e2 * sin( alpha + pi_2 );
        

        
        FTransform transform = FTransform( unreal(u), unreal(v), unreal(n), FVector::ZeroVector);
        return std::pair<FQuat, FVector>(transform.GetRotation(), unreal(n));
    }
    
    inline auto frameAtNearest(const SurfacePoint& nearest) -> std::pair<FVector,FVector>
    {
        auto face = mesh(nearest.surfaceIndex.sectionIndex).faces[nearest.surfaceIndex.faceIndex];
        
		using namespace tcods;
		using namespace DDG;

        double alpha = face.alpha;
        
        Vector e1, e2;
        face.frame(e1, e2);

        return std::pair<FVector,FVector>(unreal(e1), unreal(e2));
    }
    
public:
    bvh::BVH bvh;
    
    unsigned int sampleSubdivisions = 64;

    std::unordered_map<int, std::vector<int32>> _halfEdgeVertex_to_runtimeMeshVertices;

	std::unordered_map<uint32_t, std::unique_ptr<tcods::Mesh>> _sections;


protected:
    void _initBVH();

	void _initNearestKDTree(FBox limits);

	auto _clear();

	auto _countTriangles() -> size_t;

	auto _extractSections(tcods::MeshIO::MeshData& data)->std::vector<tcods::MeshIO::MeshData>;
	void _floodVisit(const int32 vertex, const std::vector< std::vector<int32> >& verticesConnectedToVertex, std::set<int>& unvisitedVertices, std::unordered_set<int32>& connectedVertices);
	
	void _buildRuntimeMesh(URuntimeMeshComponent * runtimeMesh, bool clipToLimits = true, FBox limits = FBox(EForceInit::ForceInitToZero));
	void _buildRuntimeMeshSection(URuntimeMeshComponent * runtimeMesh, tcods::Mesh& mesh, uint32_t section, FTransform transform, bool clipToLimits = true, FBox limits = FBox());

protected:
    std::unordered_map<FVector, VertexIndex> _vertexLookup;
    
	std::unique_ptr<vcg::KdTreeFace> _faceTree;

	FBox _worldLimits = FBox(EForceInit::ForceInitToZero);
};

struct tcodsUStaticMeshInterface : public tcodsMeshInterfaceBase
{
	auto buildMesh(UStaticMesh& uMesh_in, const FTransform& toWorld, FBox sampleLimits = FBox(), class URuntimeMeshComponent * runtimeMeshToCopyInto = nullptr) -> bool;

protected:
	struct MeshDataAndVertexLookup
	{
		tcods::MeshIO::MeshData meshData;
		std::unordered_map<FVector, VertexIndex> vertexLookup;
	};

	auto _toMeshData(UStaticMesh& uMesh, const FTransform toWorld)->MeshDataAndVertexLookup;
};

struct tcodsChunkedGridMeshInterface : public tcodsMeshInterfaceBase
{
public:
	void buildMesh(ChunkGrid<float>& chunkGrid, const FTransform& toWorld, float isoLevel, FBox limits, URuntimeMeshComponent * runtimeMesh = nullptr);


protected:
	auto _toMergedMeshData(ChunkGrid<float> &chunkGrid, float isoLevel, FTransform toWorld)->tcods::MeshIO::MeshData;

	void _buildSectionsFromMeshData(tcods::MeshIO::MeshData& data);

	tcods::MeshIO::MeshData _toTcods(SmoothMeshFactory& factory);

	SmoothMeshFactory _clean(SmoothMeshFactory& mergedFactory);
};

