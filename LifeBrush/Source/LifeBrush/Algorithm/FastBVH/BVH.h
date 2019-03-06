#pragma once

#include <vector>
#include <array>
#include <stdint.h>
#include "Object.h"
#include "Algorithm/SurfaceIndex.h"

namespace bvh
{
    struct Ray
    {
        FVector origin = FVector::ZeroVector;
        FVector direction = FVector::ZeroVector;
    };
    
    struct Triangle
    {
    public:
        Triangle() {}
        Triangle(std::array<FVector,3>& vertices, uint32_t faceIndex, uint32_t sectionIndex) : vertices(vertices), surfaceIndex(surfaceIndex) {}
        
        std::array<FVector, 3> vertices;

		FSurfaceIndex surfaceIndex;
        
        FVector getCentroid()
        {
            return (vertices[0] + vertices[1] + vertices[2]) / 3.0f;
        }
        
        FBox getBBox()
        {
            return FBox(&vertices[0], 3);
        }
    };
    
    struct IntersectionInfo
    {
        float t;                    // Intersection distance along the ray
        Triangle triangle;    // Object that was hit
    };
    
    //! Node descriptor for the flattened tree
    struct BVHFlatNode
    {
        FBox bbox;
        uint32_t start = 0;
        uint32_t nPrims = 0;
        uint32_t rightOffset = 0;
    };
    
    //! \author Brandon Pelfrey
    //! A Bounding Volume Hierarchy system for fast Ray-Object intersection tests
    class BVH
    {
        uint32_t nNodes = 0;
        uint32_t nLeafs = 0;
        uint32_t leafSize = 0;
        
        std::vector<Triangle> build_prims;
        
        // Fast Traversal System
        std::vector<BVHFlatNode> flatTree;
        
    public:
        //! Build the BVH tree out of build_prims
        void build(std::vector<Triangle>& objects, uint32_t leafSize=4);

		void clear();

		size_t size() { return build_prims.size(); }

        bool getIntersection(const Ray& ray, IntersectionInfo& intersection, bool occlusion = false, bool doubleSide = false) const ;
        std::vector<IntersectionInfo> getIntersections(const Ray& ray, bool doubleSided = false);
    };
};


