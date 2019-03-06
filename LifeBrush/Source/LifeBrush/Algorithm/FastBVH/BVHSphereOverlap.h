//#ifndef BVH_h
//#define BVH_h
//
//#include <vector>
//#include <array>
//#include <stdint.h>
//#include "Object.h"
//
//#include "Eigen/Dense"
//#include "Utility.h"
//
//namespace bvhOverlap
//{
//    struct Ray
//    {
//        FVector origin = FVector::ZeroVector;
//        FVector direction = FVector::ZeroVector;
//    };
//    
//    struct Sphere
//    {
//    public:
//        Sphere(FVector position, float radius) : position(position), radius(radius) {};
//        
//        FVector position;
//        float radius;
//        
//        FVector getCentroid()
//        {
//            return position;
//        }
//        
//        FBox getBBox()
//        {
//            FVector halfExtents(radius);
//            
//            FVector min = position - halfExtents;
//            FVector max = position + halfExtents;
//            
//            return FBox(min, max);
//        }
//    };
//    
//    struct IntersectionInfo
//    {
//        float t;                    // Intersection distance along the ray
//        Sphere triangle;    // Object that was hit
//    };
//    
//    //! Node descriptor for the flattened tree
//    struct BVHFlatNode
//    {
//        FBox bbox;
//        uint32_t start = 0;
//        uint32_t nPrims = 0;
//        uint32_t rightOffset = 0;
//    };
//    
//    //! \author Brandon Pelfrey
//    //! A Bounding Volume Hierarchy system for fast Ray-Object intersection tests
//    class BVH
//    {
//        uint32_t nNodes = 0;
//        uint32_t nLeafs = 0;
//        uint32_t leafSize = 0;
//        
//        std::vector<Sphere> build_prims;
//        
//        // Fast Traversal System
//        std::vector<BVHFlatNode> flatTree;
//        
//    public:
//        //! Build the BVH tree out of build_prims
//        void build(std::vector<Sphere>& objects, uint32_t leafSize=4);
//
//        bool getIntersection(const Ray& ray, IntersectionInfo& intersection, bool occlusion = false, bool doubleSide = false) const ;
//        std::vector<IntersectionInfo> getIntersections(const Ray& ray, bool doubleSided = false);
//    };
//}
//
//#endif
