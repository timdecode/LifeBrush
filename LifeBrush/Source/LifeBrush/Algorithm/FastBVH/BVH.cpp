#include "LifeBrush.h"

#include <algorithm>
#include <limits>

#include "BVH.h"

using namespace bvh;

// -----------------------------------------------------------------------------
/// BVH Structs
// -----------------------------------------------------------------------------

//! Node for storing state information during traversal.
struct BVHTraversal
{
    uint32_t i = 0; // Node
    float mint = 0; // Minimum hit time for this node.
    BVHTraversal() { }
    BVHTraversal(int _i, float _mint) : i(_i), mint(_mint) { }
};

struct BVHBuildEntry {
    // If non-zero then this is the index of the parent. (used in offsets)
    uint32_t parent = 0;
    // The range of objects in the object list covered by this node.
    uint32_t start = 0, end = 0;
};

// -----------------------------------------------------------------------------
/// Intersection Tests
// -----------------------------------------------------------------------------

// This code is based on GLM's intersectRayTriangle method in intersect.inl
// https://github.com/g-truc/glm/blob/master/glm/gtx/intersect.inl
// 
// Adapted by Timothy Davison, September 5, 2015 for the Unreal Engine
//
// baryPosition.z contains the time of intersect of the ray with this triangle.
bool intersects(FVector const & origin, FVector const & direction, FVector const & v0, FVector const & v1, FVector const &v2, FVector& baryPosition)
{
    auto e1 = v1 - v0;
    auto e2 = v2 - v0;
    
    auto p = FVector::CrossProduct(direction, e2);
    
    float a = FVector::DotProduct(e1, p);
    
    if( a < std::numeric_limits<float>::epsilon() )
        return false;
    
    float f = 1.0f / a;
    
    auto s = origin - v0;
    baryPosition.X = f * FVector::DotProduct(s, p);
    if( baryPosition.X < 0.0f )
        return false;
    if( baryPosition.X > 1.0f )
        return false;
    
    auto q = FVector::CrossProduct(s, e1);
    baryPosition.Y = f * FVector::DotProduct(direction, q);
    if( baryPosition.Y < 0.0f )
        return false;
    if( baryPosition.Y > 1.0f )
        return false;
    
    baryPosition.Z = f * FVector::DotProduct(e2, q);
    
    return baryPosition.Z >= 0.0f;
}

// Returns true if the ray intersects the triangle's front or back face.
// t_out will be negative if we have a backface intersection.
bool rayTriangle(FVector const & origin, FVector const & direction, FVector const & v0, FVector const & v1, FVector const &v2, float& t_out)
{
    FVector u = v1 - v0;
    FVector v = v2 - v0;
    
    FVector n = FVector::CrossProduct(u,v);
    
    if( n.Size() == 0.0f )
        return false;
    
    FVector w = origin - v0;

    float a = -FVector::DotProduct(n, w);
    float b = FVector::DotProduct(n, direction);
    
    if( FMath::Abs(b) < 0.000001f )
        return false;
    
    float r = a / b;
    
//    if( r < 0.0f )
//        return false;
//    
    FVector i = origin + r * direction;
    
    float uu, uv, vv, wu, wv, d;
    
    uu = FVector::DotProduct(u,u);
    uv = FVector::DotProduct(u,v);
    vv = FVector::DotProduct(v,v);
    
    w = i - v0;
    
    wu = FVector::DotProduct(w,u);
    wv = FVector::DotProduct(w,v);
    
    d = uv * uv - uu * vv;
    
    float s,t;
    
    s = (uv * wv - vv * wu) / d;
    if( s < 0.0f || s > 1.0f )
        return false;
    
    t = (uv * wu - uu * wv) / d;
    if( t < 0.0f || (s + t) > 1.0f )
        return false;
    
    t_out = r;
    
    return true;
}

bool intersects(FVector const & origin, FVector const & direction, const FBox& box, float& tNear, float& tFar)
{
    FVector invDirection( 1.0f / direction.X, 1.0f / direction.Y, 1.0f / direction.Z);
    
    FVector tBot = invDirection * (box.Min - origin);
    FVector tTop = invDirection * (box.Max - origin);
    
    FVector tMin = tTop.ComponentMin(tBot);
    FVector tMax = tTop.ComponentMax(tBot);
    
    tNear = tMin.GetMax();
    tFar = tMax.GetMin();
    
    return !(tNear > tFar) && tFar > 0.0f;
}

// -----------------------------------------------------------------------------
/// BVH Class
// -----------------------------------------------------------------------------


//! - Compute the nearest intersection of all objects within the tree.
//! - Return true if hit was found, false otherwise.
//! - In the case where we want to find out of there is _ANY_ intersection at all,
//!   set occlusion == true, in which case we exit on the first hit, rather
//!   than find the closest.
bool BVH::getIntersection(const Ray& ray, IntersectionInfo& intersection, bool occlusion, bool doubleSided) const
{
	if(nNodes == 0)
		return false;

    intersection.t = std::numeric_limits<float>::max();
    
    float bbhits[4];
    int32_t closer, other;
    
    // Working set
    BVHTraversal todo[64];
    int32_t stackptr = 0;
    
    // "Push" on the root node to the working set
    todo[stackptr].i = 0;
    todo[stackptr].mint = -std::numeric_limits<float>::max();
    
    bool didHitOnce = false;
    
    while(stackptr>=0)
    {
        // Pop off the next node to work on.
        int ni = todo[stackptr].i;
        float near_ = todo[stackptr].mint;
        stackptr--;
        const BVHFlatNode &node = flatTree[ ni ];
        
        // If this node is further than the closest found intersection, continue
        if(near_ > intersection.t)
            continue;
        
        // Is leaf -> Intersect
        if( node.rightOffset == 0 )
        {
            for(uint32_t o=0;o<node.nPrims;++o)
            {
                IntersectionInfo current;
                
                const Triangle& obj = build_prims[node.start+o];
                
                float t;
                bool hit = false;
                if( rayTriangle(ray.origin, ray.direction, obj.vertices[0], obj.vertices[1], obj.vertices[2], t) )
                    hit = doubleSided || t >= 0.0f;
                
                
                current.t = t;
                current.triangle = obj;
                
                if (hit)
                {
                    didHitOnce = true;
                    
                    // If we're only looking for occlusion, then any hit is good enough
                    if(occlusion)
                    {
                        intersection = current;
                        return true;
                    }
                    
                    // Otherwise, keep the closest intersection only
                    if( current.t < intersection.t )
                        intersection = current;
                }
            }
            
        }
        else
        {
            // Not a leaf
            bool hitc0 = intersects(ray.origin, ray.direction, flatTree[ni+1].bbox, bbhits[0], bbhits[1]);
            bool hitc1 = intersects(ray.origin, ray.direction, flatTree[ni+node.rightOffset].bbox, bbhits[2], bbhits[3]);
            
            
            // Did we hit both nodes?
            if(hitc0 && hitc1)
            {
                // We assume that the left child is a closer hit...
                closer = ni+1;
                other = ni+node.rightOffset;
                
                // ... If the right child was actually closer, swap the relavent values.
                if(bbhits[2] < bbhits[0])
                {
                    std::swap(bbhits[0], bbhits[2]);
                    std::swap(bbhits[1], bbhits[3]);
                    std::swap(closer,other);
                }
                
                // It's possible that the nearest object is still in the other side, but we'll
                // check the further-awar node later...
                
                // Push the farther first
                todo[++stackptr] = BVHTraversal(other, bbhits[2]);
                
                // And now the closer (with overlap test)
                todo[++stackptr] = BVHTraversal(closer, bbhits[0]);
            }
            else if (hitc0)
                todo[++stackptr] = BVHTraversal(ni+1, bbhits[0]);
            else if(hitc1)
                todo[++stackptr] = BVHTraversal(ni + node.rightOffset, bbhits[2]);
            
        }
    }
    
    return didHitOnce;
}

std::vector<IntersectionInfo> BVH::getIntersections(const Ray& ray, bool doubleSided)
{
    std::vector<IntersectionInfo> intersections;

	if(nNodes == 0)
		return intersections;
    
    IntersectionInfo intersection;
    
    float bbhits[4];
    int32_t closer, other;
    
    // Working set
    BVHTraversal todo[64];
    int32_t stackptr = 0;
    
    // "Push" on the root node to the working set
    todo[stackptr].i = 0;
    todo[stackptr].mint = -std::numeric_limits<float>::max();
    
    while( stackptr >= 0 )
    {
        // Pop off the next node to work on.
        int ni = todo[stackptr].i;
        float near_ = todo[stackptr].mint;
        stackptr--;
        const BVHFlatNode &node = flatTree[ ni ];
        
        // Is leaf -> Intersect
        if( node.rightOffset == 0 )
        {
            for(uint32_t o=0;o<node.nPrims;++o)
            {
                const Triangle& triangle = build_prims[node.start+o];
                
                float t;
                bool hit = false;
                if( rayTriangle(ray.origin, ray.direction, triangle.vertices[0], triangle.vertices[1], triangle.vertices[2], t) )
                    hit = doubleSided || t >= 0.0f;
                
                intersection.t = t;
                intersection.triangle = triangle;
                
                if( hit )
                    intersections.push_back(intersection);
            }
        }
        else
        {
            // Not a leaf
            bool hitc0 = intersects(ray.origin, ray.direction, flatTree[ni+1].bbox, bbhits[0], bbhits[1]);
            bool hitc1 = intersects(ray.origin, ray.direction, flatTree[ni+node.rightOffset].bbox, bbhits[2], bbhits[3]);
            
            // Did we hit both nodes?
            if(hitc0 && hitc1)
            {
                // We assume that the left child is a closer hit...
                closer = ni+1;
                other = ni+node.rightOffset;
                
                // ... If the right child was actually closer, swap the relavent values.
                if(bbhits[2] < bbhits[0])
                {
                    std::swap(bbhits[0], bbhits[2]);
                    std::swap(bbhits[1], bbhits[3]);
                    std::swap(closer,other);
                }
                
                // It's possible that the nearest object is still in the other side, but we'll
                // check the further-awar node later...
                
                // Push the farther first
                todo[++stackptr] = BVHTraversal(other, bbhits[2]);
                
                // And now the closer (with overlap test)
                todo[++stackptr] = BVHTraversal(closer, bbhits[0]);
            }
            else if (hitc0)
                todo[++stackptr] = BVHTraversal(ni+1, bbhits[0]);
            else if(hitc1)
                todo[++stackptr] = BVHTraversal(ni + node.rightOffset, bbhits[2]);
        }
    }

    return intersections;
}

/*! Build the BVH, given an input data set
 *  - Handling our own stack is quite a bit faster than the recursive style.
 *  - Each build stack entry's parent field eventually stores the offset
 *    to the parent of that node. Before that is finally computed, it will
 *    equal exactly three other values. (These are the magic values Untouched,
 *    Untouched-1, and TouchedTwice).
 *  - The partition here was also slightly faster than std::partition.
 */
void BVH::build(std::vector<Triangle>& objects, uint32_t leafSize_in)
{
    build_prims = objects;
    leafSize = leafSize_in;
    nNodes = 0;
    nLeafs = 0;

	if (objects.empty()) return;

    BVHBuildEntry todo[128];
    uint32_t stackptr = 0;
    const uint32_t Untouched    = 0xffffffff;
    const uint32_t TouchedTwice = 0xfffffffd;
    
    // Push the root
    todo[stackptr].start = 0;
    todo[stackptr].end = build_prims.size();
    todo[stackptr].parent = 0xfffffffc;
    stackptr++;
    
    BVHFlatNode node;
    std::vector<BVHFlatNode> buildnodes;
    buildnodes.reserve(build_prims.size()*2);
    
    while(stackptr > 0)
    {
        // Pop the next item off of the stack
        BVHBuildEntry &bnode( todo[--stackptr] );
        uint32_t start = bnode.start;
        uint32_t end = bnode.end;
        uint32_t nPrims = end - start;
        
        nNodes++;
        node.start = start;
        node.nPrims = nPrims;
        node.rightOffset = Untouched;
        
        // Calculate the bounding box for this node
        FBox bb = build_prims[start].getBBox();
        
        FVector const centroid = build_prims[start].getCentroid();
        FBox bc( centroid, centroid );
        
        for(uint32_t p = start+1; p < end; ++p)
        {
            bb += build_prims[p].getBBox();
            bc += build_prims[p].getCentroid();
        }
        node.bbox = bb;
        
        // If the number of primitives at this point is less than the leaf
        // size, then this will become a leaf. (Signified by rightOffset == 0)
        if(nPrims <= leafSize) {
            node.rightOffset = 0;
            nLeafs++;
        }
        
        buildnodes.push_back(node);
        
        // Child touches parent...
        // Special case: Don't do this for the root.
        if(bnode.parent != 0xfffffffc) {
            buildnodes[bnode.parent].rightOffset --;
            
            // When this is the second touch, this is the right child.
            // The right child sets up the offset for the flat tree.
            if( buildnodes[bnode.parent].rightOffset == TouchedTwice ) {
                buildnodes[bnode.parent].rightOffset = nNodes - 1 - bnode.parent;
            }
        }
        
        // If this is a leaf, no need to subdivide.
        if(node.rightOffset == 0)
            continue;
        
        // Set the split dimensions
        const FVector size = bc.GetSize();
        uint32_t split_dim = 0;
        if( size.Y > size.X ) split_dim = 1;
        if( size.Z > size.Y ) split_dim = 2;
        
        // Split on the center of the longest axis
        float split_coord = .5f * (bc.Min[split_dim] + bc.Max[split_dim]);
        
        // Partition the list of objects on this split
        uint32_t mid = start;
        for(uint32_t i=start;i<end;++i)
        {
            if( build_prims[i].getCentroid()[split_dim] < split_coord )
            {
                std::swap( build_prims[i], build_prims[mid] );
                ++mid;
            }
        }
        
        // If we get a bad split, just choose the center...
        if(mid == start || mid == end)
            mid = start + (end-start)/2;
        
        // Push right child
        todo[stackptr].start = mid;
        todo[stackptr].end = end;
        todo[stackptr].parent = nNodes-1;
        stackptr++;
        
        // Push left child
        todo[stackptr].start = start;
        todo[stackptr].end = mid;
        todo[stackptr].parent = nNodes-1;
        stackptr++;
    }
    
    // Copy the temp node data to a flat array
    flatTree = buildnodes;
}

void bvh::BVH::clear()
{
	build_prims.clear(); 
	flatTree.clear();

	nNodes = 0;
	nLeafs = 0;
	leafSize = 0;
}

