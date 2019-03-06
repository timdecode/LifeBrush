// Copyright(C) David W. Jeske, 2014, and released to the public domain. 
//
// Dynamic BVH (Bounding Volume Hierarchy) using incremental refit and tree-rotations
//
// initial BVH build based on: Bounding Volume Hierarchies (BVH) – A brief tutorial on what they are and how to implement them
//              http://www.3dmuve.com/3dmblog/?p=182
//
// Dynamic Updates based on: "Fast, Effective BVH Updates for Animated Scenes" (Kopta, Ize, Spjut, Brunvand, David, Kensler)
//              http://www.cs.utah.edu/~thiago/papers/rotations.pdf
//
// see also:  Space Partitioning: Octree vs. BVH
//            http://thomasdiewald.com/blog/?p=1488
//
//

#include <vector>
#include <set>

// TODO: handle merge/split when LEAF_OBJ_MAX > 1 and objects move
// TODO: add sphere traversal

namespace dynamicBVH
{
    enum Axis
    {
        X = 0,
        Y = 1,
        Z = 2,
        End = 3
    };
    
    class AABB
    {
        Eigen::Vector3f min;
        Eigen::Vector3f max;
    };

    template <typename ObjectT>
    class BVH
    {
    public:
        BVHNode<ObjectT> rootBVH;

        int LEAF_OBJ_MAX;
        int nodeCount = 0;
        
        std::set<ObjectT> refitNodes;

        /// @param LEAF_OBJ_MAX  currently this must be 1 to use dynamic BVH updates
        public BVH(std::vector<ObjectT>& objects, int LEAF_OBJ_MAX = 1)
        {
            this.LEAF_OBJ_MAX = LEAF_OBJ_MAX;
            nodeAdaptor.setBVH(this);
            this.nAda = nodeAdaptor;
            
            if( objects.size() > 0 )
                rootBVH = BVHNode<ObjectT>(objects);
            else
                rootBVH = BVHNode<ObjectT>();
        }

        // internal functional traversal...
        private void _traverse(BVHNode<ObjectT>& curNode, NodeTest& hitTest, std::vector<BVHNode<ObjectT>>& hitlist)
        {
            if (curNode == null) { return; }
            
            if (hitTest(curNode.box))
            {
                hitlist.push_back(curNode);
                
                _traverse(curNode.left,hitTest,hitlist);
                _traverse(curNode.right,hitTest,hitlist);
            }
        }

        // public interface to traversal..
        std::vector<BVHNode<ObjectT>> traverse(NodeTest& hitTest)
        {
            auto hits = std::vector<BVHNode<ObjectT>>();
            
            _traverse(rootBVH,hitTest,hits);
            
            return hits;
        }
        
        // left in for compatibility..
        std::vector<BVHNode<ObjectT>> traverseRay(SSRay ray)
        {
            float tnear = 0f, tfar = 0f;

            return traverse( box => OpenTKHelper.intersectRayAABox1(ray,box,ref tnear, ref tfar) );
        }

        public List<BVHNode<ObjectT>> traverse(SSRay ray) {
            float tnear = 0f, tfar = 0f;

            return traverse( box => OpenTKHelper.intersectRayAABox1(ray,box,ref tnear, ref tfar) );
        }
        public List<BVHNode<ObjectT>> traverse(SSAABB volume) {
            return traverse( box => box.IntersectsAABB(volume) );            
        }

        /// <summary>
        /// Call this to batch-optimize any object-changes notified through 
        /// BVHNode.refit_ObjectChanged(..). For example, in a game-loop, 
        /// call this once per frame.
        /// </summary>

        public void optimize()
        {
            assert(LEAF_OBJ_MAX == 1, "LEAF_OBJ_MAX must be 1 for dynamic updates.");
                  
            while (refitNodes.Count > 0)
            {
                int maxdepth = refitNodes.Max( n => n.depth );
            
                var sweepNodes = refitNodes.Where( n => n.depth == maxdepth ).ToList();
                sweepNodes.ForEach( n => refitNodes.Remove(n) );

                sweepNodes.ForEach( n => n.tryRotate(this) );                
            }            
        }

        public void addObject(ObjectT& newOb)
        {
            SSAABB box = SSAABB.FromSphere(nAda.objectpos(newOb),nAda.radius(newOb));
            float boxSAH = rootBVH.SAH(ref box);
            
            rootBVH.add
            
            
            rootBVH.addObject(nAda,newOb, ref box, boxSAH);
        }

        public void removeObject(GO newObj)
        {
            var leaf = nAda.getLeaf(newObj);
            leaf.removeObject(nAda,newObj);
        }

        public int countBVHNodes()
        {
            return rootBVH.countBVHNodes();
        }
    }   
}
