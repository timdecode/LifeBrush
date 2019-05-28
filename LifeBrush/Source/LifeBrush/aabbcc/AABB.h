/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>
  Copyright (c) 2019 Timothy Davison <timtimmy@gmail.com>

  Tim: Converted Lester's version to templated code. Replaced dynamic std::vector's with std::array.

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#ifndef _AABB_H
#define _AABB_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
#include <array>

/// Null node flag.
const unsigned int NULL_NODE = 0xffffffff;

namespace aabb
{
    /*! \brief The axis-aligned bounding box object.

        Axis-aligned bounding boxes (AABBs) store information for the minimum
        orthorhombic bounding-box for an object. Support is provided for
        dimensions >= 2. (In 2D the bounding box is either a rectangle,
        in 3D it is a rectangular prism.)

        Class member functions provide functionality for merging AABB objects
        and testing overlap with other AABBs.
     */
	template<unsigned int DIM, typename SCALAR>
    class AABB
    {
    public:
		typedef std::array<SCALAR, DIM> Position_t;

        /// Constructor.
        AABB();


        //! Constructor.
        /*! \param lowerBound_
                The lower bound in each dimension.

            \param upperBound_
                The upper bound in each dimension.
         */
        AABB(const float * lowerBound, const float * upperBound);

		AABB(const float * origin, float radius);

        /// Compute the surface area of the box.
        SCALAR computeSurfaceArea() const;

        /// Get the surface area of the box.
        SCALAR getSurfaceArea() const;

        //! Merge two AABBs into this one.
        /*! \param aabb1
                A reference to the first AABB.

            \param aabb2
                A reference to the second AABB.
         */
        void merge(const AABB&, const AABB&);

        //! Test whether the AABB is contained within this one.
        /*! \param aabb
                A reference to the AABB.

            \return
                Whether the AABB is fully contained.
         */
        bool contains(const AABB&) const;

        //! Test whether the AABB overlaps this one.
        /*! \param aabb
                A reference to the AABB.

            \param touchIsOverlap
                Does touching constitute an overlap?

            \return
                Whether the AABB overlaps.
         */
        bool overlaps(const AABB&, bool touchIsOverlap) const;

        //! Compute the centre of the AABB.
        /*! \returns
                The position vector of the AABB centre.
         */
		Position_t computeCentre();

        /// Lower bound of AABB in each dimension.
        Position_t lowerBound;

        /// Upper bound of AABB in each dimension.
		Position_t upperBound;

        /// The position of the AABB centre.
		Position_t centre;

        /// The AABB's surface area.
        SCALAR surfaceArea;
    };






    /*! \brief A node of the AABB tree.

        Each node of the tree contains an AABB object which corresponds to a
        particle, or a group of particles, in the simulation box. The AABB
        objects of individual particles are "fattened" before they are stored
        to avoid having to continually update and rebalance the tree when
        displacements are small.

        Nodes are aware of their position within in the tree. The isLeaf member
        function allows the tree to query whether the node is a leaf, i.e. to
        determine whether it holds a single particle.
     */
	template<unsigned int DIM, typename SCALAR>
    struct Node
    {
		typedef AABB<DIM, SCALAR> AABB_t;

        /// Constructor.
        Node();

        /// The fattened axis-aligned bounding box.
		AABB_t aabb;

        /// Index of the parent node.
        unsigned int parent;

        /// Index of the next node.
        unsigned int next;

        /// Index of the left-hand child.
        unsigned int left;

        /// Index of the right-hand child.
        unsigned int right;

        /// Height of the node. This is 0 for a leaf and -1 for a free node.
        int height;

        /// The index of the particle that the node contains (leaf nodes only).
        unsigned int particle;

        //! Test whether the node is a leaf.
        /*! \return
                Whether the node is a leaf node.
         */
        bool isLeaf() const;
    };

    /*! \brief The dynamic AABB tree.

        The dynamic AABB tree is a hierarchical data structure that can be used
        to efficiently query overlaps between objects of arbitrary shape and
        size that lie inside of a simulation box. Support is provided for
        periodic and non-periodic boxes, as well as boxes with partial
        periodicity, e.g. periodic along specific axes.
     */
    template<unsigned int DIM, typename SCALAR>
    class Tree
    {
    public:
		typedef AABB<DIM, SCALAR> AABB_t;
		typedef Node<DIM, SCALAR> Node_t;

		typedef std::array<SCALAR, DIM> Position_t;

        //! Constructor (non-periodic).
        /*! \param dimension_
                The dimensionality of the system.

            \param skinThickness_
                The skin thickness for fattened AABBs, as a fraction
                of the AABB base length.

            \param nParticles
                The number of particles (for fixed particle number systems).

            \param touchIsOverlap
                Does touching count as overlapping in query operations?
         */
        Tree(SCALAR skinThickness_ = 0.05,
            unsigned int nParticles = 16, bool touchIsOverlap=true);

        //! Constructor (custom periodicity).
        /*! \param dimension_
                The dimensionality of the system.

            \param skinThickness_
                The skin thickness for fattened AABBs, as a fraction
                of the AABB base length.

            \param periodicity_
                Whether the system is periodic in each dimension.

            \param boxSize_
                The size of the simulation box in each dimension.

            \param nParticles
                The number of particles (for fixed particle number systems).

            \param touchIsOverlap
                Does touching count as overlapping in query operations?
         */
        Tree( SCALAR, const std::array<bool,DIM>& periodicity, const float * boxSize,
            unsigned int nParticles = 16, bool touchIsOverlap=true);

        //! Set the periodicity of the simulation box.
        /*! \param periodicity_
                Whether the system is periodic in each dimension.
         */
        void setPeriodicity(const std::array<bool,DIM>& periodicity);

        //! Set the size of the simulation box.
        /*! \param boxSize_
                The size of the simulation box in each dimension.
         */
        void setBoxSize(const float * boxSize);

        //! Insert a particle into the tree (point particle).
        /*! \param index
                The index of the particle.

            \param position
                The position vector of the particle.

            \param radius
                The radius of the particle.
         */
        void insertParticle(unsigned int index, const float * position, SCALAR radius);

        //! Insert a particle into the tree (arbitrary shape with bounding box).
        /*! \param index
                The index of the particle.

            \param lowerBound
                The lower bound in each dimension.

            \param upperBound
                The upper bound in each dimension.
         */
        void insertParticle(unsigned int index, const float * lowerBound, const float * upperBound);

        /// Return the number of particles in the tree.
        unsigned int nParticles();

        //! Remove a particle from the tree.
        /*! \param particle
                The particle index (particleMap will be used to map the node).
         */
        void removeParticle(unsigned int index);

		bool containsParticle(unsigned int index);

        /// Remove all particles from the tree.
        void removeAll();

        //! Update the tree if a particle moves outside its fattened AABB.
        /*! \param particle
                The particle index (particleMap will be used to map the node).

            \param position
                The position vector of the particle.

            \param radius
                The radius of the particle.

            \param alwaysReinsert
                Always reinsert the particle, even if it's within its old AABB (default:false)

            \return
                Whether the particle was reinserted.
         */
        bool updateParticle(unsigned int index, const float * position, SCALAR radius, bool alwaysReinsert=false);

        //! Update the tree if a particle moves outside its fattened AABB.
        /*! \param particle
                The particle index (particleMap will be used to map the node).

            \param lowerBound
                The lower bound in each dimension.

            \param upperBound
                The upper bound in each dimension.

            \param alwaysReinsert
                Always reinsert the particle, even if it's within its old AABB (default: false)
         */
        bool updateParticle(unsigned int index, const float * lowerBound, const float * upperBound, bool alwaysReinsert=false);

        //! Query the tree to find candidate interactions for a particle.
        /*! \param particle
                The particle index.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> query(unsigned int index);

        //! Query the tree to find candidate interactions for an AABB.
        /*! \param particle
                The particle index.

            \param aabb
                The AABB.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> query(unsigned int index, const AABB_t& aabbQuery);

        //! Query the tree to find candidate interactions for an AABB.
        /*! \param aabb
                The AABB.

            \return particles
                A vector of particle indices.
         */
        std::vector<unsigned int> query(const AABB_t& aabbQuery);

		template<typename Callback>
		void query(const AABB_t& aabbQuery, Callback callback)
		{
			query_ref(aabbQuery, callback);
		}

		// \param callback A function that accepts a particle index (unsigned int) whose AABB overlaps
		// with the query AABB. It should return true if the query should continue, and false, if the
		// query should end.
		template<typename Callback>
		void query_ref(const AABB_t& aabbQuery, Callback& callback)
		{
			if (nParticles() == 0)
				return;

			std::vector<unsigned int> stack;

			stack.reserve(256);
			stack.push_back(root);

			while (stack.size() > 0)
			{
				unsigned int node = stack.back();
				stack.pop_back();

				// Copy the AABB.
				auto nodeAABB = nodes[node].aabb;

				if (node == NULL_NODE) continue;

				if (isPeriodic)
				{
					Position_t separation;
					Position_t shift;

					for (unsigned int i = 0; i < DIM; i++)
						separation[i] = nodeAABB.centre[i] - aabbQuery.centre[i];

					bool isShifted = minimumImage(separation, shift);

					// Shift the AABB.
					if (isShifted)
					{
						for (unsigned int i = 0; i < DIM; i++)
						{
							nodeAABB.lowerBound[i] += shift[i];
							nodeAABB.upperBound[i] += shift[i];
						}
					}
				}

				// Test for overlap between the AABBs.
				if (aabbQuery.overlaps(nodeAABB, touchIsOverlap))
				{
					// Check that we're at a leaf node.
					if (nodes[node].isLeaf())
					{
						if( !callback(nodes[node].particle) )
							return;
					}
					else
					{
						stack.push_back(nodes[node].left);
						stack.push_back(nodes[node].right);
					}
				}
			}
		}



        //! Get a particle AABB.
        /*! \param particle
                The particle index.
         */
        const AABB_t& getAABB(unsigned int index);

        //! Get the height of the tree.
        /*! \return
                The height of the binary tree.
         */
        unsigned int getHeight() const;

        //! Get the number of nodes in the tree.
        /*! \return
                The number of nodes in the tree.
         */
        unsigned int getNodeCount() const;

        //! Compute the maximum balancance of the tree.
        /*! \return
                The maximum difference between the height of two
                children of a node.
         */
        unsigned int computeMaximumBalance() const;

        //! Compute the surface area ratio of the tree.
        /*! \return
                The ratio of the sum of the node surface area to the surface
                area of the root node.
         */
        SCALAR computeSurfaceAreaRatio() const;

        /// Validate the tree.
        void validate() const;

        /// Rebuild an optimal tree.
        void rebuild();

    private:
        /// The index of the root node.
        unsigned int root;

        /// The dynamic tree.
        std::vector<Node_t> nodes;

        /// The current number of nodes in the tree.
        unsigned int nodeCount;

        /// The current node capacity.
        unsigned int nodeCapacity;

        /// The position of node at the top of the free list.
        unsigned int freeList;

        /// Whether the system is periodic along at least one axis.
        bool isPeriodic;

        /// The skin thickness of the fattened AABBs, as a fraction of the AABB base length.
        SCALAR skinThickness;

        /// Whether the system is periodic along each axis.
        std::array<bool,DIM> periodicity;

        /// The size of the system in each dimension.
        std::array<SCALAR,DIM> boxSize;

        /// The position of the negative minimum image.
        std::array<SCALAR,DIM> negMinImage;

        /// The position of the positive minimum image.
        std::array<SCALAR,DIM> posMinImage;

        /// A map between particle and node indices.
        std::map<unsigned int, unsigned int> particleMap;

        /// Does touching count as overlapping in tree queries?
        bool touchIsOverlap;

        //! Allocate a new node.
        /*! \return
                The index of the allocated node.
         */
        unsigned int allocateNode();

        //! Free an existing node.
        /*! \param node
                The index of the node to be freed.
         */
        void freeNode(unsigned int);

        //! Insert a leaf into the tree.
        /*! \param leaf
                The index of the leaf node.
         */
        void insertLeaf(unsigned int);

        //! Remove a leaf from the tree.
        /*! \param leaf
                The index of the leaf node.
         */
        void removeLeaf(unsigned int);

        //! Balance the tree.
        /*! \param node
                The index of the node.
         */
        unsigned int balance(unsigned int);

        //! Compute the height of the tree.
        /*! \return
                The height of the entire tree.
         */
        unsigned int computeHeight() const;

        //! Compute the height of a sub-tree.
        /*! \param node
                The index of the root node.

            \return
                The height of the sub-tree.
         */
        unsigned int computeHeight(unsigned int) const;

        //! Assert that the sub-tree has a valid structure.
        /*! \param node
                The index of the root node.
         */
        void validateStructure(unsigned int) const;

        //! Assert that the sub-tree has valid metrics.
        /*! \param node
                The index of the root node.
         */
        void validateMetrics(unsigned int) const;

        //! Apply periodic boundary conditions.
        /* \param position
                The position vector.
         */
        void periodicBoundaries(float * position);

        //! Compute minimum image separation.
        /*! \param separation
                The separation vector.

            \param shift
                The shift vector.

            \return
                Whether a periodic shift has been applied.
         */
        bool minimumImage(Position_t&, Position_t&);
    };

	template<unsigned int DIM, typename SCALAR>
	AABB<DIM, SCALAR>::AABB()
	{

	}


	template<unsigned int DIM, typename SCALAR>
	AABB<DIM, SCALAR>::AABB(const float * lowerBound_, const float * upperBound_) :
		lowerBound({ lowerBound_[0], lowerBound_[1], lowerBound_[2] }), 
		upperBound({ upperBound_[0], upperBound_[1], upperBound_[2] })
	{
		// Validate that the upper bounds exceed the lower bounds.
		for (unsigned int i = 0; i < lowerBound.size(); i++)
		{
			// Validate the bound.
			if (lowerBound[i] > upperBound[i])
			{
				throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
			}
		}

		surfaceArea = computeSurfaceArea();
		centre = computeCentre();
	}

	template<unsigned int DIM, typename SCALAR>
	AABB<DIM, SCALAR>::AABB(const float * origin, float radius)
		:centre({ origin[0],origin[1],origin[2] })
	{
		for (int c = 0; c < DIM; ++c)
		{
			lowerBound[c] = origin[c] - radius;
			upperBound[c] = origin[c] + radius;
		}

		surfaceArea = computeSurfaceArea();
	}

	template<unsigned int DIM, typename SCALAR>
	SCALAR AABB<DIM, SCALAR>::computeSurfaceArea() const
	{
		// Sum of "area" of all the sides.
		SCALAR sum = 0;

		// General formula for one side: hold one dimension constant
		// and multiply by all the other ones.
		for (unsigned int d1 = 0; d1 < lowerBound.size(); d1++)
		{
			// "Area" of current side.
			SCALAR product = 1;

			for (unsigned int d2 = 0; d2 < lowerBound.size(); d2++)
			{
				if (d1 == d2)
					continue;

				SCALAR dx = upperBound[d2] - lowerBound[d2];
				product *= dx;
			}

			// Update the sum.
			sum += product;
		}

		return 2.0 * sum;
	}

	template<unsigned int DIM, typename SCALAR>
	SCALAR AABB<DIM, SCALAR>::getSurfaceArea() const
	{
		return surfaceArea;
	}

	template<unsigned int DIM, typename SCALAR>
	void AABB<DIM, SCALAR>::merge(const AABB& aabb1, const AABB& aabb2)
	{
		for (unsigned int i = 0; i < lowerBound.size(); i++)
		{
			lowerBound[i] = std::min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
			upperBound[i] = std::max(aabb1.upperBound[i], aabb2.upperBound[i]);
		}

		surfaceArea = computeSurfaceArea();
		centre = computeCentre();
	}

	template<unsigned int DIM, typename SCALAR>
	bool AABB<DIM, SCALAR>::contains(const AABB& aabb) const
	{
		assert(aabb.lowerBound.size() == lowerBound.size());

		for (unsigned int i = 0; i < lowerBound.size(); i++)
		{
			if (aabb.lowerBound[i] < lowerBound[i]) return false;
			if (aabb.upperBound[i] > upperBound[i]) return false;
		}

		return true;
	}

	template<unsigned int DIM, typename SCALAR>
	bool AABB<DIM, SCALAR>::overlaps(const AABB& aabb, bool touchIsOverlap) const
	{
		assert(aabb.lowerBound.size() == lowerBound.size());

		bool rv = true;

		if (touchIsOverlap)
		{
			for (unsigned int i = 0; i < lowerBound.size(); ++i)
			{
				if (aabb.upperBound[i] < lowerBound[i] || aabb.lowerBound[i] > upperBound[i])
				{
					rv = false;
					break;
				}
			}
		}
		else
		{
			for (unsigned int i = 0; i < lowerBound.size(); ++i)
			{
				if (aabb.upperBound[i] <= lowerBound[i] || aabb.lowerBound[i] >= upperBound[i])
				{
					rv = false;
					break;
				}
			}
		}

		return rv;
	}


	template<unsigned int DIM, typename SCALAR>
	typename AABB<DIM, SCALAR>::Position_t AABB<DIM, SCALAR>::computeCentre()
	{
		Position_t position;

		for (unsigned int i = 0; i < position.size(); i++)
			position[i] = 0.5 * (lowerBound[i] + upperBound[i]);

		return position;
	}

	template<unsigned int DIM, typename SCALAR>
	Node<DIM, SCALAR>::Node()
	{
	}

	template<unsigned int DIM, typename SCALAR>
	bool Node<DIM, SCALAR>::isLeaf() const
	{
		return (left == NULL_NODE);
	}

	template<unsigned int DIM, typename SCALAR>
	Tree<DIM, SCALAR>::Tree(
		SCALAR skinThickness_,
		unsigned int nParticles,
		bool touchIsOverlap_) :
		isPeriodic(false), 
		skinThickness(skinThickness_),
		touchIsOverlap(touchIsOverlap_)
	{
		// Initialise the periodicity vector.
		std::fill(periodicity.begin(), periodicity.end(), false);

		// Initialise the tree.
		root = NULL_NODE;
		nodeCount = 0;
		nodeCapacity = nParticles;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (unsigned int i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		// Assign the index of the first free node.
		freeList = 0;
	}

	template<unsigned int DIM, typename SCALAR>
	Tree<DIM, SCALAR>::Tree(
		SCALAR skinThickness_,
		const std::array<bool,DIM>& periodicity_,
		const float * boxSize_,
		unsigned int nParticles,
		bool touchIsOverlap_) :
		skinThickness(skinThickness_),
		periodicity(periodicity_),
		boxSize({boxSize_[0],boxSize_[1],boxSize_[2]}),
		touchIsOverlap(touchIsOverlap_)
	{
		// Initialise the tree.
		root = NULL_NODE;
		touchIsOverlap = true;
		nodeCount = 0;
		nodeCapacity = nParticles;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (unsigned int i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		// Assign the index of the first free node.
		freeList = 0;

		// Check periodicity.
		isPeriodic = false;
		for (unsigned int i = 0; i < DIM; i++)
		{
			posMinImage[i] = 0.5*boxSize[i];
			negMinImage[i] = -0.5*boxSize[i];

			if (periodicity[i])
				isPeriodic = true;
		}
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::setPeriodicity(const std::array<bool,DIM>& periodicity_)
	{
		periodicity = periodicity_;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::setBoxSize(const float * boxSize_)
	{
		boxSize = { boxSize_[0],boxSize_[1],boxSize_[2] };
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::allocateNode()
	{
		// Exand the node pool as needed.
		if (freeList == NULL_NODE)
		{
			assert(nodeCount == nodeCapacity);

			// The free list is empty. Rebuild a bigger pool.
			nodeCapacity *= 2;
			nodes.resize(nodeCapacity);

			// Build a linked list for the list of free nodes.
			for (unsigned int i = nodeCount; i < nodeCapacity - 1; i++)
			{
				nodes[i].next = i + 1;
				nodes[i].height = -1;
			}
			nodes[nodeCapacity - 1].next = NULL_NODE;
			nodes[nodeCapacity - 1].height = -1;

			// Assign the index of the first free node.
			freeList = nodeCount;
		}

		// Peel a node off the free list.
		unsigned int node = freeList;
		freeList = nodes[node].next;
		nodes[node].parent = NULL_NODE;
		nodes[node].left = NULL_NODE;
		nodes[node].right = NULL_NODE;
		nodes[node].height = 0;
		nodeCount++;

		return node;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::freeNode(unsigned int node)
	{
		assert(node < nodeCapacity);
		assert(0 < nodeCount);

		nodes[node].next = freeList;
		nodes[node].height = -1;
		freeList = node;
		nodeCount--;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::insertParticle(unsigned int particle, const float * position, SCALAR radius)
	{
		// Make sure the particle doesn't already exist.
		if (particleMap.count(particle) != 0)
		{
			throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
		}

		// Allocate a new node for the particle.
		unsigned int node = allocateNode();

		// AABB size in each dimension.
		std::array<SCALAR,DIM> size;

		// Compute the AABB limits.
		for (unsigned int i = 0; i < DIM; i++)
		{
			nodes[node].aabb.lowerBound[i] = position[i] - radius;
			nodes[node].aabb.upperBound[i] = position[i] + radius;
			size[i] = nodes[node].aabb.upperBound[i] - nodes[node].aabb.lowerBound[i];
		}

		// Fatten the AABB.
		for (unsigned int i = 0; i < DIM; i++)
		{
			nodes[node].aabb.lowerBound[i] -= skinThickness * size[i];
			nodes[node].aabb.upperBound[i] += skinThickness * size[i];
		}
		nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
		nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

		// Zero the height.
		nodes[node].height = 0;

		// Insert a new leaf into the tree.
		insertLeaf(node);

		// Add the new particle to the map.
		particleMap.insert(std::map<unsigned int, unsigned int>::value_type(particle, node));

		// Store the particle index.
		nodes[node].particle = particle;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::insertParticle(unsigned int particle, const float * lowerBound, const float *upperBound)
	{
		// Make sure the particle doesn't already exist.
		if (particleMap.count(particle) != 0)
		{
			throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
		}

		// Allocate a new node for the particle.
		unsigned int node = allocateNode();

		// AABB size in each dimension.
		Position_t size;

		// Compute the AABB limits.
		for (unsigned int i = 0; i < DIM; i++)
		{
			// Validate the bound.
			if (lowerBound[i] > upperBound[i])
			{
				throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
			}

			nodes[node].aabb.lowerBound[i] = lowerBound[i];
			nodes[node].aabb.upperBound[i] = upperBound[i];
			size[i] = upperBound[i] - lowerBound[i];
		}

		// Fatten the AABB.
		for (unsigned int i = 0; i < DIM; i++)
		{
			nodes[node].aabb.lowerBound[i] -= skinThickness * size[i];
			nodes[node].aabb.upperBound[i] += skinThickness * size[i];
		}
		nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
		nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

		// Zero the height.
		nodes[node].height = 0;

		// Insert a new leaf into the tree.
		insertLeaf(node);

		// Add the new particle to the map.
		particleMap.insert(std::map<unsigned int, unsigned int>::value_type(particle, node));

		// Store the particle index.
		nodes[node].particle = particle;
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::nParticles()
	{
		return particleMap.size();
	}

	template<unsigned int DIM, typename SCALAR>
	bool Tree<DIM, SCALAR>::containsParticle(unsigned int index)
	{
		return particleMap.find(index) != particleMap.end();
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::removeParticle(unsigned int particle)
	{
		// Map iterator.
		std::map<unsigned int, unsigned int>::iterator it;

		// Find the particle.
		it = particleMap.find(particle);

		// The particle doesn't exist.
		assert(it != particleMap.end());

		// Extract the node index.
		unsigned int node = it->second;

		// Erase the particle from the map.
		particleMap.erase(it);

		assert(node < nodeCapacity);
		assert(nodes[node].isLeaf());

		removeLeaf(node);
		freeNode(node);
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::removeAll()
	{
		// Initialise the tree.
		root = NULL_NODE;
		nodeCount = 0;
		nodeCapacity = 16;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (unsigned int i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		particleMap.clear();

		freeList = 0;
	}

	template<unsigned int DIM, typename SCALAR>
	bool Tree<DIM, SCALAR>::updateParticle(unsigned int particle, const float * position, SCALAR radius,
		bool alwaysReinsert)
	{
		// AABB bounds vectors.
		Position_t lowerBound;
		Position_t upperBound;

		// Compute the AABB limits.
		for (unsigned int i = 0; i < DIM; i++)
		{
			lowerBound[i] = position[i] - radius;
			upperBound[i] = position[i] + radius;
		}

		// Update the particle.
		return updateParticle(particle, &lowerBound[0], &upperBound[0], alwaysReinsert);
	}

	template<unsigned int DIM, typename SCALAR>
	bool Tree<DIM, SCALAR>::updateParticle(
		unsigned int particle, 
		const float * lowerBound,
		const float * upperBound, 
		bool alwaysReinsert)
	{
		// Map iterator.
		std::map<unsigned int, unsigned int>::iterator it;

		// Find the particle.
		it = particleMap.find(particle);

		// The particle doesn't exist.
		if (it == particleMap.end())
		{
			throw std::invalid_argument("[ERROR]: Invalid particle index!");
		}

		// Extract the node index.
		unsigned int node = it->second;

		assert(node < nodeCapacity);
		assert(nodes[node].isLeaf());

		// AABB size in each dimension.
		Position_t size;

		// Compute the AABB limits.
		for (unsigned int i = 0; i < DIM; i++)
		{
			// Validate the bound.
			if (lowerBound[i] > upperBound[i])
			{
				throw std::invalid_argument("[ERROR]: AABB lower bound is greater than the upper bound!");
			}

			size[i] = upperBound[i] - lowerBound[i];
		}

		// Create the new AABB.
		AABB_t aabb(lowerBound, upperBound);

		// No need to update if the particle is still within its fattened AABB.
		if (!alwaysReinsert && nodes[node].aabb.contains(aabb)) return false;

		// Remove the current leaf.
		removeLeaf(node);

		// Fatten the new AABB.
		for (unsigned int i = 0; i < DIM; i++)
		{
			aabb.lowerBound[i] -= skinThickness * size[i];
			aabb.upperBound[i] += skinThickness * size[i];
		}

		// Assign the new AABB.
		nodes[node].aabb = aabb;

		// Update the surface area and centroid.
		nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();
		nodes[node].aabb.centre = nodes[node].aabb.computeCentre();

		// Insert a new leaf node.
		insertLeaf(node);

		return true;
	}

	template<unsigned int DIM, typename SCALAR>
	std::vector<unsigned int> Tree<DIM, SCALAR>::query(unsigned int particle)
	{
		// Make sure that this is a valid particle.
		if (particleMap.count(particle) == 0)
		{
			throw std::invalid_argument("[ERROR]: Invalid particle index!");
		}

		// Test overlap of particle AABB against all other particles.
		return query(particle, nodes[particleMap.find(particle)->second].aabb);
	}


	template<unsigned int DIM, typename SCALAR>
	std::vector<unsigned int> Tree<DIM, SCALAR>::query(unsigned int particleSelf, const AABB_t& aabb)
	{
		std::vector<unsigned int> particles;

		query(aabb, [&](unsigned int particle) {
			if (particle != particleSelf)
				particles.push_back(particle);

			return true; // continue the query
		});

		return particles;
	}

	template<unsigned int DIM, typename SCALAR>
	std::vector<unsigned int> Tree<DIM, SCALAR>::query(const typename Tree<DIM,SCALAR>::AABB_t& aabb)
	{
		// Make sure the tree isn't empty.
		if (particleMap.size() == 0)
		{
			return std::vector<unsigned int>();
		}

		// Test overlap of AABB against all particles.
		return query(std::numeric_limits<unsigned int>::max(), aabb);
	}

	template<unsigned int DIM, typename SCALAR>
	const typename Tree<DIM,SCALAR>::AABB_t& Tree<DIM, SCALAR>::getAABB(unsigned int particle)
	{
		return nodes[particleMap[particle]].aabb;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::insertLeaf(unsigned int leaf)
	{
		if (root == NULL_NODE)
		{
			root = leaf;
			nodes[root].parent = NULL_NODE;
			return;
		}

		// Find the best sibling for the node.

		AABB_t leafAABB = nodes[leaf].aabb;
		unsigned int index = root;

		while (!nodes[index].isLeaf())
		{
			// Extract the children of the node.
			unsigned int left = nodes[index].left;
			unsigned int right = nodes[index].right;

			SCALAR surfaceArea = nodes[index].aabb.getSurfaceArea();

			AABB_t combinedAABB;
			combinedAABB.merge(nodes[index].aabb, leafAABB);
			SCALAR combinedSurfaceArea = combinedAABB.getSurfaceArea();

			// Cost of creating a new parent for this node and the new leaf.
			SCALAR cost = 2.0 * combinedSurfaceArea;

			// Minimum cost of pushing the leaf further down the tree.
			SCALAR inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

			// Cost of descending to the left.
			SCALAR costLeft;
			if (nodes[left].isLeaf())
			{
				AABB_t aabb;
				aabb.merge(leafAABB, nodes[left].aabb);
				costLeft = aabb.getSurfaceArea() + inheritanceCost;
			}
			else
			{
				AABB_t aabb;
				aabb.merge(leafAABB, nodes[left].aabb);
				SCALAR oldArea = nodes[left].aabb.getSurfaceArea();
				SCALAR newArea = aabb.getSurfaceArea();
				costLeft = (newArea - oldArea) + inheritanceCost;
			}

			// Cost of descending to the right.
			SCALAR costRight;
			if (nodes[right].isLeaf())
			{
				AABB_t aabb;
				aabb.merge(leafAABB, nodes[right].aabb);
				costRight = aabb.getSurfaceArea() + inheritanceCost;
			}
			else
			{
				AABB_t aabb;
				aabb.merge(leafAABB, nodes[right].aabb);
				SCALAR oldArea = nodes[right].aabb.getSurfaceArea();
				SCALAR newArea = aabb.getSurfaceArea();
				costRight = (newArea - oldArea) + inheritanceCost;
			}

			// Descend according to the minimum cost.
			if ((cost < costLeft) && (cost < costRight)) break;

			// Descend.
			if (costLeft < costRight) index = left;
			else                      index = right;
		}

		unsigned int sibling = index;

		// Create a new parent.
		unsigned int oldParent = nodes[sibling].parent;
		unsigned int newParent = allocateNode();
		nodes[newParent].parent = oldParent;
		nodes[newParent].aabb.merge(leafAABB, nodes[sibling].aabb);
		nodes[newParent].height = nodes[sibling].height + 1;

		// The sibling was not the root.
		if (oldParent != NULL_NODE)
		{
			if (nodes[oldParent].left == sibling) nodes[oldParent].left = newParent;
			else                                  nodes[oldParent].right = newParent;

			nodes[newParent].left = sibling;
			nodes[newParent].right = leaf;
			nodes[sibling].parent = newParent;
			nodes[leaf].parent = newParent;
		}
		// The sibling was the root.
		else
		{
			nodes[newParent].left = sibling;
			nodes[newParent].right = leaf;
			nodes[sibling].parent = newParent;
			nodes[leaf].parent = newParent;
			root = newParent;
		}

		// Walk back up the tree fixing heights and AABBs.
		index = nodes[leaf].parent;
		while (index != NULL_NODE)
		{
			index = balance(index);

			unsigned int left = nodes[index].left;
			unsigned int right = nodes[index].right;

			assert(left != NULL_NODE);
			assert(right != NULL_NODE);

			nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);
			nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);

			index = nodes[index].parent;
		}
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::removeLeaf(unsigned int leaf)
	{
		if (leaf == root)
		{
			root = NULL_NODE;
			return;
		}

		unsigned int parent = nodes[leaf].parent;
		unsigned int grandParent = nodes[parent].parent;
		unsigned int sibling;

		if (nodes[parent].left == leaf) sibling = nodes[parent].right;
		else                            sibling = nodes[parent].left;

		// Destroy the parent and connect the sibling to the grandparent.
		if (grandParent != NULL_NODE)
		{
			if (nodes[grandParent].left == parent) nodes[grandParent].left = sibling;
			else                                   nodes[grandParent].right = sibling;

			nodes[sibling].parent = grandParent;
			freeNode(parent);

			// Adjust ancestor bounds.
			unsigned int index = grandParent;
			while (index != NULL_NODE)
			{
				index = balance(index);

				unsigned int left = nodes[index].left;
				unsigned int right = nodes[index].right;

				nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);
				nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);

				index = nodes[index].parent;
			}
		}
		else
		{
			root = sibling;
			nodes[sibling].parent = NULL_NODE;
			freeNode(parent);
		}
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::balance(unsigned int node)
	{
		assert(node != NULL_NODE);

		if (nodes[node].isLeaf() || (nodes[node].height < 2))
			return node;

		unsigned int left = nodes[node].left;
		unsigned int right = nodes[node].right;

		assert(left < nodeCapacity);
		assert(right < nodeCapacity);

		int currentBalance = nodes[right].height - nodes[left].height;

		// Rotate right branch up.
		if (currentBalance > 1)
		{
			unsigned int rightLeft = nodes[right].left;
			unsigned int rightRight = nodes[right].right;

			assert(rightLeft < nodeCapacity);
			assert(rightRight < nodeCapacity);

			// Swap node and its right-hand child.
			nodes[right].left = node;
			nodes[right].parent = nodes[node].parent;
			nodes[node].parent = right;

			// The node's old parent should now point to its right-hand child.
			if (nodes[right].parent != NULL_NODE)
			{
				if (nodes[nodes[right].parent].left == node) nodes[nodes[right].parent].left = right;
				else
				{
					assert(nodes[nodes[right].parent].right == node);
					nodes[nodes[right].parent].right = right;
				}
			}
			else root = right;

			// Rotate.
			if (nodes[rightLeft].height > nodes[rightRight].height)
			{
				nodes[right].right = rightLeft;
				nodes[node].right = rightRight;
				nodes[rightRight].parent = node;
				nodes[node].aabb.merge(nodes[left].aabb, nodes[rightRight].aabb);
				nodes[right].aabb.merge(nodes[node].aabb, nodes[rightLeft].aabb);

				nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightRight].height);
				nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightLeft].height);
			}
			else
			{
				nodes[right].right = rightRight;
				nodes[node].right = rightLeft;
				nodes[rightLeft].parent = node;
				nodes[node].aabb.merge(nodes[left].aabb, nodes[rightLeft].aabb);
				nodes[right].aabb.merge(nodes[node].aabb, nodes[rightRight].aabb);

				nodes[node].height = 1 + std::max(nodes[left].height, nodes[rightLeft].height);
				nodes[right].height = 1 + std::max(nodes[node].height, nodes[rightRight].height);
			}

			return right;
		}

		// Rotate left branch up.
		if (currentBalance < -1)
		{
			unsigned int leftLeft = nodes[left].left;
			unsigned int leftRight = nodes[left].right;

			assert(leftLeft < nodeCapacity);
			assert(leftRight < nodeCapacity);

			// Swap node and its left-hand child.
			nodes[left].left = node;
			nodes[left].parent = nodes[node].parent;
			nodes[node].parent = left;

			// The node's old parent should now point to its left-hand child.
			if (nodes[left].parent != NULL_NODE)
			{
				if (nodes[nodes[left].parent].left == node) nodes[nodes[left].parent].left = left;
				else
				{
					assert(nodes[nodes[left].parent].right == node);
					nodes[nodes[left].parent].right = left;
				}
			}
			else root = left;

			// Rotate.
			if (nodes[leftLeft].height > nodes[leftRight].height)
			{
				nodes[left].right = leftLeft;
				nodes[node].left = leftRight;
				nodes[leftRight].parent = node;
				nodes[node].aabb.merge(nodes[right].aabb, nodes[leftRight].aabb);
				nodes[left].aabb.merge(nodes[node].aabb, nodes[leftLeft].aabb);

				nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftRight].height);
				nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftLeft].height);
			}
			else
			{
				nodes[left].right = leftRight;
				nodes[node].left = leftLeft;
				nodes[leftLeft].parent = node;
				nodes[node].aabb.merge(nodes[right].aabb, nodes[leftLeft].aabb);
				nodes[left].aabb.merge(nodes[node].aabb, nodes[leftRight].aabb);

				nodes[node].height = 1 + std::max(nodes[right].height, nodes[leftLeft].height);
				nodes[left].height = 1 + std::max(nodes[node].height, nodes[leftRight].height);
			}

			return left;
		}

		return node;
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::computeHeight() const
	{
		return computeHeight(root);
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::computeHeight(unsigned int node) const
	{
		assert(node < nodeCapacity);

		if (nodes[node].isLeaf()) return 0;

		unsigned int height1 = computeHeight(nodes[node].left);
		unsigned int height2 = computeHeight(nodes[node].right);

		return 1 + std::max(height1, height2);
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::getHeight() const
	{
		if (root == NULL_NODE) return 0;
		return nodes[root].height;
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::getNodeCount() const
	{
		return nodeCount;
	}

	template<unsigned int DIM, typename SCALAR>
	unsigned int Tree<DIM, SCALAR>::computeMaximumBalance() const
	{
		unsigned int maxBalance = 0;
		for (unsigned int i = 0; i < nodeCapacity; i++)
		{
			if (nodes[i].height <= 1)
				continue;

			assert(nodes[i].isLeaf() == false);

			unsigned int balance = std::abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
			maxBalance = std::max(maxBalance, balance);
		}

		return maxBalance;
	}

	template<unsigned int DIM, typename SCALAR>
	SCALAR Tree<DIM, SCALAR>::computeSurfaceAreaRatio() const
	{
		if (root == NULL_NODE) return 0.0;

		SCALAR rootArea = nodes[root].aabb.computeSurfaceArea();
		SCALAR totalArea = 0.0;

		for (unsigned int i = 0; i < nodeCapacity; i++)
		{
			if (nodes[i].height < 0) continue;

			totalArea += nodes[i].aabb.computeSurfaceArea();
		}

		return totalArea / rootArea;
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::validate() const
	{
#ifndef NDEBUG
		validateStructure(root);
		validateMetrics(root);

		unsigned int freeCount = 0;
		unsigned int freeIndex = freeList;

		while (freeIndex != NULL_NODE)
		{
			assert(freeIndex < nodeCapacity);
			freeIndex = nodes[freeIndex].next;
			freeCount++;
		}

		assert(getHeight() == computeHeight());
		assert((nodeCount + freeCount) == nodeCapacity);
#endif
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::rebuild()
	{
		std::vector<unsigned int> nodeIndices(nodeCount);
		unsigned int count = 0;

		for (unsigned int i = 0; i < nodeCapacity; i++)
		{
			// Free node.
			if (nodes[i].height < 0) continue;

			if (nodes[i].isLeaf())
			{
				nodes[i].parent = NULL_NODE;
				nodeIndices[count] = i;
				count++;
			}
			else freeNode(i);
		}

		while (count > 1)
		{
			SCALAR minCost = std::numeric_limits<SCALAR>::max();
			int iMin = -1, jMin = -1;

			for (unsigned int i = 0; i < count; i++)
			{
				AABB_t aabbi = nodes[nodeIndices[i]].aabb;

				for (unsigned int j = i + 1; j < count; j++)
				{
					AABB_t aabbj = nodes[nodeIndices[j]].aabb;
					AABB_t aabb;
					aabb.merge(aabbi, aabbj);
					SCALAR cost = aabb.getSurfaceArea();

					if (cost < minCost)
					{
						iMin = i;
						jMin = j;
						minCost = cost;
					}
				}
			}

			unsigned int index1 = nodeIndices[iMin];
			unsigned int index2 = nodeIndices[jMin];

			unsigned int parent = allocateNode();
			nodes[parent].left = index1;
			nodes[parent].right = index2;
			nodes[parent].height = 1 + std::max(nodes[index1].height, nodes[index2].height);
			nodes[parent].aabb.merge(nodes[index1].aabb, nodes[index2].aabb);
			nodes[parent].parent = NULL_NODE;

			nodes[index1].parent = parent;
			nodes[index2].parent = parent;

			nodeIndices[jMin] = nodeIndices[count - 1];
			nodeIndices[iMin] = parent;
			count--;
		}

		root = nodeIndices[0];

		validate();
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::validateStructure(unsigned int node) const
	{
		if (node == NULL_NODE) return;

		if (node == root) assert(nodes[node].parent == NULL_NODE);

		unsigned int left = nodes[node].left;
		unsigned int right = nodes[node].right;

		if (nodes[node].isLeaf())
		{
			assert(left == NULL_NODE);
			assert(right == NULL_NODE);
			assert(nodes[node].height == 0);
			return;
		}

		assert(left < nodeCapacity);
		assert(right < nodeCapacity);

		assert(nodes[left].parent == node);
		assert(nodes[right].parent == node);

		validateStructure(left);
		validateStructure(right);
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::validateMetrics(unsigned int node) const
	{
		if (node == NULL_NODE) return;

		unsigned int left = nodes[node].left;
		unsigned int right = nodes[node].right;

		if (nodes[node].isLeaf())
		{
			assert(left == NULL_NODE);
			assert(right == NULL_NODE);
			assert(nodes[node].height == 0);
			return;
		}

		assert(left < nodeCapacity);
		assert(right < nodeCapacity);

		int height1 = nodes[left].height;
		int height2 = nodes[right].height;
		int height = 1 + std::max(height1, height2);
		(void)height; // Unused variable in Release build
		assert(nodes[node].height == height);

		AABB_t aabb;
		aabb.merge(nodes[left].aabb, nodes[right].aabb);

		for (unsigned int i = 0; i < DIM; i++)
		{
			assert(aabb.lowerBound[i] == nodes[node].aabb.lowerBound[i]);
			assert(aabb.upperBound[i] == nodes[node].aabb.upperBound[i]);
		}

		validateMetrics(left);
		validateMetrics(right);
	}

	template<unsigned int DIM, typename SCALAR>
	void Tree<DIM, SCALAR>::periodicBoundaries(float * position)
	{
		for (unsigned int i = 0; i < DIM; i++)
		{
			if (position[i] < 0)
			{
				position[i] += boxSize[i];
			}
			else
			{
				if (position[i] >= boxSize[i])
				{
					position[i] -= boxSize[i];
				}
			}
		}
	}

	template<unsigned int DIM, typename SCALAR>
	bool Tree<DIM, SCALAR>::minimumImage(typename Tree<DIM, SCALAR>::Position_t& separation, typename Tree<DIM, SCALAR>::Position_t& shift)
	{
		bool isShifted = false;

		for (unsigned int i = 0; i < DIM; i++)
		{
			if (separation[i] < negMinImage[i])
			{
				separation[i] += periodicity[i] * boxSize[i];
				shift[i] = periodicity[i] * boxSize[i];
				isShifted = true;
			}
			else
			{
				if (separation[i] >= posMinImage[i])
				{
					separation[i] -= float(periodicity[i]) * boxSize[i];
					shift[i] = -float(periodicity[i]) * boxSize[i];
					isShifted = true;
				}
			}
		}

		return isShifted;
	}



	template class AABB<3, float>;
	template struct Node<3, float>;
	template class Tree<3, float>;
}


#endif /* _AABB_H */