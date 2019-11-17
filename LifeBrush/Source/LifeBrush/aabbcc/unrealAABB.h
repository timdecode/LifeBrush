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

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>
#include <array>

namespace unrealAABB
{
	/// Null node flag.
	const uint32_t NULL_NODE = 0xffffffff;

	/*! \brief Maps particle indices to node indices.

	The assumption is that the particle indices are fairly compact between 0 and the number of particles.
	Behind the scenes, the map is a vector. So, if you have particle indices in the billions, you are going
	to pay the price in memory for the sparsity of the indices.
	*/
	struct ParticleMap
	{
		std::vector<int> particleMap;

		size_t _size = 0;

		void clear() {
			particleMap.clear();
		}

		void erase(uint32_t particleIndex)
		{
			if (particleIndex < particleMap.size() && particleMap[particleIndex] >= 0)
			{
				particleMap[particleIndex] = -1;

				_size--;
			}
		}

		bool contains(uint32_t particleIndex)
		{
			if (particleIndex >= particleMap.size())
				return false;
			else
				return particleMap[particleIndex] >= 0;
		}

		uint32_t find(uint32_t particleIndex)
		{
			return particleMap[particleIndex];
		}

		void insert(uint32_t particleIndex, uint32_t nodeIndex)
		{
			if (particleIndex >= particleMap.size())
				particleMap.resize(particleIndex + 1, -1);

			if (particleMap[particleIndex] < 0)
				_size++;

			particleMap[particleIndex] = nodeIndex;
		}

		size_t size()
		{
			return _size;
		}
	};

	struct ParticleMap_TMapBased
	{
		TMap<uint32_t, uint32_t> particleMap;

		void clear() {
			particleMap.Empty();
		}

		void erase(uint32_t particleIndex)
		{
			particleMap.Remove(particleIndex);
		}

		bool contains(uint32_t particleIndex)
		{
			return particleMap.Contains(particleIndex);
		}

		uint32_t find(uint32_t particleIndex)
		{
			return *particleMap.Find(particleIndex);
		}

		void insert(uint32_t particleIndex, uint32_t nodeIndex)
		{
			particleMap.Emplace(particleIndex, nodeIndex);
		}

		size_t size()
		{
			return particleMap.Num();
		}
	};

	/*! \brief The axis-aligned bounding box object.

		Axis-aligned bounding boxes (AABBs) store information for the minimum
		orthorhombic bounding-box for an object. Support is provided for
		dimensions >= 2. (In 2D the bounding box is either a rectangle,
		in 3D it is a rectangular prism.)

		Class member functions provide functionality for merging AABB objects
		and testing overlap with other AABBs.
	 */
	class AABB
	{
	public:
		/// Constructor.
		AABB();


		//! Constructor.
		/*! \param lowerBound_
				The lower bound in each dimension.

			\param upperBound_
				The upper bound in each dimension.
		 */
		AABB(FVector lowerBound, FVector upperBound);

		AABB(FVector origin, float radius);

		/// Compute the surface area of the box.
		float computeSurfaceArea() const;

		/// Get the surface area of the box.
		float getSurfaceArea() const;

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

		//! The position of the AABB centre.
		FVector centre() { return (lowerBound + upperBound) * 0.5f; }

		/// Lower bound of AABB in each dimension.
		FVector lowerBound;

		/// Upper bound of AABB in each dimension.
		FVector upperBound;



		/// The AABB's surface area.
		float surfaceArea;
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
	struct Node
	{
		/// Constructor.
		Node();

		/// The fattened axis-aligned bounding box.
		AABB aabb;

		/// Index of the parent node.
		uint32_t parent;

		/// Index of the next node.
		uint32_t next;

		/// Index of the left-hand child.
		uint32_t left;

		/// Index of the right-hand child.
		uint32_t right;

		/// Height of the node. This is 0 for a leaf and -1 for a free node.
		int height;

		/// The index of the particle that the node contains (leaf nodes only).
		uint32_t particle;

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
	class Tree
	{
	public:
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
		Tree(float skinThickness_ = 0.05,
			uint32_t nParticles = 16, bool touchIsOverlap = true);

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
		Tree(float, const float * boxSize,
			uint32_t nParticles = 16, bool touchIsOverlap = true);

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
		void insertParticle(uint32_t index, FVector position, float radius);

		//! Insert a particle into the tree (arbitrary shape with bounding box).
		/*! \param index
				The index of the particle.

			\param lowerBound
				The lower bound in each dimension.

			\param upperBound
				The upper bound in each dimension.
		 */
		void insertParticle(uint32_t index, FVector lowerBound, FVector upperBound);

		/// Return the number of particles in the tree.
		uint32_t nParticles();

		//! Remove a particle from the tree.
		/*! \param particle
				The particle index (particleMap will be used to map the node).
		 */
		void removeParticle(uint32_t index);

		bool containsParticle(uint32_t index);

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
		bool updateParticle(uint32_t index, const FVector position, float radius);

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
		bool updateParticle(uint32_t index, const FVector lowerBound, const FVector upperBound);
		bool updateParticleReinstert(uint32_t index, const FVector lowerBound, const FVector upperBound);

		//! Query the tree to find candidate interactions for a particle.
		/*! \param particle
				The particle index.

			\return particles
				A vector of particle indices.
		 */
		std::vector<uint32_t> query(uint32_t index);

		//! Query the tree to find candidate interactions for an AABB.
		/*! \param particle
				The particle index.

			\param aabb
				The AABB.

			\return particles
				A vector of particle indices.
		 */
		std::vector<uint32_t> query(uint32_t index, const AABB& aabbQuery);

		//! Query the tree to find candidate interactions for an AABB.
		/*! \param aabb
				The AABB.

			\return particles
				A vector of particle indices.
		 */
		std::vector<uint32_t> query(const AABB& aabbQuery);

		template<typename Callback>
		void query(const AABB& aabbQuery, Callback callback)
		{
			query_ref(aabbQuery, callback);
		}

		// \param callback A function that accepts a particle index (uint32_t) whose AABB overlaps
		// with the query AABB. It should return true if the query should continue, and false, if the
		// query should end.
		template<typename Callback>
		void query_ref(const AABB& aabbQuery, Callback& callback)
		{
			if (nParticles() == 0)
				return;

			std::vector<uint32_t> stack;

			stack.reserve(16);
			stack.push_back(root);

			while (stack.size() > 0)
			{
				uint32_t node = stack.back();
				stack.pop_back();

				// Copy the AABB.
				auto nodeAABB = nodes[node].aabb;

				if (node == NULL_NODE) continue;

				// Test for overlap between the AABBs.
				if (aabbQuery.overlaps(nodeAABB, touchIsOverlap))
				{
					// Check that we're at a leaf node.
					if (nodes[node].isLeaf())
					{
						if (!callback(nodes[node].particle))
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
		const AABB& getAABB(uint32_t index);

		//! Get the height of the tree.
		/*! \return
				The height of the binary tree.
		 */
		uint32_t getHeight() const;

		//! Get the number of nodes in the tree.
		/*! \return
				The number of nodes in the tree.
		 */
		uint32_t getNodeCount() const;

		//! Compute the maximum balancance of the tree.
		/*! \return
				The maximum difference between the height of two
				children of a node.
		 */
		uint32_t computeMaximumBalance() const;

		//! Compute the surface area ratio of the tree.
		/*! \return
				The ratio of the sum of the node surface area to the surface
				area of the root node.
		 */
		float computeSurfaceAreaRatio() const;

		/// Validate the tree.
		void validate() const;

		/// Rebuild an optimal tree.
		void rebuild();

	private:
		/// The index of the root node.
		uint32_t root;

		/// The dynamic tree.
		std::vector<Node> nodes;

		/// The current number of nodes in the tree.
		uint32_t nodeCount;

		/// The current node capacity.
		uint32_t nodeCapacity;

		/// The position of node at the top of the free list.
		uint32_t freeList;

		/// Whether the system is periodic along at least one axis.
		bool isPeriodic;

		/// The skin thickness of the fattened AABBs, as a fraction of the AABB base length.
		float skinThickness;

		/// The size of the system in each dimension.
		FVector boxSize;

		/// A map between particle and node indices.
		ParticleMap_TMapBased particleMap;

		/// Does touching count as overlapping in tree queries?
		bool touchIsOverlap;

		//! Allocate a new node.
		/*! \return
				The index of the allocated node.
		 */
		uint32_t allocateNode();

		//! Free an existing node.
		/*! \param node
				The index of the node to be freed.
		 */
		void freeNode(uint32_t);

		//! Insert a leaf into the tree.
		/*! \param leaf
				The index of the leaf node.
		 */
		void insertLeaf(uint32_t);

		//! Remove a leaf from the tree.
		/*! \param leaf
				The index of the leaf node.
		 */
		void removeLeaf(uint32_t);

		//! Balance the tree.
		/*! \param node
				The index of the node.
		 */
		uint32_t balance(uint32_t);

		//! Compute the height of the tree.
		/*! \return
				The height of the entire tree.
		 */
		uint32_t computeHeight() const;

		//! Compute the height of a sub-tree.
		/*! \param node
				The index of the root node.

			\return
				The height of the sub-tree.
		 */
		uint32_t computeHeight(uint32_t) const;

		//! Assert that the sub-tree has a valid structure.
		/*! \param node
				The index of the root node.
		 */
		void validateStructure(uint32_t) const;

		//! Assert that the sub-tree has valid metrics.
		/*! \param node
				The index of the root node.
		 */
		void validateMetrics(uint32_t) const;

	};
}


