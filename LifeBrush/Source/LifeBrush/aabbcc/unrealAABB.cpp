/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>
  Copyright (c) 2019 Timothy Davison 
     I ported this to Unreal and got rid of the crazy dynamic allocation all over the place!

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

#include "LifeBrush.h"

#include "unrealAABB.h"

namespace unrealAABB
{

	AABB::AABB()
	{

	}


	AABB::AABB(FVector lowerBound, FVector upperBound) :
		lowerBound(lowerBound),
		upperBound(upperBound)
	{
		// Validate that the upper bounds exceed the lower bounds.
		for (uint32_t i = 0; i < 3; i++)
		{
			// Validate the bound.
			assert(lowerBound[i] <= upperBound[i] && "[ERROR]: AABB lower bound is greater than the upper bound!");
		}

		surfaceArea = computeSurfaceArea();
	}

	AABB::AABB(const FVector origin, float radius)
	{
		lowerBound = origin - FVector(radius);
		upperBound = origin + FVector(radius);

		surfaceArea = computeSurfaceArea();
	}

	float AABB::computeSurfaceArea() const
	{
		// Sum of "area" of all the sides.
		float sum = 0;

		// General formula for one side: hold one dimension constant
		// and multiply by all the other ones.
		for (uint32_t d1 = 0; d1 < 3; d1++)
		{
			// "Area" of current side.
			float product = 1;

			for (uint32_t d2 = 0; d2 < 3; d2++)
			{
				if (d1 == d2)
					continue;

				float dx = upperBound[d2] - lowerBound[d2];
				product *= dx;
			}

			// Update the sum.
			sum += product;
		}

		return 2.0 * sum;
	}

	float AABB::getSurfaceArea() const
	{
		return surfaceArea;
	}

	void AABB::merge(const AABB& aabb1, const AABB& aabb2)
	{
		lowerBound = aabb1.lowerBound.ComponentMin(aabb2.lowerBound);
		upperBound = aabb1.upperBound.ComponentMax(aabb2.upperBound);

		surfaceArea = computeSurfaceArea();
	}

	bool AABB::contains(const AABB& aabb) const
	{
		const FVector& l = aabb.lowerBound;
		const FVector& u = aabb.upperBound;

		return l.X >= lowerBound.X && l.Y >= lowerBound.Y && l.Z >= lowerBound.Z &&
			u.X <= upperBound.X && u.Y <= upperBound.Y && u.Z <= upperBound.Z;
	}

	bool AABB::overlaps(const AABB& aabb, bool touchIsOverlap) const
	{
		assert(aabb.lowerBound.size() == lowerBound.size());

		bool rv = true;

		if (touchIsOverlap)
		{
			for (uint32_t i = 0; i < 3; ++i)
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
			for (uint32_t i = 0; i < 3; ++i)
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

	Node::Node()
	{
	}

	bool Node::isLeaf() const
	{
		return (left == NULL_NODE);
	}

	Tree::Tree(
		float skinThickness_,
		uint32_t nParticles,
		bool touchIsOverlap_) :
		skinThickness(skinThickness_),
		touchIsOverlap(touchIsOverlap_)
	{
		// Initialise the tree.
		root = NULL_NODE;
		nodeCount = 0;
		nodeCapacity = nParticles;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (uint32_t i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		// Assign the index of the first free node.
		freeList = 0;
	}

	Tree::Tree(
		float skinThickness_,
		const float * boxSize_,
		uint32_t nParticles,
		bool touchIsOverlap_) :
		skinThickness(skinThickness_),
		boxSize(boxSize_[0], boxSize_[1], boxSize_[2]),
		touchIsOverlap(touchIsOverlap_)
	{
		// Initialise the tree.
		root = NULL_NODE;
		touchIsOverlap = true;
		nodeCount = 0;
		nodeCapacity = nParticles;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (uint32_t i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		// Assign the index of the first free node.
		freeList = 0;

	}

	void Tree::setBoxSize(const float * boxSize_)
	{
		boxSize = { boxSize_[0],boxSize_[1],boxSize_[2] };
	}

	uint32_t Tree::allocateNode()
	{
		// Exand the node pool as needed.
		if (freeList == NULL_NODE)
		{
			assert(nodeCount == nodeCapacity);

			// The free list is empty. Rebuild a bigger pool.
			nodeCapacity *= 2;
			nodes.resize(nodeCapacity);

			// Build a linked list for the list of free nodes.
			for (uint32_t i = nodeCount; i < nodeCapacity - 1; i++)
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
		uint32_t node = freeList;
		freeList = nodes[node].next;
		nodes[node].parent = NULL_NODE;
		nodes[node].left = NULL_NODE;
		nodes[node].right = NULL_NODE;
		nodes[node].height = 0;
		nodeCount++;

		return node;
	}

	void Tree::freeNode(uint32_t node)
	{
		assert(node < nodeCapacity);
		assert(0 < nodeCount);

		nodes[node].next = freeList;
		nodes[node].height = -1;
		freeList = node;
		nodeCount--;
	}

	void Tree::insertParticle(uint32_t particle, FVector position, float radius)
	{
		// Make sure the particle doesn't already exist.
		assert(!particleMap.contains(particle) && "The particle already exists.");

		// Allocate a new node for the particle.
		uint32_t node = allocateNode();

		Node& theNode = nodes[node];


		// Compute the AABB limits.
		theNode.aabb.lowerBound = position - radius;
		theNode.aabb.upperBound = position + radius;

		// AABB size in each dimension.
		FVector size = nodes[node].aabb.upperBound - nodes[node].aabb.lowerBound;

		// Fatten the AABB.
		theNode.aabb.lowerBound -= skinThickness * size;
		theNode.aabb.upperBound += skinThickness * size;

		theNode.aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();

		// Zero the height.
		theNode.height = 0;
		
		// Store the particle index.
		theNode.particle = particle;

		// Insert a new leaf into the tree.
		insertLeaf(node);

		// Add the new particle to the map.
		particleMap.insert(particle, node);
	}

	void Tree::insertParticle(uint32_t particle, FVector lowerBound, FVector upperBound)
	{
		// Make sure the particle doesn't already exist.
		assert(!particleMap.contains(particle) && "The particle already exists.");

		// Allocate a new node for the particle.
		uint32_t node = allocateNode();

		Node& theNode = nodes[node];



		// Compute the AABB limits.
		for (uint32_t i = 0; i < 3; i++)
		{
			// Validate the bound.
			assert(lowerBound[i] <= upperBound[i] && "[ERROR]: AABB lower bound is greater than the upper bound!");
		}

		theNode.aabb.lowerBound = lowerBound;
		theNode.aabb.upperBound = upperBound;

		// AABB size in each dimension.
		FVector size = upperBound - lowerBound;


		// Fatten the AABB.
		theNode.aabb.lowerBound -= skinThickness * size;
		theNode.aabb.upperBound += skinThickness * size;

		theNode.aabb.surfaceArea = theNode.aabb.computeSurfaceArea();

		// Zero the height.
		theNode.height = 0;

		// Store the particle index.
		theNode.particle = particle;

		// Insert a new leaf into the tree.
		insertLeaf(node);

		// Add the new particle to the map.
		particleMap.insert(particle, node);
	}

	uint32_t Tree::nParticles()
	{
		return particleMap.size();
	}

	bool Tree::containsParticle(uint32_t index)
	{
		return particleMap.contains(index);
	}

	void Tree::removeParticle(uint32_t particle)
	{
		assert(!particleMap.contains(particle) && "The particle was not in the tree.");

		// Extract the node index.
		uint32_t node = particleMap.find(particle);

		// Erase the particle from the map.
		particleMap.erase(particle);

		assert(node < nodeCapacity);
		assert(nodes[node].isLeaf());

		removeLeaf(node);
		freeNode(node);
	}

	void Tree::removeAll()
	{
		// Initialise the tree.
		root = NULL_NODE;
		nodeCount = 0;
		nodeCapacity = 16;
		nodes.resize(nodeCapacity);

		// Build a linked list for the list of free nodes.
		for (uint32_t i = 0; i < nodeCapacity - 1; i++)
		{
			nodes[i].next = i + 1;
			nodes[i].height = -1;
		}
		nodes[nodeCapacity - 1].next = NULL_NODE;
		nodes[nodeCapacity - 1].height = -1;

		particleMap.clear();

		freeList = 0;
	}

	bool Tree::updateParticle(uint32_t particle, const FVector position, float radius)
	{
		// AABB bounds vectors.
		FVector lowerBound = position - radius;
		FVector upperBound = position + radius;

		// Update the particle.
		return updateParticle(particle, lowerBound, upperBound);
	}


	bool Tree::updateParticle(
		uint32_t particle,
		const FVector lowerBound,
		const FVector upperBound)
	{
		assert(particleMap.contains(particle) && "The particle was not in the tree.");

		// Extract the node index.
		uint32_t node = particleMap.find(particle);

		assert(node < nodeCapacity);
		assert(nodes[node].isLeaf());

		// AABB size in each dimension.
		FVector size = upperBound - lowerBound;

		// Compute the AABB limits.
		for (uint32_t i = 0; i < 3; i++)
		{
			// Validate the bound.
			assert(lowerBound[i] <= upperBound[i] && "[ERROR]: AABB lower bound is greater than the upper bound!");
		}


		// Create the new AABB.
		AABB aabb(lowerBound, upperBound);

		// No need to update if the particle is still within its fattened AABB.
		if (nodes[node].aabb.contains(aabb)) return false;

		// Remove the current leaf.
		removeLeaf(node);

		// Fatten the new AABB.
		aabb.lowerBound -= skinThickness * size;
		aabb.upperBound += skinThickness * size;

		// Assign the new AABB.
		nodes[node].aabb = aabb;

		// Update the surface area and centroid.
		nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();

		// Insert a new leaf node.
		insertLeaf(node);

		return true;
	}
	bool Tree::updateParticleReinstert(
		uint32_t particle,
		const FVector lowerBound,
		const FVector upperBound)
	{
		assert(particleMap.contains(particle) && "The particle was not in the tree.");

		// Extract the node index.
		uint32_t node = particleMap.find(particle);

		assert(node < nodeCapacity);
		assert(nodes[node].isLeaf());

		// AABB size in each dimension.
		FVector size = upperBound - lowerBound;

		// Compute the AABB limits.
		for (uint32_t i = 0; i < 3; i++)
		{
			// Validate the bound.
			assert(lowerBound[i] <= upperBound[i] && "[ERROR]: AABB lower bound is greater than the upper bound!");
		}

		// Create the new AABB.
		AABB aabb(lowerBound, upperBound);

		// Remove the current leaf.
		removeLeaf(node);

		// Fatten the new AABB.
		aabb.lowerBound -= skinThickness * size;
		aabb.upperBound += skinThickness * size;

		// Assign the new AABB.
		nodes[node].aabb = aabb;

		// Update the surface area and centroid.
		nodes[node].aabb.surfaceArea = nodes[node].aabb.computeSurfaceArea();

		// Insert a new leaf node.
		insertLeaf(node);

		return true;
	}

	std::vector<uint32_t> Tree::query(uint32_t particle)
	{
		// Make sure that this is a valid particle.
		assert(particleMap.contains(particle) && "The particle was not in the tree.");

		// Test overlap of particle AABB against all other particles.
		return query(particle, nodes[particleMap.find(particle)].aabb);
	}


	std::vector<uint32_t> Tree::query(uint32_t particleSelf, const AABB& aabb)
	{
		std::vector<uint32_t> particles;

		query(aabb, [&](uint32_t particle) {
			if (particle != particleSelf)
				particles.push_back(particle);

			return true; // continue the query
		});

		return particles;
	}

	std::vector<uint32_t> Tree::query(const AABB& aabb)
	{
		// Make sure the tree isn't empty.
		if (particleMap.size() == 0)
		{
			return std::vector<uint32_t>();
		}

		// Test overlap of AABB against all particles.
		return query(std::numeric_limits<uint32_t>::max(), aabb);
	}

	const AABB& Tree::getAABB(uint32_t particle)
	{
		return nodes[particleMap.find(particle)].aabb;
	}

	void Tree::insertLeaf(uint32_t leaf)
	{
		if (root == NULL_NODE)
		{
			root = leaf;
			nodes[root].parent = NULL_NODE;
			return;
		}

		// Find the best sibling for the node.

		AABB leafAABB = nodes[leaf].aabb;
		uint32_t index = root;

		while (!nodes[index].isLeaf())
		{
			// Extract the children of the node.
			uint32_t left = nodes[index].left;
			uint32_t right = nodes[index].right;

			float surfaceArea = nodes[index].aabb.getSurfaceArea();

			AABB combinedAABB;
			combinedAABB.merge(nodes[index].aabb, leafAABB);
			float combinedSurfaceArea = combinedAABB.getSurfaceArea();

			// Cost of creating a new parent for this node and the new leaf.
			float cost = 2.0 * combinedSurfaceArea;

			// Minimum cost of pushing the leaf further down the tree.
			float inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

			// Cost of descending to the left.
			float costLeft;
			if (nodes[left].isLeaf())
			{
				AABB aabb;
				aabb.merge(leafAABB, nodes[left].aabb);
				costLeft = aabb.getSurfaceArea() + inheritanceCost;
			}
			else
			{
				AABB aabb;
				aabb.merge(leafAABB, nodes[left].aabb);
				float oldArea = nodes[left].aabb.getSurfaceArea();
				float newArea = aabb.getSurfaceArea();
				costLeft = (newArea - oldArea) + inheritanceCost;
			}

			// Cost of descending to the right.
			float costRight;
			if (nodes[right].isLeaf())
			{
				AABB aabb;
				aabb.merge(leafAABB, nodes[right].aabb);
				costRight = aabb.getSurfaceArea() + inheritanceCost;
			}
			else
			{
				AABB aabb;
				aabb.merge(leafAABB, nodes[right].aabb);
				float oldArea = nodes[right].aabb.getSurfaceArea();
				float newArea = aabb.getSurfaceArea();
				costRight = (newArea - oldArea) + inheritanceCost;
			}

			// Descend according to the minimum cost.
			if ((cost < costLeft) && (cost < costRight)) break;

			// Descend.
			if (costLeft < costRight) index = left;
			else                      index = right;
		}

		uint32_t sibling = index;

		// Create a new parent.
		uint32_t oldParent = nodes[sibling].parent;
		uint32_t newParent = allocateNode();
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

			uint32_t left = nodes[index].left;
			uint32_t right = nodes[index].right;

			assert(left != NULL_NODE);
			assert(right != NULL_NODE);

			nodes[index].height = 1 + std::max(nodes[left].height, nodes[right].height);
			nodes[index].aabb.merge(nodes[left].aabb, nodes[right].aabb);

			index = nodes[index].parent;
		}
	}

	void Tree::removeLeaf(uint32_t leaf)
	{
		if (leaf == root)
		{
			root = NULL_NODE;
			return;
		}

		uint32_t parent = nodes[leaf].parent;
		uint32_t grandParent = nodes[parent].parent;
		uint32_t sibling;

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
			uint32_t index = grandParent;
			while (index != NULL_NODE)
			{
				index = balance(index);

				uint32_t left = nodes[index].left;
				uint32_t right = nodes[index].right;

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

	uint32_t Tree::balance(uint32_t node)
	{
		assert(node != NULL_NODE);

		if (nodes[node].isLeaf() || (nodes[node].height < 2))
			return node;

		uint32_t left = nodes[node].left;
		uint32_t right = nodes[node].right;

		assert(left < nodeCapacity);
		assert(right < nodeCapacity);

		int currentBalance = nodes[right].height - nodes[left].height;

		// Rotate right branch up.
		if (currentBalance > 1)
		{
			uint32_t rightLeft = nodes[right].left;
			uint32_t rightRight = nodes[right].right;

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
			uint32_t leftLeft = nodes[left].left;
			uint32_t leftRight = nodes[left].right;

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

	uint32_t Tree::computeHeight() const
	{
		return computeHeight(root);
	}

	uint32_t Tree::computeHeight(uint32_t node) const
	{
		assert(node < nodeCapacity);

		if (nodes[node].isLeaf()) return 0;

		uint32_t height1 = computeHeight(nodes[node].left);
		uint32_t height2 = computeHeight(nodes[node].right);

		return 1 + std::max(height1, height2);
	}

	uint32_t Tree::getHeight() const
	{
		if (root == NULL_NODE) return 0;
		return nodes[root].height;
	}

	uint32_t Tree::getNodeCount() const
	{
		return nodeCount;
	}

	uint32_t Tree::computeMaximumBalance() const
	{
		uint32_t maxBalance = 0;
		for (uint32_t i = 0; i < nodeCapacity; i++)
		{
			if (nodes[i].height <= 1)
				continue;

			assert(nodes[i].isLeaf() == false);

			uint32_t balance = std::abs(nodes[nodes[i].left].height - nodes[nodes[i].right].height);
			maxBalance = std::max(maxBalance, balance);
		}

		return maxBalance;
	}

	float Tree::computeSurfaceAreaRatio() const
	{
		if (root == NULL_NODE) return 0.0;

		float rootArea = nodes[root].aabb.computeSurfaceArea();
		float totalArea = 0.0;

		for (uint32_t i = 0; i < nodeCapacity; i++)
		{
			if (nodes[i].height < 0) continue;

			totalArea += nodes[i].aabb.computeSurfaceArea();
		}

		return totalArea / rootArea;
	}

	void Tree::validate() const
	{
#ifndef NDEBUG
		validateStructure(root);
		validateMetrics(root);

		uint32_t freeCount = 0;
		uint32_t freeIndex = freeList;

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

	void Tree::rebuild()
	{
		std::vector<uint32_t> nodeIndices(nodeCount);
		uint32_t count = 0;

		for (uint32_t i = 0; i < nodeCapacity; i++)
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
			float minCost = std::numeric_limits<float>::max();
			int iMin = -1, jMin = -1;

			for (uint32_t i = 0; i < count; i++)
			{
				AABB aabbi = nodes[nodeIndices[i]].aabb;

				for (uint32_t j = i + 1; j < count; j++)
				{
					AABB aabbj = nodes[nodeIndices[j]].aabb;
					AABB aabb;
					aabb.merge(aabbi, aabbj);
					float cost = aabb.getSurfaceArea();

					if (cost < minCost)
					{
						iMin = i;
						jMin = j;
						minCost = cost;
					}
				}
			}

			uint32_t index1 = nodeIndices[iMin];
			uint32_t index2 = nodeIndices[jMin];

			uint32_t parent = allocateNode();
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

	void Tree::validateStructure(uint32_t node) const
	{
		if (node == NULL_NODE) return;

		if (node == root) assert(nodes[node].parent == NULL_NODE);

		uint32_t left = nodes[node].left;
		uint32_t right = nodes[node].right;

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

	void Tree::validateMetrics(uint32_t node) const
	{
		if (node == NULL_NODE) return;

		uint32_t left = nodes[node].left;
		uint32_t right = nodes[node].right;

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

		AABB aabb;
		aabb.merge(nodes[left].aabb, nodes[right].aabb);

		for (uint32_t i = 0; i < 3; i++)
		{
			assert(aabb.lowerBound[i] == nodes[node].aabb.lowerBound[i]);
			assert(aabb.upperBound[i] == nodes[node].aabb.upperBound[i]);
		}

		validateMetrics(left);
		validateMetrics(right);
	}
}
