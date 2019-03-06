/****************************************************************************
* VCGLib                                                            o o     *
* Visual and Computer Graphics Library                            o     o   *
*                                                                _   O  _   *
* Copyright(C) 2004-2016                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

// Ported to the Unreal Engine by Timothy Davison, 2019.
// Based on: https://github.com/cnr-isti-vclab/vcglib/blob/master/vcg/space/index/kdtree/kdtree_face.h

#pragma  once

#include <vector>
#include <array>
#include <unordered_set>

#include <Vector.h>
#include "Algorithm/SurfaceIndex.h"
#include <Box.h>

namespace vcg {

	struct TriangleSample
	{
		std::array<FVector, 3> triangle;

		FSurfaceIndex surfaceIndex;

		// \return true if closer than the passed distance, false otherwise.
		bool closestPoint(const FVector& query, float& distance, FVector& nearestPoint) 
		{
			const float in_distance = distance;

			nearestPoint = FMath::ClosestPointOnTriangleToPoint(query, triangle[0], triangle[1], triangle[2]);

			distance = FVector::Dist(nearestPoint, query);

			if (distance > in_distance)
				return false;
			else
				return true; 
		}

		bool hackIntersectsBox(const FBox& box)
		{
			return box.IsInsideOrOn(triangle[0]) || box.IsInsideOrOn(triangle[1]) || box.IsInsideOrOn(triangle[2]);
		}
	};

	/**
	* This class allows to create a Kd-Tree thought to perform the neighbour query using the mesh faces (closest search).
	* The class implemetantion is thread-safe.
	*/
	class KdTreeFace
	{
	public:

		typedef typename float Scalar;

		typedef typename TriangleSample TriangleSample;

		class Node
		{
		public:
			Scalar splitValue;
			unsigned int firstChildId : 24;
			unsigned int dim : 2;
			unsigned int leaf : 1;
			FBox aabb = FBox(EForceInit::ForceInitToZero);
			std::vector<TriangleSample*> list;
		};
		typedef std::vector<Node> NodeListPointer;

	public:
		std::vector<TriangleSample> triangles;

		KdTreeFace(std::vector<TriangleSample>& trianglesIin, FBox aabb, unsigned int maxObjPerCell = 64, unsigned int maxDepth = 64) : epsilon(std::numeric_limits<Scalar>::epsilon())
		{
			triangles = trianglesIin;

			targetCellSize = maxObjPerCell;
			targetMaxDepth = maxDepth;
			mNodes.resize(1);
			Node& node = mNodes.back();
			node.leaf = 0;
			node.aabb = aabb;
			node.aabb = node.aabb.ExpandBy(FVector(epsilon));

			std::transform(triangles.begin(), triangles.end(), std::back_inserter(node.list),
				[](TriangleSample &s) { return &s; }
			);

			numLevel = createTree(0, 1);
		};

		~KdTreeFace()
		{

		};

		struct ObjectMarker
		{
			std::unordered_set<TriangleSample*> marks;

			bool IsMarked(TriangleSample * sample) 
			{ 
				return marks.find(sample) != marks.end();
			}

			void Mark(TriangleSample * sample)
			{
				marks.insert(sample);
			}
		};

		TriangleSample* doQueryClosest(const FVector& queryPoint, FVector& narestPoint, Scalar& dist, ObjectMarker& marker, Scalar maxDist = std::numeric_limits<Scalar>::max())
		{
			float maxDistSqrd = maxDist * maxDist;

			{
				float distSqrd = mNodes[0].aabb.ComputeSquaredDistanceToPoint(queryPoint);

				if (maxDist < std::numeric_limits<Scalar>::max() &&
					!mNodes[0].aabb.IsInside(queryPoint) &&
					distSqrd >= maxDistSqrd
					)
				{
					dist = std::sqrt(distSqrd);
					return nullptr;
				}
			}


			std::vector<QueryNode> mNodeStack(numLevel + 1);
			mNodeStack[0].nodeId = 0;
			mNodeStack[0].sq = 0.f;
			unsigned int count = 1;

			float minDist = maxDist;
			FVector p;
			TriangleSample* face = nullptr;
			while (count)
			{
				QueryNode& qnode = mNodeStack[count - 1];
				Node& node = mNodes[qnode.nodeId];

				if (qnode.sq < minDist)
				{
					if (node.leaf)
					{
						--count; // pop
						for (int i = 0; i < node.list.size(); i++)
						{
							if (!marker.IsMarked(node.list[i]))
							{
								marker.Mark(node.list[i]);
								Scalar tempDist = minDist;
								FVector tempP;

								if (node.list[i]->closestPoint(queryPoint, tempDist, tempP))
								{
									if (tempDist < minDist)
									{
										minDist = tempDist;
										p = tempP;
										face = node.list[i];
									}
								}
							}
						}
					}
					else
					{
						// replace the stack top by the farthest and push the closest
						float new_off = queryPoint[node.dim] - node.splitValue;
						float abs_off = abs(new_off);
						if (abs_off < minDist)
						{
							if (new_off < 0.)
							{
								mNodeStack[count].nodeId = node.firstChildId;
								qnode.nodeId = node.firstChildId + 1;

								const float distSqrd = mNodes[node.firstChildId + 1].aabb.ComputeSquaredDistanceToPoint(queryPoint);

								new_off = std::sqrt(distSqrd);
							}
							else
							{
								mNodeStack[count].nodeId = node.firstChildId + 1;
								qnode.nodeId = node.firstChildId;

								const float distSqrd = mNodes[node.firstChildId].aabb.ComputeSquaredDistanceToPoint(queryPoint);

								new_off = std::sqrt(distSqrd);
							}
							mNodeStack[count].sq = qnode.sq;
							qnode.sq = new_off;
							++count;
						}
						else
						{
							if (new_off < 0.)
								qnode.nodeId = node.firstChildId;
							else
								qnode.nodeId = node.firstChildId + 1;
						}
					}
				}
				else
				{
					// pop
					--count;
				}
			}
			dist = minDist;
			narestPoint = p;
			return face;
		}

	protected:

		// element of the stack
		struct QueryNode
		{
			QueryNode() {}
			QueryNode(unsigned int id) : nodeId(id) {}
			unsigned int nodeId;  // id of the next node
			Scalar sq;            // distance to the next node
		};


		int createTree(unsigned int nodeId, unsigned int level)
		{
			Node& node = mNodes[nodeId];
			FVector diag = node.aabb.Max - node.aabb.Min;
			unsigned int dim;
			if (diag.X > diag.Y)
				dim = diag.X > diag.Z ? 0 : 2;
			else
				dim = diag.Y > diag.Z ? 1 : 2;

			node.splitValue = Scalar(0.5*(node.aabb.Max[dim] + node.aabb.Min[dim]));
			node.dim = dim;

			FBox leftBox(EForceInit::ForceInitToZero);
			FBox rightBox(EForceInit::ForceInitToZero);

			leftBox += node.aabb.Min;
			rightBox += node.aabb.Max;

			if (node.dim == 0)
			{
				leftBox += FVector(node.splitValue, node.aabb.Max[1], node.aabb.Max[2]);
				rightBox += FVector(node.splitValue, node.aabb.Min[1], node.aabb.Min[2]);
			}
			else if (node.dim == 1)
			{
				leftBox += FVector(node.aabb.Max[0], node.splitValue, node.aabb.Max[2]);
				rightBox += FVector(node.aabb.Min[0], node.splitValue, node.aabb.Min[2]);
			}
			else if (node.dim == 2)
			{
				leftBox += FVector(node.aabb.Max[0], node.aabb.Max[1], node.splitValue);
				rightBox += FVector(node.aabb.Min[0], node.aabb.Min[1], node.splitValue);
			}

			leftBox = leftBox.ExpandBy(FVector(epsilon, epsilon, epsilon));
			rightBox = rightBox.ExpandBy(FVector(epsilon, epsilon, epsilon));

			node.firstChildId = mNodes.size();
			int firstChildId = node.firstChildId;
			mNodes.resize(mNodes.size() + 2);
			Node& parent = mNodes[nodeId];
			Node& leftChild = mNodes[firstChildId];
			Node& rightChild = mNodes[firstChildId + 1];

			leftChild.aabb.Init(); 
			rightChild.aabb.Init(); 

			for (int i = 0; i < parent.list.size(); i++)
			{
				unsigned int state = 0;
				TriangleSample * fp = parent.list[i];
				for (int j = 0; j < 3; j++)
				{
					if (fp->triangle[j][dim] < parent.splitValue)
						state |= (1 << 0);
					else if (fp->triangle[j][dim] > parent.splitValue)
						state |= (1 << 1);
					else
					{
						state |= (1 << 0);
						state |= (1 << 1);
					}
				}
				if (state & (1 << 0))
				{
					leftChild.list.push_back(fp);
					leftChild.aabb += fp->triangle[0];
					leftChild.aabb += fp->triangle[1];
					leftChild.aabb += fp->triangle[2];
				}
				if (state & (1 << 1))
				{
					rightChild.list.push_back(fp);
					rightChild.aabb += fp->triangle[0];
					rightChild.aabb += fp->triangle[1];
					rightChild.aabb += fp->triangle[2];
				}
			}
			parent.list.clear();
			leftChild.aabb = leftChild.aabb.Overlap(leftBox);
			rightChild.aabb = rightChild.aabb.Overlap(rightBox);

			int leftLevel, rightLevel;
			{
				if (leftChild.list.size() <= targetCellSize || level >= targetMaxDepth)
				{
					leftChild.leaf = 1;
					leftLevel = level;
				}
				else
				{
					leftChild.leaf = 0;
					leftLevel = createTree(firstChildId, level + 1);
				}
			}

			{
				Node& rightChild = mNodes[firstChildId + 1];
				if (rightChild.list.size() <= targetCellSize || level >= targetMaxDepth)
				{
					rightChild.leaf = 1;
					rightLevel = level;
				}
				else
				{
					rightChild.leaf = 0;
					rightLevel = createTree(firstChildId + 1, level + 1);
				}
			}
			if (leftLevel > rightLevel)
				return leftLevel;
			return rightLevel;
		};


	protected:

		NodeListPointer mNodes; //kd-tree nodes
		unsigned int numLevel;
		const Scalar epsilon;
		unsigned int targetCellSize;
		unsigned int targetMaxDepth;
	};

}

