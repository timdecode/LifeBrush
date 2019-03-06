// Copyright 2016 Code Monkey Castle, all rights reserved.

#include "LifeBrush.h"

#include <algorithm>
#include <set>
#include <map>
#include <vector>
#include <array>

#include "Utility.h"

#include "Algorithm/dbscan.hpp"

#include "MeshFactory.h"

void SmoothMeshFactory::sortTriangleIndices()
{
	int32 n = indices.Num();

	if (n == 0)
		return;

	check(n % 3 == 0);

	for (int i = 0; i < n; i += 3)
	{
		int32* triangle = &indices[i];

		// find the smallest index in this triangle
		int32 min_c = 0;

		for (int32 c = 1; c < 3; ++c)
		{
			if (triangle[c] < triangle[min_c])
				min_c = c;
		}

		if (min_c == 0)
			continue;

		// spin the triangle to the smallest
		std::rotate(triangle, triangle + min_c, triangle + 3);
	}
}

void SmoothMeshFactory::removeDuplicateFaces()
{
	int32 n = indices.Num();

	if (n == 0)
		return;

	check(n % 3 == 0);

	sortTriangleIndices();

	// if we were clever, we could do this with just indices. But I am lazy today.
	std::set<FIntVector > uniqueFaces;

	TArray<int32> newIndices;

	for (int i = 0; i < n; i += 3)
	{
		int32* triangle = &indices[i];

		FIntVector face(triangle[0], triangle[1], triangle[2]);

		bool isUnique = uniqueFaces.find(face) == uniqueFaces.end();

		if (isUnique)
		{
			// add a unique face
			newIndices.Add(triangle[0]);
			newIndices.Add(triangle[1]);
			newIndices.Add(triangle[2]);

			uniqueFaces.insert(face);
		}
	}

	indices = newIndices;
}

void SmoothMeshFactory::append(SmoothMeshFactory& other)
{
	const int32 offset = vertices.Num();

	vertices.Append(other.vertices);
	normals.Append(other.normals);
	uvs.Append(other.uvs);
	tangents.Append(other.tangents);

	indices.Reserve(indices.Num() + other.indices.Num());
	for (auto other_i : other.indices)
		indices.Add(other_i + offset);
}

SmoothMeshFactory SmoothMeshFactory::createMesh_withoutDuplicateVertices()
{
	SmoothMeshFactory result;

	// track duplicate vertices
	std::map<FVector, uint32> vertexLookup;
	std::vector<uint32> deduplicatedVertexIndices;

	// copy positions
	const auto nv = vertices.Num();

	for (int vi = 0; vi < nv; ++vi)
	{
		const FVector& v = vertices[vi];

		auto found = vertexLookup.find(v);

		// it's new
		if (found == vertexLookup.end())
		{
			int32 new_vi = result.vertices.Num();

			vertexLookup[v] = new_vi;

			deduplicatedVertexIndices.push_back(new_vi);

			result.vertices.Add(v);

			if( normals.Num() > vi ) result.normals.Add(normals[vi]);
			if( tangents.Num() > vi ) result.tangents.Add(tangents[vi]);
			if (uvs.Num() > vi) result.uvs.Add(uvs[vi]);
		}
		else
			deduplicatedVertexIndices.push_back(found->second);
	}

	// build the new index buffer
	const auto ni = indices.Num();

	for (int32 i = 0; i < ni; ++i)
	{
		result.indices.Add(deduplicatedVertexIndices[i]);
	}

	return result;
}

SmoothMeshFactory SmoothMeshFactory::createMesh_mergeNearbyVertices(float minVertexDistance /*= 0.1f*/)
{
	SmoothMeshFactory result;

	if (vertices.Num() < 3)
		return result;

	// cluster nearby vertices

	using DBScanElements = rg::DBScan<float, FVector>;

	std::vector< std::vector<DBScanElements::PointAndIndex> > clusters_out;
	std::vector<DBScanElements::PointAndIndex> noise_out;

	std::vector<FVector> vertices_;
	vertices_.reserve(vertices.Num());
	for (auto& v : vertices)
		vertices_.push_back(v);

	FVector_kdTree kdTree(vertices_);
	DBScanElements::dbscan(kdTree, minVertexDistance, 1, clusters_out, noise_out);

	// find the centroids
	std::vector<FVector> centroids;
	{
		centroids.reserve(clusters_out.size());

		for (auto& cluster : clusters_out)
		{
			if (cluster.size() == 0)
				continue;

			FVector centroid = FVector::ZeroVector;

			for (auto& c : cluster)
			{
				centroid += c.point;
			}

			centroid /= float(cluster.size());

			centroids.push_back(centroid);
		}
	}

	// nuke the guys too far from the centroid
	// track duplicate vertices
	std::vector<uint32> deduplicatedVertexIndices(vertices.Num());

	for (int i = 0; i < clusters_out.size(); ++i)
	{
		auto& cluster = clusters_out[i];
		FVector& centroid = centroids[i];

		size_t nearestIndex;
		float distance;

		kdTree.kdTree.knnSearch(&centroid[0], 1, &nearestIndex, &distance);

		for (auto& c : cluster)
		{
			deduplicatedVertexIndices[c.index] = i;
		}

		result.vertices.Add(centroid);

		if (normals.Num() > nearestIndex) result.normals.Add(normals[nearestIndex]);
		if (tangents.Num() > nearestIndex) result.tangents.Add(tangents[nearestIndex]);
		if (uvs.Num() > nearestIndex) result.uvs.Add(uvs[nearestIndex]);
	}

	// build the result
	// build the new index buffer
	const auto ni = indices.Num();

	for (int32 i = 0; i < ni; ++i)
	{
		result.indices.Add(deduplicatedVertexIndices[i]);
	}

	return result;
}

void SmoothMeshFactory::transform(const FTransform transform)
{
	FQuat normalRotation = transform.GetRotation();

	for (auto& v : vertices)
		v = transform.TransformPosition(v);

	for (auto& n : normals)
		n = normalRotation.RotateVector(n);

	for (auto& t : tangents)
		t.TangentX = normalRotation.RotateVector(t.TangentX);
}
