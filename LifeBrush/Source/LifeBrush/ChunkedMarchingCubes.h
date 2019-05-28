// Copyright 2017 Code Monkey Castle, all rights reserved.

#pragma once

#include <vector>
#include <map>

#include "MeshFactory.h"
#include "Utility.h"

#include "UniformGrid.h"

// This code is adapted from Paul Bourke's Marching Cubes Guide 
// http://paulbourke.net/geometry/polygonise/
// Accessed August 27, 2016





class ChunkedMarchingCubes
{
public:


	struct GridCell
	{
	public:
		FVector p[8];
		FVector4 v[8];
	};

	struct Vertex
	{
	public:
		FVector p;	// point
		FVector n;	// normal
	};

	struct Triangle
	{
	public:
		Vertex p0;
		Vertex p1;
		Vertex p2;
	};

public:

	struct ChunkGeometry
	{
		SmoothMeshFactory factory;

		// vector is index by triangle indices (i_0 / 3 of a triangle)
		std::vector<FIntVector> cellIndices;
	};

	// the return map is indexed by gridIndices
	template<typename ElementType>
	static std::unordered_map<FIntVector, ChunkGeometry> marchingCubes_wholeChunkGridAccessSlow(
		float isoLevel,
		ChunkGrid<ElementType>& grid,
		FIntVector gridMin,
		FIntVector gridMax, // inclusive (that is, up to and including)
		const FVector2D& uvScale,
		bool reverseWinding = false,
		bool flatShading = false)
	{
		const FVector cellSize = grid.cellSize();

		FIntVector const vertexIndices[8] = {
			FIntVector(0,0,0),
			FIntVector(1,0,0),
			FIntVector(1,1,0),
			FIntVector(0,1,0),

			FIntVector(0,0,1),
			FIntVector(1,0,1),
			FIntVector(1,1,1),
			FIntVector(0,1,1)
		};

		FVector const vertexOffsets[] = {
			FVector(vertexIndices[0]) * cellSize,
			FVector(vertexIndices[1]) * cellSize,
			FVector(vertexIndices[2]) * cellSize,
			FVector(vertexIndices[3]) * cellSize,

			FVector(vertexIndices[4]) * cellSize,
			FVector(vertexIndices[5]) * cellSize,
			FVector(vertexIndices[6]) * cellSize,
			FVector(vertexIndices[7]) * cellSize
		};

		typedef std::remove_reference<decltype(grid)>::type GridType;


		std::unordered_map<FIntVector, ChunkGeometry> result;

		float normalSign = reverseWinding ? -1.0f : 1.0f;

		GridCell cell;
		Triangle triangles[5];


		FIntVector si;
		for (si.Z = gridMin.Z; si.Z <= gridMax.Z; si.Z += 1)
		{
			for (si.Y = gridMin.Y; si.Y <= gridMax.Y; si.Y += 1)
			{
				for (si.X = gridMin.X; si.X <= gridMax.X; si.X += 1)
				{
					PaddedUniformGrid<ElementType>& chunk = grid.chunkAtGridIndex(si);

					FIntVector chunkIndex = grid.indexOfChunk(si);
					FIntVector cellIndex = grid.indexOfCell(si);

					ChunkGeometry& geometry = result[chunkIndex];
					SmoothMeshFactory& builder = geometry.factory;
					std::vector<int32>& cellIndices = geometry.cellIndices;

					FVector p = grid.samplePoint(si);

					// cell position
					for (int c = 0; c < 8; ++c)
						cell.p[c] = p + vertexOffsets[c];

					// indices
					FIntVector indices[8];
					for (int i = 0; i < 8; ++i)
						indices[i] = si + vertexIndices[i];

					// normal
					FVector g[8];
					for (int c = 0; c < 8; ++c)
						g[c] = normalSign * grid.normalAt(indices[c]);

					// cell normal and value
					for (int c = 0; c < 8; ++c)
						cell.v[c] = FVector4(g[c].X, g[c].Y, g[c].Z, grid(indices[c]));

					int n = polygonise(cell, isoLevel, triangles);

					for (int i = 0; i < n; i++)
					{
						Triangle& t = triangles[i];

						int ux = 0;
						int uy = 2;

						// face normal
						FVector normal = FVector::CrossProduct(t.p2.p - t.p0.p, t.p1.p - t.p0.p);
						normal.Normalize();

						std::vector<float> alignment{
							FMath::Abs(normal.X),
							FMath::Abs(normal.Y),
							FMath::Abs(normal.Z) };

						auto max = std::max_element(alignment.begin(), alignment.end());
						int maxI = std::distance(alignment.begin(), max);

						if (maxI == 0)
						{
							ux = 1;
							uy = 2;
						}
						else if (maxI == 1)
						{
							ux = 0;
							uy = 2;
						}
						else if (maxI == 2)
						{
							ux = 0;
							uy = 1;
						}

						FVector2D u0(t.p0.p[ux], t.p0.p[uy]);
						FVector2D u1(t.p1.p[ux], t.p1.p[uy]);
						FVector2D u2(t.p2.p[ux], t.p2.p[uy]);

						u0 *= uvScale;
						u1 *= uvScale;
						u2 *= uvScale;

						FVector n0, n1, n2;

						if (flatShading)
						{
							n0 = -normal;
							n1 = -normal;
							n2 = -normal;
						}
						else
						{
							n0 = t.p0.n;
							n1 = t.p1.n;
							n2 = t.p2.n;
						}

						// don't forget that unity uses a clockwise winding order
						if (reverseWinding)
						{

							builder.pushTriangle(
								t.p0.p, t.p2.p, t.p1.p,
								n0, n2, n1,
								u0, u2, u1);

						}
						else
						{

							builder.pushTriangle(
								t.p0.p, t.p1.p, t.p2.p,
								n0, n1, n2,
								u0, u1, u2);
						}

						cellIndices.push_back(si);
					}
				}
			}
		}

		return result;
	}

	// the return map is indexed by gridIndices
	template<typename ElementType>
	static ChunkGeometry marchingCubes(
		float isoLevel,
		PaddedUniformGrid<ElementType>& grid,
		FIntVector gridMin,
		FIntVector gridMax, // inclusive (that is, up to and including)
		const FVector2D& uvScale,
		bool reverseWinding = false,
		bool flatShading = false)
	{
		const FVector cellSize = grid.cellSize();

		FIntVector const vertexIndices[8] = {
			FIntVector(0,0,0),
			FIntVector(1,0,0),
			FIntVector(1,1,0),
			FIntVector(0,1,0),

			FIntVector(0,0,1),
			FIntVector(1,0,1),
			FIntVector(1,1,1),
			FIntVector(0,1,1)
		};

		FVector const vertexOffsets[] = {
			FVector(vertexIndices[0]) * cellSize,
			FVector(vertexIndices[1]) * cellSize,
			FVector(vertexIndices[2]) * cellSize,
			FVector(vertexIndices[3]) * cellSize,

			FVector(vertexIndices[4]) * cellSize,
			FVector(vertexIndices[5]) * cellSize,
			FVector(vertexIndices[6]) * cellSize,
			FVector(vertexIndices[7]) * cellSize
		};

		typedef std::remove_reference<decltype(grid)>::type GridType;


		ChunkGeometry geometry;

		float normalSign = reverseWinding ? -1.0f : 1.0f;

		GridCell cell;
		Triangle triangles[5];


		FIntVector si;
		for (si.Z = gridMin.Z; si.Z <= gridMax.Z; si.Z += 1)
		{
			for (si.Y = gridMin.Y; si.Y <= gridMax.Y; si.Y += 1)
			{
				for (si.X = gridMin.X; si.X <= gridMax.X; si.X += 1)
				{
					//if (si.X == gridMax.X && si.Y == gridMax.Y && si.Z == gridMax.Z)
					//{
					//	float bob = 2;

					//	bob += 2;
					//}

					SmoothMeshFactory& builder = geometry.factory;
					std::vector<FIntVector>& cellIndices = geometry.cellIndices;

					FVector p = grid.samplePoint(si);

					// cell position
					for (int c = 0; c < 8; ++c)
						cell.p[c] = p + vertexOffsets[c];

					// indices
					FIntVector indices[8];
					for (int i = 0; i < 8; ++i)
						indices[i] = si + vertexIndices[i];

					// normal
					FVector g[8];
					for (int c = 0; c < 8; ++c)
						g[c] = normalSign * grid.normalAt(indices[c]);

					// cell normal and value
					for (int c = 0; c < 8; ++c)
						cell.v[c] = FVector4(g[c].X, g[c].Y, g[c].Z, grid(indices[c]));

					int n = polygonise(cell, isoLevel, triangles);

					for (int i = 0; i < n; i++)
					{
						Triangle& t = triangles[i];

						int ux = 0;
						int uy = 2;

						// face normal
						FVector normal = FVector::CrossProduct(t.p2.p - t.p0.p, t.p1.p - t.p0.p);
						normal.Normalize();

						std::vector<float> alignment{
							FMath::Abs(normal.X),
							FMath::Abs(normal.Y),
							FMath::Abs(normal.Z) };

						auto max = std::max_element(alignment.begin(), alignment.end());
						int maxI = std::distance(alignment.begin(), max);

						if (maxI == 0)
						{
							ux = 1;
							uy = 2;
						}
						else if (maxI == 1)
						{
							ux = 0;
							uy = 2;
						}
						else if (maxI == 2)
						{
							ux = 0;
							uy = 1;
						}

						FVector2D u0(t.p0.p[ux], t.p0.p[uy]);
						FVector2D u1(t.p1.p[ux], t.p1.p[uy]);
						FVector2D u2(t.p2.p[ux], t.p2.p[uy]);

						u0 *= uvScale;
						u1 *= uvScale;
						u2 *= uvScale;

						FVector n0, n1, n2;

						if (flatShading)
						{
							n0 = -normal;
							n1 = -normal;
							n2 = -normal;
						}
						else
						{
							n0 = t.p0.n;
							n1 = t.p1.n;
							n2 = t.p2.n;
						}

						// don't forget that unity uses a clockwise winding order
						if (reverseWinding)
						{

							builder.pushTriangle(
								t.p0.p, t.p2.p, t.p1.p,
								n0, n2, n1,
								u0, u2, u1);

						}
						else
						{

							builder.pushTriangle(
								t.p0.p, t.p1.p, t.p2.p,
								n0, n1, n2,
								u0, u1, u2);
						}

						cellIndices.push_back(si);
					}
				}
			}
		}

		return geometry;
	}

	// the return map is indexed by gridIndices
	template<typename ElementType>
	static std::unordered_map<FIntVector, ChunkGeometry> marchingCubes_padded(
		float isoLevel,
		ChunkGrid<ElementType>& grid,
		FIntVector gridMin,
		FIntVector gridMax,
		const FVector2D& uvScale,
		bool reverseWinding = false,
		bool flatShading = false)
	{
		const FVector cellSize = grid.cellSize();

		FIntVector const vertexIndices[8] = {
			FIntVector(0,0,0),
			FIntVector(1,0,0),
			FIntVector(1,1,0),
			FIntVector(0,1,0),

			FIntVector(0,0,1),
			FIntVector(1,0,1),
			FIntVector(1,1,1),
			FIntVector(0,1,1)
		};

		FVector const vertexOffsets[] = {
			FVector(vertexIndices[0]) * cellSize,
			FVector(vertexIndices[1]) * cellSize,
			FVector(vertexIndices[2]) * cellSize,
			FVector(vertexIndices[3]) * cellSize,

			FVector(vertexIndices[4]) * cellSize,
			FVector(vertexIndices[5]) * cellSize,
			FVector(vertexIndices[6]) * cellSize,
			FVector(vertexIndices[7]) * cellSize
		};

		typedef std::remove_reference<decltype(grid)>::type GridType;

		std::unordered_map<FIntVector, ChunkGeometry> result;

		float normalSign = reverseWinding ? -1.0f : 1.0f;

		FIntVector fromChunkIndex = grid.indexOfChunk(gridMin);
		FIntVector toChunkIndex = grid.indexOfChunk(gridMax);

		FIntVector fromCellIndex = grid.indexOfCell(gridMin);
		FIntVector toCellIndex = grid.indexOfCell(gridMax);

		FIntVector chunkIndex = fromChunkIndex;

		for (chunkIndex.Z = fromChunkIndex.Z; chunkIndex.Z <= toChunkIndex.Z; chunkIndex.Z++)
		{
			for (chunkIndex.Y = fromChunkIndex.Y; chunkIndex.Y <= toChunkIndex.Y; chunkIndex.Y++)
			{
				for (chunkIndex.X = fromChunkIndex.X; chunkIndex.X <= toChunkIndex.X; chunkIndex.X++)
				{
					PaddedUniformGrid<float>& chunk = grid.chunkAtChunkIndex(chunkIndex);

					// calculate the local cell index start and end
					FIntVector localFromCell;

					for (int c = 0; c < 3; ++c)
					{
						if (chunkIndex[c] == fromChunkIndex[c])
							localFromCell[c] = fromCellIndex[c];
						else
							localFromCell[c] = 0;
					}

					FIntVector localToCell;

					for (int c = 0; c < 3; ++c)
					{
						int maxDim = grid.chunkDimensions()[c];

						if (chunkIndex[c] == toChunkIndex[c])
							localToCell[c] = toCellIndex[c] >= maxDim ? maxDim - 1 : toCellIndex[c];
						else
							localToCell[c] = maxDim - 1;
					}

					result[chunkIndex] = marchingCubes(isoLevel, chunk, localFromCell, localToCell, uvScale, reverseWinding, flatShading);
				}
			}
		}

		return result;
	}

	// the return map is indexed by gridIndices
	template<typename ElementType>
	static ChunkGeometry marchingCubes_byChunkAccessFast(
		float isoLevel,
		ChunkGrid<ElementType>& grid,
		PaddedUniformGrid<ElementType>& chunk,
		FIntVector chunkIndex,
		FIntVector cellMin,
		FIntVector cellMax, // inclusive (that is, up to and including)
		const FVector2D& uvScale,
		bool reverseWinding = false,
		bool flatShading = false)
	{
		const FVector cellSize = grid.cellSize();

		FIntVector const vertexIndices[8] = {
			FIntVector(0,0,0),
			FIntVector(1,0,0),
			FIntVector(1,1,0),
			FIntVector(0,1,0),

			FIntVector(0,0,1),
			FIntVector(1,0,1),
			FIntVector(1,1,1),
			FIntVector(0,1,1)
		};

		FVector const vertexOffsets[] = {
			FVector(vertexIndices[0]) * cellSize,
			FVector(vertexIndices[1]) * cellSize,
			FVector(vertexIndices[2]) * cellSize,
			FVector(vertexIndices[3]) * cellSize,

			FVector(vertexIndices[4]) * cellSize,
			FVector(vertexIndices[5]) * cellSize,
			FVector(vertexIndices[6]) * cellSize,
			FVector(vertexIndices[7]) * cellSize
		};

		typedef std::remove_reference<decltype(grid)>::type GridType;


		ChunkGeometry geometry;

		float normalSign = reverseWinding ? -1.0f : 1.0f;

		GridCell cell;
		Triangle triangles[5];

		const FIntVector dim = chunk.dimensions();
		const FIntVector chunkGridStart = grid.gridIndexFromChunkIndex(chunkIndex);

		FIntVector si;
		for (si.Z = cellMin.Z; si.Z <= cellMax.Z; si.Z += 1)
		{
			for (si.Y = cellMin.Y; si.Y <= cellMax.Y; si.Y += 1)
			{
				for (si.X = cellMin.X; si.X <= cellMax.X; si.X += 1)
				{
					SmoothMeshFactory& builder = geometry.factory;
					std::vector<FIntVector>& cellIndices = geometry.cellIndices;

					FVector p = chunk.samplePoint(si);

					// cell position
					for (int c = 0; c < 8; ++c)
						cell.p[c] = p + vertexOffsets[c];

					// indices
					FIntVector indices[8];
					for (int i = 0; i < 8; ++i)
						indices[i] = si + vertexIndices[i];

					// normal
					FVector g[8];
					for (int c = 0; c < 8; ++c)
					{
						const FIntVector& index = indices[c];

						if (index.X < dim.X - 1 && index.Y < dim.Y - 1 && index.Z < dim.Z - 1 &&
							index.X > 0 && index.Y > 0 && index.Z > 0 )
						//if (index.X <= dim.X && index.Y <= dim.Y && index.Z <= dim.Z &&
						//	index.X >= 0 && index.Y >= 0 && index.Z >= 0) 
						{
							g[c] = normalSign * chunk.normalAt(index);
							cell.v[c] = FVector4(g[c].X, g[c].Y, g[c].Z, chunk(indices[c]));
						}
						else
						{
							FIntVector gridIndex = chunkGridStart + index;

							g[c] = normalSign * grid.normalAt(gridIndex);
							cell.v[c] = FVector4(g[c].X, g[c].Y, g[c].Z, grid(gridIndex));
						}
					}

					int n = polygonise(cell, isoLevel, triangles);

					for (int i = 0; i < n; i++)
					{
						Triangle& t = triangles[i];

						int ux = 0;
						int uy = 2;

						// face normal
						FVector normal = FVector::CrossProduct(t.p2.p - t.p0.p, t.p1.p - t.p0.p);
						normal.Normalize();

						std::vector<float> alignment{
							FMath::Abs(normal.X),
							FMath::Abs(normal.Y),
							FMath::Abs(normal.Z) };

						auto max = std::max_element(alignment.begin(), alignment.end());
						int maxI = std::distance(alignment.begin(), max);

						if (maxI == 0)
						{
							ux = 1;
							uy = 2;
						}
						else if (maxI == 1)
						{
							ux = 0;
							uy = 2;
						}
						else if (maxI == 2)
						{
							ux = 0;
							uy = 1;
						}

						FVector2D u0(t.p0.p[ux], t.p0.p[uy]);
						FVector2D u1(t.p1.p[ux], t.p1.p[uy]);
						FVector2D u2(t.p2.p[ux], t.p2.p[uy]);

						u0 *= uvScale;
						u1 *= uvScale;
						u2 *= uvScale;

						FVector n0, n1, n2;

						if (flatShading)
						{
							n0 = -normal;
							n1 = -normal;
							n2 = -normal;
						}
						else
						{
							n0 = t.p0.n;
							n1 = t.p1.n;
							n2 = t.p2.n;
						}

						// don't forget that unity uses a clockwise winding order
						if (reverseWinding)
						{

							builder.pushTriangle(
								t.p0.p, t.p2.p, t.p1.p,
								n0, n2, n1,
								u0, u2, u1);

						}
						else
						{

							builder.pushTriangle(
								t.p0.p, t.p1.p, t.p2.p,
								n0, n1, n2,
								u0, u1, u2);
						}

						cellIndices.push_back(si);
					}
				}
			}
		}

		return geometry;
	}

	// the return map is indexed by gridIndices
	template<typename ElementType>
	static std::map<FIntVector, ChunkGeometry> marchingCubes_byChunk(
		float isoLevel,
		ChunkGrid<ElementType>& grid,
		FIntVector gridMin,
		FIntVector gridMax,
		const FVector2D& uvScale,
		bool reverseWinding = false,
		bool flatShading = false)
	{
		const FVector cellSize = grid.cellSize();

		FIntVector const vertexIndices[8] = {
			FIntVector(0,0,0),
			FIntVector(1,0,0),
			FIntVector(1,1,0),
			FIntVector(0,1,0),

			FIntVector(0,0,1),
			FIntVector(1,0,1),
			FIntVector(1,1,1),
			FIntVector(0,1,1)
		};

		FVector const vertexOffsets[] = {
			FVector(vertexIndices[0]) * cellSize,
			FVector(vertexIndices[1]) * cellSize,
			FVector(vertexIndices[2]) * cellSize,
			FVector(vertexIndices[3]) * cellSize,

			FVector(vertexIndices[4]) * cellSize,
			FVector(vertexIndices[5]) * cellSize,
			FVector(vertexIndices[6]) * cellSize,
			FVector(vertexIndices[7]) * cellSize
		};

		typedef std::remove_reference<decltype(grid)>::type GridType;

		std::map<FIntVector, ChunkGeometry> result;

		float normalSign = reverseWinding ? -1.0f : 1.0f;

		FIntVector fromChunkIndex = grid.indexOfChunk(gridMin);
		FIntVector toChunkIndex = grid.indexOfChunk(gridMax);

		FIntVector fromCellIndex = grid.indexOfCell(gridMin);
		FIntVector toCellIndex = grid.indexOfCell(gridMax);

		FIntVector chunkIndex = fromChunkIndex;

		struct IndexAndChunk
		{
			FIntVector index;
			ChunkGeometry geometry;
		};

		std::vector<TFuture<IndexAndChunk>> futures;

		// preallocate chunks so we don't do it on an async thread
		for (chunkIndex.Z = fromChunkIndex.Z; chunkIndex.Z <= toChunkIndex.Z; chunkIndex.Z++)
		{
			for (chunkIndex.Y = fromChunkIndex.Y; chunkIndex.Y <= toChunkIndex.Y; chunkIndex.Y++)
			{
				for (chunkIndex.X = fromChunkIndex.X; chunkIndex.X <= toChunkIndex.X; chunkIndex.X++)
				{
					// gotta access this on this thread because it can create a chunk
					auto& chunk = grid.chunkAtChunkIndex(chunkIndex);
				}
			}
		}

		for (chunkIndex.Z = fromChunkIndex.Z; chunkIndex.Z <= toChunkIndex.Z; chunkIndex.Z++)
		{
			for (chunkIndex.Y = fromChunkIndex.Y; chunkIndex.Y <= toChunkIndex.Y; chunkIndex.Y++)
			{
				for (chunkIndex.X = fromChunkIndex.X; chunkIndex.X <= toChunkIndex.X; chunkIndex.X++)
				{
					// calculate the local cell index start and end
					FIntVector localFromCell;

					for (int c = 0; c < 3; ++c)
					{
						if (chunkIndex[c] == fromChunkIndex[c])
							localFromCell[c] = fromCellIndex[c];
						else
							localFromCell[c] = 0;
					}

					FIntVector localToCell;

					for (int c = 0; c < 3; ++c)
					{
						int maxDim = grid.chunkDimensions()[c];

						if (chunkIndex[c] == toChunkIndex[c])
							localToCell[c] = toCellIndex[c] >= maxDim ? maxDim - 1 : toCellIndex[c];
						else
							localToCell[c] = maxDim - 1;
					}

					futures.emplace_back(Async<IndexAndChunk>(EAsyncExecution::ThreadPool, [=,&grid]() mutable
					{
						IndexAndChunk chunkAndIndex;

						auto& chunk = grid.chunkAtChunkIndex(chunkIndex);

						chunkAndIndex.geometry = marchingCubes_byChunkAccessFast(isoLevel, grid, chunk, chunkIndex, localFromCell, localToCell, uvScale, reverseWinding, flatShading);
						chunkAndIndex.index = chunkIndex;

						return chunkAndIndex;
					}));
				}
			}
		}

		// join the results from the async tasks
		for (auto& future : futures)
		{
			const IndexAndChunk& chunkAndIndex = future.Get();

			result[chunkAndIndex.index] = std::move(chunkAndIndex.geometry);
		}
			
		return result;
	}

	/*
	Given a grid cell and an isolevel, calculate the triangular
	facets required to represent the isosurface through the cell.
	Return the number of triangular facets, the array "triangles"
	will be loaded up with the vertices at most 5 triangular facets.
	0 will be returned if the grid cell is either totally above
	of totally below the isolevel.
	*/
	static int polygonise(const GridCell& grid, float isolevel, Triangle * triangles /* array of size 5 */)
	{
		int i, ntriang;
		int cubeindex;
		Vertex vertlist[12];

		/*
		Determine the index into the edge table which
		tells us which vertices are inside of the surface
		*/
		cubeindex = 0;
		if (grid.v[0].W < isolevel) cubeindex |= 1;
		if (grid.v[1].W < isolevel) cubeindex |= 2;
		if (grid.v[2].W < isolevel) cubeindex |= 4;
		if (grid.v[3].W < isolevel) cubeindex |= 8;
		if (grid.v[4].W < isolevel) cubeindex |= 16;
		if (grid.v[5].W < isolevel) cubeindex |= 32;
		if (grid.v[6].W < isolevel) cubeindex |= 64;
		if (grid.v[7].W < isolevel) cubeindex |= 128;

		/* Cube is entirely in/out of the surface */
		if (edgeTable[cubeindex] == 0)
			return(0);

		/* Find the vertices where the surface intersects the cube */
		if ((edgeTable[cubeindex] & 1) != 0)
			vertlist[0] = vertexInterpolate(isolevel, grid.p[0], grid.p[1], grid.v[0], grid.v[1]);

		if ((edgeTable[cubeindex] & 2) != 0)
			vertlist[1] = vertexInterpolate(isolevel, grid.p[1], grid.p[2], grid.v[1], grid.v[2]);

		if ((edgeTable[cubeindex] & 4) != 0)
			vertlist[2] = vertexInterpolate(isolevel, grid.p[2], grid.p[3], grid.v[2], grid.v[3]);

		if ((edgeTable[cubeindex] & 8) != 0)
			vertlist[3] = vertexInterpolate(isolevel, grid.p[3], grid.p[0], grid.v[3], grid.v[0]);

		if ((edgeTable[cubeindex] & 16) != 0)
			vertlist[4] = vertexInterpolate(isolevel, grid.p[4], grid.p[5], grid.v[4], grid.v[5]);

		if ((edgeTable[cubeindex] & 32) != 0)
			vertlist[5] = vertexInterpolate(isolevel, grid.p[5], grid.p[6], grid.v[5], grid.v[6]);

		if ((edgeTable[cubeindex] & 64) != 0)
			vertlist[6] = vertexInterpolate(isolevel, grid.p[6], grid.p[7], grid.v[6], grid.v[7]);

		if ((edgeTable[cubeindex] & 128) != 0)
			vertlist[7] = vertexInterpolate(isolevel, grid.p[7], grid.p[4], grid.v[7], grid.v[4]);

		if ((edgeTable[cubeindex] & 256) != 0)
			vertlist[8] = vertexInterpolate(isolevel, grid.p[0], grid.p[4], grid.v[0], grid.v[4]);

		if ((edgeTable[cubeindex] & 512) != 0)
			vertlist[9] = vertexInterpolate(isolevel, grid.p[1], grid.p[5], grid.v[1], grid.v[5]);

		if ((edgeTable[cubeindex] & 1024) != 0)
			vertlist[10] = vertexInterpolate(isolevel, grid.p[2], grid.p[6], grid.v[2], grid.v[6]);

		if ((edgeTable[cubeindex] & 2048) != 0)
			vertlist[11] = vertexInterpolate(isolevel, grid.p[3], grid.p[7], grid.v[3], grid.v[7]);

		/* Create the triangle */
		ntriang = 0;
		for (i = 0; triTable[cubeindex][i] != -1; i += 3)
		{
			triangles[ntriang].p0 = vertlist[triTable[cubeindex][i]];
			triangles[ntriang].p1 = vertlist[triTable[cubeindex][i + 1]];
			triangles[ntriang].p2 = vertlist[triTable[cubeindex][i + 2]];
			ntriang++;
		}

		return(ntriang);
	}

	/*
	Linearly interpolate the position where an isosurface cuts
	an edge between two vertices, each with their own scalar value
	*/
	static Vertex vertexInterpolate(
		float isoLevel,
		const FVector& point0,
		const FVector& point1,
		const FVector4& value0,
		const FVector4& value1)
	{
		FVector normal;

		float mu = (isoLevel - value0.W) / (value1.W - value0.W);

		FVector p = point0 + mu * (point1 - point0);
		FVector4 normal4 = value0 + mu * (value1 - value0);

		// flip the normal
		normal = -1.0f * FVector(normal4.X, normal4.Y, normal4.Z);
		normal.Normalize();

		return{ p, normal };
	}

private:
	static const int8 triTable[256][16];

	// 256 elements
	static const uint16 edgeTable[256];
};

