// Copyright 2018 Code Monkey Castle, all rights reserved.

#pragma once

#include <vector>
#include <unordered_map>
#include <array>
#include <memory>

#include "Utility.h"

namespace lb
{
	template<typename ElementType>
	class BasicGrid
	{
	protected:
		std::vector<ElementType> samples;

		FIntVector _dimensions;

	public:
		void init(size_t xn, size_t yn, size_t zn)
		{
			_dimensions = FIntVector(xn, yn, zn);

			samples.resize(_dimensions.X * _dimensions.Y * _dimensions.Z);
		}

		FIntVector dimensions() { return _dimensions; }


		ElementType& operator()(size_t x, size_t y, size_t z)
		{
			return samples[x + y * _dimensions.X + z * _dimensions.X * _dimensions.Y];
		}

		ElementType& operator()(const FIntVector& i)
		{
			return samples[i.X + i.Y * _dimensions.X + i.Z * _dimensions.X * _dimensions.Y];
		}

		ElementType& operator()(size_t i)
		{
			return samples[i];
		}
	};
};




template<typename ElementType>
class UniformGrid
{
protected:
	FBox _worldBounds;

	// cached properties of the world bounds
	FVector _worldMin = FVector::ZeroVector;
	FVector _worldSize = FVector::ZeroVector;

	FVector _cellSize = FVector::ZeroVector;
	FVector _invCellSize = FVector::ZeroVector;

	FIntVector _dimensions = FIntVector::ZeroValue;

protected:
	void _setup( FIntVector dimensions, FVector cellSize, FVector worldMin )
	{
		_worldMin = worldMin;

		for(int c = 0; c < 3; ++c)
			_worldSize[c] = dimensions[c] * cellSize[c];

		_cellSize = cellSize;
		_invCellSize = FVector( 1.0f / cellSize.X, 1.0f / cellSize.Y, 1.0f / cellSize.Z );
		_dimensions = dimensions;

		_worldBounds = FBox( _worldMin, _worldMin + _worldSize );
	}

public:
	FIntVector dimensions() { return _dimensions; }

	FVector cellSize() { return _cellSize; }
	FBox bounds() { return _worldBounds; }

	std::vector<ElementType> samples;

	static const FIntVector vertexIndices[8];

	void init( FIntVector dimensions, FVector cellSize, FVector worldMin, const std::vector<ElementType>& samples_in )
	{
		_setup( dimensions, cellSize, worldMin );

		if(samples_in.size() == 0)
			samples.resize( _dimensions.X * _dimensions.Y * _dimensions.Z );
		else
			samples = samples_in;
	}

	// ---------------------------------------------------
	// Sample Data Accessors
	// ---------------------------------------------------
	ElementType& operator()( size_t x, size_t y, size_t z )
	{
		return samples[x + y * _dimensions.X + z * _dimensions.X * _dimensions.Y];
	}

	ElementType& operator()( const FIntVector& i )
	{
		return samples[i.X + i.Y * _dimensions.X + i.Z * _dimensions.X * _dimensions.Y];
	}

	ElementType& operator()( size_t i )
	{
		return samples[i];
	}


	ElementType& clamped( int x, int y, int z )
	{
		x = x < 0 ? 0 : x >= _dimensions.X ? _dimensions.X - 1 : x;
		y = y < 0 ? 0 : y >= _dimensions.Y ? _dimensions.Y - 1 : y;
		z = z < 0 ? 0 : z >= _dimensions.Z ? _dimensions.Z - 1 : z;

		return samples[x + y * _dimensions.X + z * _dimensions.X * _dimensions.Y];
	}

	ElementType& clamped( const FIntVector& i )
	{
		const FIntVector ci = clampedIndex( i );

		return samples[ci.X + ci.Y * _dimensions.X + ci.Z * _dimensions.X * _dimensions.Y];
	}

	ElementType& clampedToZero(int x, int y, int z)
	{
		if (x < 0 || y < 0 || z < 0)
			return 0;
		else if (x >= _dimensions.X || y >= _dimensions.Y || z >= _dimensions.Z)
			return 0;
		else
			return samples[x + y * _dimensions.X + z * _dimensions.X * _dimensions.Y];
	}

	ElementType& clampedToZero(const FIntVector& i)
	{
		return clampedToZero(i.X, i.Y, i.Z);
	}

	FIntVector clampedIndex( const FIntVector& index )
	{
		FIntVector clamped;

		clamped.X = FMath::Clamp( index.X, 0, _dimensions.X - 1 );
		clamped.Y = FMath::Clamp( index.Y, 0, _dimensions.Y - 1 );
		clamped.Z = FMath::Clamp( index.Z, 0, _dimensions.Z - 1 );

		return clamped;
	}

	// ---------------------------------------------------
	// Point Transforms
	// ---------------------------------------------------
	FVector toGridSpace( const FVector& worldPoint )
	{
		return (worldPoint - _worldMin) * _invCellSize;
	}

	FVector toReferenceSpace( const FVector& worldPoint )
	{
		FVector gridPoint = toGridSpace( worldPoint );

		FVector cellPoint;
		for(int c = 0; c < 3; ++c)
			cellPoint[c] = FMath::FloorToFloat( gridPoint[c] );

		FVector r = gridPoint - cellPoint;

		return r;
	}

	FIntVector sampleIndex( const FVector& worldPoint )
	{
		FVector gridSpace = toGridSpace( worldPoint );
		FIntVector index;
		for(int c = 0; c < 3; ++c)
			index[c] = int( gridSpace[c] );

		return index;
	}

	FVector samplePoint( const FIntVector& index )
	{
		return Utility::scale( index, _cellSize ) + _worldMin;
	}

	int flatten( const FIntVector& i )
	{
		return i.X + i.Y * _dimensions.X + i.Z * _dimensions.X * _dimensions.Y;
	}

	FIntVector inflate( int i )
	{
		FIntVector vec;

		const int nLayer = _dimensions.X * _dimensions.Y;
		const int remainder = i % nLayer;

		vec.Z = i / (nLayer);
		vec.Y = remainder / _dimensions.X;
		vec.X = remainder % _dimensions.X;

		return vec;
	}


	//	// ---------------------------------------------------
	//	// Interpolation
	//	// ---------------------------------------------------
	float trilinearValueAt( FVector worldPoint )
	{
		FIntVector pi = sampleIndex( worldPoint );

		FVector r = toReferenceSpace( worldPoint );

		float c_00 = clamped( pi.X + 0, pi.Y + 0, pi.Z + 0 ) * (1.0f - r.X) + clamped( pi.X + 1, pi.Y + 0, pi.Z + 0 ) * r.X;
		float c_10 = clamped( pi.X + 0, pi.Y + 1, pi.Z + 0 ) * (1.0f - r.X) + clamped( pi.X + 1, pi.Y + 1, pi.Z + 0 ) * r.X;
		float c_01 = clamped( pi.X + 0, pi.Y + 0, pi.Z + 1 ) * (1.0f - r.X) + clamped( pi.X + 1, pi.Y + 0, pi.Z + 1 ) * r.X;
		float c_11 = clamped( pi.X + 0, pi.Y + 1, pi.Z + 1 ) * (1.0f - r.X) + clamped( pi.X + 1, pi.Y + 1, pi.Z + 1 ) * r.X;

		float c_0 = c_00 * (1.0f - r.Y) + c_10 * r.Y;
		float c_1 = c_01 * (1.0f - r.Y) + c_11 * r.Y;

		float c = c_0 * (1.0f - r.Z) + c_1 * r.Z;

		return c;
	}

	FVector normalAt( FVector worldPoint )
	{
		FVector n;

		FVector dx = FVector::ZeroVector;
		dx.X = _cellSize.X;

		FVector dy = FVector::ZeroVector;
		dy.Y = _cellSize.Y;

		FVector dz = FVector::ZeroVector;
		dz.Z = _cellSize.Z;

		n.X = trilinearValueAt( worldPoint - dx ) - trilinearValueAt( worldPoint + dx );
		n.Y = trilinearValueAt( worldPoint - dy ) - trilinearValueAt( worldPoint + dy );
		n.Z = trilinearValueAt( worldPoint - dz ) - trilinearValueAt( worldPoint + dz );

		return n.GetSafeNormal();
	}

	FVector normalAt( const FIntVector& i )
	{
		FVector n;

		n.X = clamped( i.X - 1, i.Y, i.Z ) - clamped( i.X + 1, i.Y, i.Z );
		n.Y = clamped( i.X, i.Y - 1, i.Z ) - clamped( i.X, i.Y + 1, i.Z );
		n.Z = clamped( i.X, i.Y, i.Z - 1 ) - clamped( i.X, i.Y, i.Z + 1 );

		return n.GetSafeNormal();
	}
};

template<typename ElementType>
class PaddedUniformGrid
{
protected:
	FBox _worldBounds;

	// cached properties of the world bounds
	FVector _worldMin = FVector::ZeroVector;
	FVector _worldSize = FVector::ZeroVector;

	FVector _cellSize = FVector::ZeroVector;
	FVector _invCellSize = FVector::ZeroVector;

	FIntVector _paddedDimensions = FIntVector::ZeroValue;

	FIntVector _padding = FIntVector(2, 2, 2);

protected:
	void _setup(FIntVector dimensions, FVector cellSize, FVector worldMin)
	{
		_worldMin = worldMin;

		for (int c = 0; c < 3; ++c)
			_worldSize[c] = dimensions[c] * cellSize[c];

		_cellSize = cellSize;
		_invCellSize = FVector(1.0f / cellSize.X, 1.0f / cellSize.Y, 1.0f / cellSize.Z);
		_paddedDimensions = dimensions + _padding + _padding;

		_worldBounds = FBox(_worldMin, _worldMin + _worldSize);
	}

public:
	FIntVector dimensions() { return _paddedDimensions - (_padding + _padding); }

	const FIntVector& padding() { return _padding; }

	FVector cellSize() { return _cellSize; }
	FBox bounds() { return _worldBounds; }

	std::vector<ElementType> samples;

	static const FIntVector vertexIndices[8];

	void init(FIntVector dimensions, FVector cellSize, FVector worldMin)
	{
		_setup(dimensions, cellSize, worldMin);

		samples.resize(_paddedDimensions.X * _paddedDimensions.Y * _paddedDimensions.Z);
	}

	// ---------------------------------------------------
	// Sample Data Accessors
	// ---------------------------------------------------
	ElementType& operator()(size_t x, size_t y, size_t z)
	{
		x += _padding.X;
		y += _padding.Y;
		z += _padding.Z;

		return samples[x + y * _paddedDimensions.X + z * _paddedDimensions.X * _paddedDimensions.Y];
	}

	ElementType& operator()(FIntVector i)
	{
		i += _padding;

		return samples[i.X + i.Y * _paddedDimensions.X + i.Z * _paddedDimensions.X * _paddedDimensions.Y];
	}

	// ---------------------------------------------------
	// Point Transforms
	// ---------------------------------------------------
	FVector toGridSpace(const FVector& worldPoint)
	{
		return (worldPoint - _worldMin) * _invCellSize;
	}

	FVector toReferenceSpace(const FVector& worldPoint)
	{
		FVector gridPoint = toGridSpace(worldPoint);

		FVector cellPoint;
		for (int c = 0; c < 3; ++c)
			cellPoint[c] = FMath::FloorToFloat(gridPoint[c]);

		FVector r = gridPoint - cellPoint;

		return r;
	}

	FIntVector sampleIndex(const FVector& worldPoint)
	{
		FVector gridSpace = toGridSpace(worldPoint);
		FIntVector index;
		for (int c = 0; c < 3; ++c)
			index[c] = int(gridSpace[c]);

		return index;
	}

	FVector samplePoint(const FIntVector& index)
	{
		return Utility::scale(index, _cellSize) + _worldMin;
	}

	int flatten(const FIntVector& i)
	{
		const FIntVector dim = _paddedDimensions - _padding;

		return i.X + i.Y * _paddedDimensions.X + i.Z * _paddedDimensions.X * _paddedDimensions.Y;
	}

	FIntVector inflate(int i)
	{
		const FIntVector dim = _paddedDimensions - _padding;

		FIntVector vec;

		const int nLayer = dim.X * dim.Y;
		const int remainder = i % nLayer;

		vec.Z = i / (nLayer);
		vec.Y = remainder / dim.X;
		vec.X = remainder % dim.X;

		return vec;
	}

	FVector normalAt(const FIntVector& i)
	{
		FVector n;

		n.X = (*this)(i.X - 1, i.Y, i.Z) - (*this)(i.X + 1, i.Y, i.Z);
		n.Y = (*this)(i.X, i.Y - 1, i.Z) - (*this)(i.X, i.Y + 1, i.Z);
		n.Z = (*this)(i.X, i.Y, i.Z - 1) - (*this)(i.X, i.Y, i.Z + 1);

		return n.GetSafeNormal();
	}
};

template<typename ElementType>
class ChunkGrid
{
protected:
	FIntVector _chunkDimensions; // the number of columns, rows and slices in a chunk

	FVector _cellSize; // the world-size of a cell 
	FVector _invCellSize;

	FVector _chunkSize; // the world-size of a chunk

public:
	FIntVector chunkDimensions() { return _chunkDimensions; }
	FVector cellSize() { return _cellSize; }
	FVector chunkSize() { return _chunkSize; }

public:
	void init(FIntVector chunkDimensions, FVector cellSize)
	{
		this->_chunkDimensions = chunkDimensions;
		this->_cellSize = cellSize;

		for (int c = 0; c < 3; ++c)
			_chunkSize[c] = chunkDimensions[c] * cellSize[c];

		for (int c = 0; c < 3; ++c)
			_invCellSize[c] = 1.0f / cellSize[c];

		_chunks.Empty();
	}

	bool didInit()
	{
		return chunkDimensions().X > 0;
	}

	FBox bounds()
	{
		FVector min(std::numeric_limits<float>::max());
		FVector max(std::numeric_limits<float>::min());

		for (auto& pair : _chunks)
		{
			FIntVector chunk = pair.Key;

			FVector chunkStart = gridIndexFromChunkIndex(chunk);
			FVector chunkEnd = chunkStart + _chunkSize;

			for (int c = 0; c < 3; ++c)
			{
				if (chunkStart[c] < min[c]) min[c] = chunkStart[c];
				if (chunkEnd[c] > max[c]) max[c] = chunkEnd[c];
			}
		}

		return FBox(min, max);
	}

	FIntVector minChunkndex()
	{
		FIntVector min(std::numeric_limits<int>::max());

		for (auto& pair : _chunks)
		{
			FIntVector index = pair.Key;

			for (int c = 0; c < 3; ++c)
			{
				if (index[c] < min[c]) min[c] = index[c];
			}
		}

		return min;
	}

	FIntVector maxChunkIndex()
	{
		FIntVector max(std::numeric_limits<int>::min());

		for (auto& pair : _chunks)
		{
			FIntVector index = pair.Key;

			for (int c = 0; c < 3; ++c)
			{
				if (index[c] > max[c]) max[c] = index[c];
			}
		}

		return max;
	}

	// Will return positive remainders. Examples:
	// mod_negative(-5, 64) => 59
	// mod_negative( 5, 64) => 5
	static int mod_forChunks(int dividend, int divisor)
	{
		int r = std::abs(dividend) % divisor;

		if (dividend < 0 && r != 0)
			r = divisor - r;

		return r;
	}

	// will also overwrite an existing chunk
	PaddedUniformGrid<float>& _createChunkAt(FIntVector chunkIndex)
	{
		_chunks.Add(chunkIndex, std::make_unique<PaddedUniformGrid<float>>());

		PaddedUniformGrid<float>& chunk = *_chunks[chunkIndex].get();

		FVector base;
		for (int c = 0; c < 3; ++c)
		{
			base[c] = _chunkDimensions[c] * chunkIndex[c] * _cellSize[c];
		}

		chunk.init(_chunkDimensions, _cellSize, base);


		FIntVector chunkGridBase = gridIndexFromChunkIndex(chunkIndex);

		// reordering array for the planes component indices
		const FIntVector planeComponents[6] = {
			FIntVector(0,1,2),
			FIntVector(0,1,2),
			FIntVector(1,0,2),
			FIntVector(1,0,2),
			FIntVector(2,0,1),
			FIntVector(2,0,1)
		};

		// copy borders
		// idea: the control what plane we are enumerate with the components array
		//       the indices in the array are the components of the plane to access
		int planeIndex = 0;
		for (const FIntVector& componentVec = planeComponents[planeIndex]; planeIndex < 6; ++planeIndex)
		{
			FIntVector index;

			const std::array<int, 2> planeBaseIndices = { {-1, _chunkDimensions[0]} };

			// choose the plane base (offset along the normal)
			for (int baseIndex : planeBaseIndices)
			{
				index[componentVec[0]] = baseIndex;

				// use i, the plane index. if even, we'll do the negative normal, else the positive normal
				FIntVector chunkOffset = FIntVector::ZeroValue;
				chunkOffset[componentVec[0]] = planeIndex % 2 == 0 ? -1 : 1;

				FIntVector neighbourChunkIndex = chunkIndex + chunkOffset;

				// if we don't have a neighbour, there is nothing to copy!
				if(!_chunks.Contains(neighbourChunkIndex)) continue;

				FIntVector neighbourGridBase = gridIndexFromChunkIndex(neighbourChunkIndex);


				PaddedUniformGrid<float>& neighbour = *_chunks[neighbourChunkIndex].get();

				// rows of the plane (from the non-border values)
				for (index[componentVec[1]] = 0; index[componentVec[1]] < _chunkDimensions[componentVec[1]]; index[componentVec[1]]++)
				{
					// columns of the plane
					for (index[componentVec[2]] = 0; index[componentVec[2]] < _chunkDimensions[componentVec[2]]; index[componentVec[2]]++)
					{

						FIntVector neighbourIndex = chunkGridBase + index - neighbourGridBase;

						chunk(index) = neighbour(neighbourIndex);
					}
				}
			}
		}

		return chunk;
	}

	PaddedUniformGrid<ElementType>& chunkAtChunkIndex(FIntVector chunkIndex)
	{
		// do we need a new chunk?
		auto found = _chunks.Find(chunkIndex);

		// not found
		if (!found)
		{
			// we need a new chunk
			return _createChunkAt(chunkIndex);
		}
		// found
		else
		{
			return *found->get();
		}
	}

	PaddedUniformGrid<ElementType>& chunkAtGridIndex(FIntVector gridIndex)
	{
		FIntVector ci = indexOfChunk(gridIndex);

		return chunkAtChunkIndex(ci);
	}

	FVector normalAt(const FIntVector& i)
	{
		FVector n;

		n.X = (*this)(i.X - 1, i.Y, i.Z) - (*this)(i.X + 1, i.Y, i.Z);
		n.Y = (*this)(i.X, i.Y - 1, i.Z) - (*this)(i.X, i.Y + 1, i.Z);
		n.Z = (*this)(i.X, i.Y, i.Z - 1) - (*this)(i.X, i.Y, i.Z + 1);

		return n.GetSafeNormal();
	}

	FIntVector sampleIndex(const FVector& samplePoint)
	{
		return FIntVector(Utility::scale(samplePoint, _invCellSize));
	}

	FVector samplePoint(const FIntVector& index)
	{
		return Utility::scale(index, _cellSize);
	}

	// can contain negative indices.
	// origin is at 0,0,0
	__declspec(noinline)
	ElementType operator()(const FIntVector& gridIndex)
	{
		FIntVector chunkIndex = indexOfChunk(gridIndex);

		// do we need a new chunk?
		auto found = _chunks.Find(chunkIndex);

		// not found
		if (!found )
		{
			// fake it
			return 0.0f;
		}
		// found
		else
		{
			PaddedUniformGrid<ElementType>& grid = *found->get();

			FIntVector ci = indexOfCell(gridIndex);

			return grid(ci);
		}
	}

	__declspec(noinline)
	ElementType operator()(const int x, const int y, const int z)
	{
		FIntVector gridIndex(x, y, z);

		return (*this)(gridIndex);
	}

	void set(const FIntVector& gridIndex, ElementType value)
	{
		FIntVector cellIndex = indexOfCell(gridIndex);
		FIntVector chunkIndex = indexOfChunk(gridIndex);

		PaddedUniformGrid<ElementType>& grid = chunkAtChunkIndex(chunkIndex);

		// luckily the border offset will only access up to 3 adjacent planes of the neighbouring grid
		FIntVector borderOffsets[3] = { FIntVector::ZeroValue, FIntVector::ZeroValue, FIntVector::ZeroValue };

		int n = 0;

		const int padding = grid.padding()[0];

		// if we are in a padded region, set a border offset
		for (int c = 0; c < 3; ++c)
		{
			if (cellIndex[c] < padding)
				borderOffsets[n++][c] = -1;
			else if (cellIndex[c] >= _chunkDimensions[c] - padding)
				borderOffsets[n++][c] = 1;
		}

		// we have to set this grid value before we set borders as there could be copying
		// involved if we create a new grid as the result of a chunkAt call

		grid(cellIndex) = value;


		for (int j = 0; j < n; ++j)
		{
			FIntVector borderOffset = borderOffsets[j];


			for (int p = 0; p < padding; p++)
			{
				FIntVector gi = gridIndex + borderOffset * float(p);

				FIntVector chunki = indexOfChunk(gi);

				PaddedUniformGrid<ElementType>& borderGrid = chunkAtGridIndex(gi);

				FIntVector ciOffset = indexOfCellRelativeToChunk(gridIndex, chunki);

				borderGrid(ciOffset) = value;
			}


		}


	}

	FIntVector inflateCellIndex(int i)
	{
		FIntVector vec;

		const int nLayer = _chunkDimensions.X * _chunkDimensions.Y;
		const int remainder = i % nLayer;

		vec.Z = i / (nLayer);
		vec.Y = remainder / _chunkDimensions.X;
		vec.X = remainder % _chunkDimensions.X;

		return vec;
	}

	FIntVector indexOfChunk(FIntVector gridIndex)
	{
		FIntVector index;
		for (int c = 0; c < 3; ++c)
		{
			index[c] = std::floor(float(gridIndex[c]) / float(_chunkDimensions[c]));
		}

		return index;
	}

	FIntVector gridIndexFromChunkIndex(FIntVector chunkIndex)
	{
		FIntVector gridIndex;

		for( int c = 0; c < 3; ++c )
		{
			gridIndex[c] = chunkIndex[c] * _chunkDimensions[c];
		}

		return gridIndex;
	}

	FIntVector indexOfCell(FIntVector gridIndex)
	{
		FIntVector index;

		FIntVector chunkIndex;

		for (int c = 0; c < 3; ++c)
		{
			index[c] = mod_forChunks(gridIndex[c], _chunkDimensions[c]);
		}

		return index;
	}

	FIntVector indexOfCellRelativeToChunk(FIntVector gridIndex_in, FIntVector chunkIndex_in)
	{
		FIntVector chunkGridStart = gridIndexFromChunkIndex(chunkIndex_in);

		return gridIndex_in - chunkGridStart;
	}

	template<typename Lambda>
	void enumerate( FIntVector from, FIntVector to, Lambda& enumerator)
	{
		FIntVector fromChunkIndex = chunkIndex(from);
		FIntVector toChunkIndex = chunkIndex(to);

		FIntVector fromCellIndex = indexOfCell(from);
		FIntVector toCellIndex = indexOfCell(to);

		FIntVector chunkIndex = fromChunkIndex;

		for (chunkIndex.Z = fromChunkIndex.Z; chunkIndex.Z <= toChunkIndex.Z; chunkIndex.Z++)
		{
			for (chunkIndex.Y = fromChunkIndex.Y; chunkIndex.Y <= toChunkIndex.Y; chunkIndex.Y++)
			{
				for (chunkIndex.X = fromChunkIndex.X; chunkIndex.X <= toChunkIndex.X; chunkIndex.X++)
				{
					// we'll enumerate the relative indices inside of each chunk
					PaddedUniformGrid<ElementType>& chunk = chunkAtChunkIndex(chunkIndex);

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
						if (chunkIndex[c] == toChunkIndex[c])
							localToCell = toCellIndex[c];
						else
							localToCell = _chunkDimensions[c];
					}

					FIntVector cellIndex;

					for (cellIndex.Z = localFromCell.Z; cellIndex.Z <= localToCell.Z; cellIndex.Z++)
					{
						for (cellIndex.Y = localFromCell.Y; cellIndex.Y <= localToCell.Y; cellIndex.Y++)
						{
							for (cellIndex.X = localFromCell.X; cellIndex.X <= localToCell.X; cellIndex.X++)
							{
								enumerator(cellIndex, chunkIndex, chunk);
							}
						}
					}


				}
			}
		}
	}


protected:
	// actually, let not use FVector as in the index. we'll use an int vector
	/*
	std::size_t operator()(std::vector<uint32_t> const& vec) const {
  std::size_t seed = vec.size();
  for(auto& i : vec) {
	seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
  return seed;
}*/
	TMap< FIntVector, std::unique_ptr<PaddedUniformGrid<float>> > _chunks;
};