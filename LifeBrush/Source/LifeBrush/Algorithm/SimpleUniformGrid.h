// Copyright 2018 Code Monkey Castle, all rights reserved.

#pragma once

#include <vector>

#include "Utility.h"

template<typename ElementType>
class SimpleUniformGrid
{
protected:
	Eigen::AlignedBox3f _worldBounds;

	// cached properties of the world bounds
	Eigen::Vector3f _worldMin = Eigen::Vector3f::Zero();
	Eigen::Vector3f _worldSize = Eigen::Vector3f::Zero();

	Eigen::Vector3f _cellSize = Eigen::Vector3f::Zero();
	Eigen::Vector3f _invCellSize = Eigen::Vector3f::Zero();

	Eigen::Vector3i _numSamples = Eigen::Vector3i::Zero();

protected:
	void _setup( Eigen::Vector3i dimensions, Eigen::Vector3f cellSize, Eigen::Vector3f worldMin )
	{
		_worldMin = worldMin;

		for(int c = 0; c < 3; ++c)
			_worldSize[c] = dimensions[c] * cellSize[c];

		_cellSize = cellSize;
		_invCellSize = Eigen::Vector3f( 1.0f / cellSize.x(), 1.0f / cellSize.y(), 1.0f / cellSize.z() );
		_numSamples = dimensions;

		_worldBounds = AlignedBox3f( _worldMin, _worldMin + _worldSize );
	}

public:
	Eigen::Vector3i numSamples() { return _numSamples; }

	Eigen::Vector3f cellSize() { return _cellSize; }
	FBox bounds() { return _worldBounds; }

	std::vector<ElementType> samples;

	void init( Eigen::Vector3i dimensions, Eigen::Vector3f cellSize, Eigen::Vector3f worldMin )
	{
		_setup( dimensions, cellSize, worldMin );

		samples.resize( _numSamples.x() * _numSamples.y() * _numSamples.z() );
	}

	// ---------------------------------------------------
	// Sample Data Accessors
	// ---------------------------------------------------
	ElementType& operator()( size_t x, size_t y, size_t z )
	{
		return samples[x + y * _numSamples.x() + z * _numSamples.x() * _numSamples.y()];
	}

	ElementType& operator()( const Eigen::Vector3i& i )
	{
		return samples[i.x() + i.y() * _numSamples.x() + i.z() * _numSamples.x() * _numSamples.y()];
	}

	ElementType& operator()( size_t i )
	{
		return samples[i];
	}


	ElementType& clamped( int x, int y, int z )
	{
		x = x < 0 ? 0 : x >= _numSamples.x() ? _numSamples.x() - 1 : x;
		y = y < 0 ? 0 : y >= _numSamples.y() ? _numSamples.y() - 1 : y;
		z = z < 0 ? 0 : z >= _numSamples.z() ? _numSamples.z() - 1 : z;

		return samples[x + y * _numSamples.x() + z * _numSamples.x() * _numSamples.y()];
	}

	ElementType& clamped( const Eigen::Vector3i& i )
	{
		const Eigen::Vector3i ci = clampedIndex( i );

		return samples[ci.x() + ci.y() * _numSamples.x() + ci.z() * _numSamples.x() * _numSamples.y()];
	}

	Eigen::Vector3i clampedIndex( const Eigen::Vector3i& index )
	{
		Eigen::Vector3i clamped;

		clamped.x() = FMath::Clamp( index.x(), 0, _numSamples.x() - 1 );
		clamped.y() = FMath::Clamp( index.y(), 0, _numSamples.y() - 1 );
		clamped.z() = FMath::Clamp( index.z(), 0, _numSamples.z() - 1 );

		return clamped;
	}

	// ---------------------------------------------------
	// Point Transforms
	// ---------------------------------------------------
	Eigen::Vector3f toGridSpace( const Eigen::Vector3f& worldPoint )
	{
		return (worldPoint - _worldMin).cwiseProduct(_invCellSize);
	}

	// The relative position within a cell of the world point
	Eigen::Vector3f toReferenceSpace( const Eigen::Vector3f& worldPoint )
	{
		Eigen::Vector3f gridPoint = toGridSpace( worldPoint );

		Eigen::Vector3f cellPoint;
		for(int c = 0; c < 3; ++c)
			cellPoint[c] = FMath::FloorToFloat( gridPoint[c] );

		Eigen::Vector3f r = gridPoint - cellPoint;

		return r;
	}

	Eigen::Vector3i sampleIndex( const Eigen::Vector3f& worldPoint )
	{
		return toGridSpace( worldPoint ).cast<int>();
	}

	Eigen::Vector3f samplePoint( const Eigen::Vector3i& index )
	{
		return index.cast<float>().cwiseProduct( _cellSize ) + _worldMin;
	}

	int flatten( const Eigen::Vector3i& i )
	{
		return i.x() + i.y() * _numSamples.x() + i.z() * _numSamples.x() * _numSamples.y();
	}

	Eigen::Vector3i inflate( int i )
	{
		Eigen::Vector3i vec;

		const int nLayer = _numSamples.x() * _numSamples.y();
		const int remainder = i % nLayer;

		vec.z() = i / (nLayer);
		vec.y() = remainder / _numSamples.x();
		vec.x() = remainder % _numSamples.x();

		return vec;
	}
};