// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include <vector>
#include <cinttypes>
#include <cassert>

#include "Pool.generated.h"

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FPoolBase
{
	GENERATED_BODY()

protected:

public:
	UPROPERTY()

	uint32 _size = 0;           // number of elements we have

	UPROPERTY()

	TArray<uint8> _data;

public:
	FPoolBase() {}
	
	FPoolBase& operator=( const FPoolBase& other )
	{
		_size = other._size; 
		_data = TArray<uint8>( other._data );

		return *this;
	}

	bool operator==( const FPoolBase& other ) const 
	{
		return (_size == other._size) && (_data == other._data);
	}

	void resize( size_t nElements, size_t elementSize )
	{
		_data.SetNumUninitialized( nElements * elementSize, true );

		_size = nElements;
	}

	void* get( size_t i, size_t elementSize )
	{
		assert( i < _size );

		return &(_data[i * elementSize]);
	}

	void* begin( size_t elementSize )
	{
		return &(_data.GetData()[0]);
	}

	void* end(size_t elementSize)
	{
		return &(_data.GetData()[_size * elementSize]);
	}
};

template<class TElement>
class Pool : FPoolBase
{
public:
	// won't call destructor on TElement
	void resize( size_t nElements )
	{
		FPoolBase::resize( nElements, sizeof( TElement ) );
	}

	TElement* get( size_t i )
	{
		return static_cast<TElement*>(FPoolBase::get( i, sizeof( TElement ) ));
	}

	TElement* begin()
	{
		return static_cast<TElement*>(FPoolBase::begin( sizeof( TElement ) ));
	}

	TElement* end()
	{
		return static_cast<TElement*>(FPoolBase::end( sizeof( TElement ) ));
	}

	void destroy( size_t i )
	{
		assert( i < _size );

		TElement * element = static_cast<TElement*>get( i );
		element->~TElement();
	}

	TElement& operator[](size_t i)
	{
		return *get(i);
	}
};



