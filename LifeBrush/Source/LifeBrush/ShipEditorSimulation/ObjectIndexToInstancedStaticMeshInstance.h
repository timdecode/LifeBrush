// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

struct ObjectIndexToInstancedStaticMeshInstance
{
	struct ObjectIndex
	{
		ObjectIndex( int32 objectIndex_in ) : nodeIndex( objectIndex_in ) {}

		int32 nodeIndex;
	};

	TArray<ObjectIndex> _objectIndices;
	MappedArray<ObjectIndex> _map;

	ObjectIndexToInstancedStaticMeshInstance()
	{
		_map.init( _objectIndices );
	}

	void emplace( int32 objectIndex, int32 instanceIndex )
	{
		_map.emplace( objectIndex, instanceIndex );
	}

	void erase( int32 objectIndex )
	{
		_map.erase( objectIndex );
	}

	int32 instanceAt( int32 objectIndex )
	{
		return _map.elementIndex( objectIndex );
	}

	int32 objectAt( int32 instanceIndex )
	{
		return _objectIndices[instanceIndex].nodeIndex;
	}

	bool contains( int32 objectIndex )
	{
		return _map.contains( objectIndex );
	}

	void clear()
	{
		_map.clear();
		_objectIndices.Empty();
	}
};