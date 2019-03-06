// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include <map>

struct InstanceManager
{
public:
	void clear()
	{
		for(auto& pair : elementTypesToInstancedMeshes)
			pair.second->ClearInstances();

		for(auto& available : _availableInstancesBuffer)
			available.clear();
	}

	void clearMeshTypes()
	{
		elementTypesToInstancedMeshes.clear();
	}

	void markDirty()
	{
		for(auto & pair : elementTypesToInstancedMeshes)
			pair.second->MarkRenderStateDirty();
	}

	auto meshForType( int16 type ) -> UInstancedStaticMeshComponent*
	{
		return elementTypesToInstancedMeshes[type];
	}

	void setMeshForType( UInstancedStaticMeshComponent * mesh, int16 type )
	{
		elementTypesToInstancedMeshes[type] = mesh;
	}

	void updateInstance( int32 instance, int16 type, FTransform& transform, bool worldSpace = true )
	{
		auto mesh = meshForType( type );

		if(mesh == nullptr)
			std::logic_error( "we need a mesh" );

		mesh->UpdateInstanceTransform( instance, transform, worldSpace );
	}

	void removeInstance( int32 instance, int16 type )
	{
		_nextAvailable[type].push_back( instance );

		hideInstance( instance, type );

		bob();
	}

	int32 addInstance( int16 type, FTransform& transform, bool worldSpace = true )
	{
		int32 instance;

		UInstancedStaticMeshComponent * mesh = elementTypesToInstancedMeshes[type];

		if (!mesh)
			return -1;

		if(mesh == nullptr)
			std::logic_error( "we need a mesh" );

		auto& available = _currentAvailable[type];
		// pull from the current available buffer
		if(available.size() > 0)
		{
			instance = available.back();
			available.pop_back();

			mesh->UpdateInstanceTransform( instance, transform, worldSpace );
		}
		else
		{
			if(worldSpace)
				instance = mesh->AddInstanceWorldSpace( transform );
			else
				instance = mesh->AddInstance( transform );
		}

		return instance;
	}

	void hideInstance( int32 instance, int16 type )
	{
		UInstancedStaticMeshComponent * mesh = elementTypesToInstancedMeshes[type];

		if(mesh == nullptr)
			std::logic_error( "we need a mesh" );

		FTransform t;
		mesh->GetInstanceTransform( instance, t, true );
		t.SetTranslation( FVector( 100000.0f, 0.0f, 0.0f ) );

		mesh->UpdateInstanceTransform( instance, t, true );
	}

	void rotateBuffers()
	{
		if(&_currentAvailable == &_availableInstancesBuffer[0])
		{
			_currentAvailable = _availableInstancesBuffer[1];
			_nextAvailable = _availableInstancesBuffer[0];
		}
		else
		{
			_currentAvailable = _availableInstancesBuffer[0];
			_nextAvailable = _availableInstancesBuffer[1];
		}
	}

	void setVisibility( bool visibile )
	{
		for(auto& pair : elementTypesToInstancedMeshes)
		{
			pair.second->SetVisibility( visibile );
		}
	}

	std::map<int16, UInstancedStaticMeshComponent*> elementTypesToInstancedMeshes;

protected:

	// a list of available instances in the instanced static meshes (due to deletion or change of type)
	typedef std::map<int16, std::vector<int32>> AvailableInstances;

	std::array< AvailableInstances, 2 > _availableInstancesBuffer;

	AvailableInstances& _currentAvailable = _availableInstancesBuffer[0];
	AvailableInstances& _nextAvailable = _availableInstancesBuffer[0];
private:

	void bob()
	{
	}


};