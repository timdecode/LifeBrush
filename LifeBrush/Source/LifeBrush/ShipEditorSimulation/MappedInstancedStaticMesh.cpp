// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"

#include "Components/InstancedStaticMeshComponent.h"

#include "MappedInstancedStaticMesh.h"

bool UMappedInstancedStaticMesh::RemoveInstance( int32 InstanceIndex )
{


	// Save the render index
	const int32 RemovedRenderIndex = InstanceReorderTable[InstanceIndex];
	if(RemovedRenderIndex != INDEX_NONE)
	{
		RemovedInstances.Add( RemovedRenderIndex );
	}

	// Remove the instance
	PerInstanceSMData.RemoveAtSwap( InstanceIndex );
	InstanceReorderTable.RemoveAtSwap( InstanceIndex );

#if WITH_EDITOR
	if(SelectedInstances.Num())
	{
		if(GetWorld()->WorldType != EWorldType::Editor)
			SelectedInstances.RemoveAtSwap( InstanceIndex );
	}
#endif

	// update the physics state
	if(bPhysicsStateCreated)
	{

	}

	RemovedInstances.Reset();

	// Indicate we need to update render state to reflect changes
	MarkRenderStateDirty();

	return true;

}

