// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "Components/InstancedStaticMeshComponent.h"


#include "MappedInstancedStaticMesh.generated.h"

class FPrimitiveSceneProxy;

/**
 * 
 */
UCLASS( ClassGroup = (ShipEditor), meta = (BlueprintSpawnableComponent) )
class LIFEBRUSH_API UMappedInstancedStaticMesh : public UInstancedStaticMeshComponent
{
	GENERATED_BODY()
public:
	/** Remove the instance specified. Returns True on success. The array will be out of order (calls RemoveSwapAt). 
	Code adapted from UHierarchicalInstancedStaticMesh.*/
	virtual bool RemoveInstance( int32 InstanceIndex ) override;
};

