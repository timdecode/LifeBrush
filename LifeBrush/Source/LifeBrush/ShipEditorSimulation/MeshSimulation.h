// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "Graph.h"
#include "ObjectSimulation.h"

#include <unordered_map>
#include "ObjectIndexToInstancedStaticMeshInstance.h"

#include "MeshSimulation.generated.h"

class UInstancedStaticMeshComponent;

struct FGraphMeshKey
{
	UStaticMesh * staticMesh;
	UMaterialInterface * material;

	bool operator==( const FGraphMeshKey& other ) const
	{
		return staticMesh == other.staticMesh && material == other.material;
	}
};



USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraphMesh : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UStaticMesh * staticMesh = nullptr;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UMaterialInterface * material = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	bool visible = true;

	operator FGraphMeshKey() { return FGraphMeshKey{ staticMesh, material }; }

	UInstancedStaticMeshComponent * _transientInstanceMesh = nullptr;

	void markDirty()
	{
		_transientInstanceMesh = nullptr;
	}
};



namespace std
{
	template<> struct hash<FGraphMeshKey>
	{
		std::size_t operator()( const FGraphMeshKey& key ) const
		{
			using std::hash;

			return hash<UStaticMesh*>()(key.staticMesh) ^ hash<UMaterialInterface*>()(key.material);
		}
	};
}

/**
 * 
 */
UCLASS()
class LIFEBRUSH_API UMeshSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	// ---------------------------------------------
	// UMeshSimulation
	// ---------------------------------------------
	virtual void attach() override;
	virtual void detach() override;

	virtual void begin() override;

	virtual void tick(float deltaT) override;

	virtual void tick_paused(float deltaT) override;

	virtual bool canRunInEditor() override;

	virtual void snapshotToActor(AActor * actor) override;

	void updateInstances();

	TArray<UInstancedStaticMeshComponent*> instancedStaticMeshes();

	FGraphNodeHandle nodeForInstance(UInstancedStaticMeshComponent * mesh, int32 instance);

	FBox boundsBoxForNode(FGraphNodeHandle node);

	// ---------------------------------------------
	// Component Listener Interface
	// ---------------------------------------------
	virtual void componentAdded( FGraphNodeHandle nodeHandle, ComponentType type ) override;

	virtual void componentRemoved( FGraphNodeHandle node, ComponentType type ) override;

	virtual void componentUpdated( FGraphNodeHandle node, ComponentType type ) override;


protected:
	UInstancedStaticMeshComponent * getOrCreateInstancedMesh( FGraphMesh& graphMesh );

	void _clearTransients();

	std::unordered_map<FGraphMeshKey, UInstancedStaticMeshComponent*> _instancedStaticeMeshes;
	std::unordered_map<UInstancedStaticMeshComponent*, std::vector<FGraphNodeHandle>> _instanceToNodeHandle;
};
