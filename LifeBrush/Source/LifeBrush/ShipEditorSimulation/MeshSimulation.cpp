// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"

#include "Components/HierarchicalInstancedStaticMeshComponent.h"

#include "MeshSimulation.h"

void UMeshSimulation::attach()
{
	_clearTransients();
}

void UMeshSimulation::begin()
{
	_clearTransients();
}

void UMeshSimulation::detach()
{
	// nuke the ISMCs
	for (auto& pair : _instancedStaticeMeshes)
	{
		UInstancedStaticMeshComponent * mesh = pair.second;

		mesh->ClearInstances();

		mesh->DestroyComponent();
	}

	_instancedStaticeMeshes.clear();

	_clearTransients();
}

void UMeshSimulation::tick( float deltaT )
{
	updateInstances();
}

void UMeshSimulation::tick_paused(float deltaT)
{
	updateInstances();
}

bool UMeshSimulation::canRunInEditor()
{
	return true;
}

void UMeshSimulation::snapshotToActor(AActor * actor)
{
	auto copyISMCs = [&](auto& pair)
	{
		FGraphMeshKey graphMeshKey = pair.first;
		UInstancedStaticMeshComponent * instance = pair.second;

		UInstancedStaticMeshComponent * newInstance = NewObject<UInstancedStaticMeshComponent>(actor);

		newInstance->SetStaticMesh(instance->GetStaticMesh());

		for (int32 mi = 0; mi < instance->GetNumMaterials(); mi++)
			newInstance->SetMaterial(mi, instance->GetMaterial(mi));

		newInstance->SetCollisionEnabled(ECollisionEnabled::NoCollision);
		newInstance->SetCollisionProfileName(TEXT("NoCollision"));
		newInstance->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
		newInstance->SetRenderCustomDepth(instance->bRenderCustomDepth);
		newInstance->SetVisibility(instance->bVisible);

		actor->AddInstanceComponent(newInstance);

		newInstance->RegisterComponent();

		// copy the instances
		int32 n = instance->GetInstanceCount();
		for (int32 i = 0; i < n; i++)
		{
			FTransform transform;

			instance->GetInstanceTransform(i, transform, true);

			newInstance->AddInstanceWorldSpace(transform);
		}
	};


	for (auto& pair : _instancedStaticeMeshes)
		copyISMCs(pair);

	for (auto& pair : _desaturatedInstancedStaticMeshes)
		copyISMCs(pair);
}

void UMeshSimulation::updateInstances()
{
	auto& meshes = graph->componentStorage<FGraphMesh>();

	USceneComponent * root = actor->GetRootComponent();
	FTransform rootTransform = root->GetComponentToWorld();

	// count number of visible then resize ISMCs
	std::unordered_map<UInstancedStaticMeshComponent*, int> visibleCount;
	TMap<UInstancedStaticMeshComponent*, int> desaturatedCount;

	for (FGraphMesh& mesh : meshes)
	{
		// if a mesh doesn't have a static mesh, or its invisible, we will skip it in the following loops
		// as well
		if (mesh.visible && mesh.staticMesh && mesh.isValid())
		{
			UInstancedStaticMeshComponent * ismc = mesh._transientInstanceMesh;

			// it's a new instance
			if (!ismc || ismc->GetStaticMesh() != mesh.staticMesh || ismc->GetMaterial(0) != mesh.material)
			{
				ismc = getOrCreateInstancedMesh(mesh);

				// assign the transient mesh
				mesh._transientInstanceMesh = ismc;
			}
		}
		else 
			continue;

		if (!mesh.desaturated)
		{
			UInstancedStaticMeshComponent * ismc = mesh._transientInstanceMesh;

			int& count = visibleCount[ismc];

			count++;
		}

		if (mesh.desaturated)
		{
			UInstancedStaticMeshComponent * ismc = _desaturatedInstancedStaticMeshes[mesh];

			if (ismc)
			{
				int& count = desaturatedCount.FindOrAdd(ismc);

				count++;
			}
		}



	}

	// make sure all ismcs are in the count table
	for (auto& pair : _instancedStaticeMeshes)
	{
		visibleCount[pair.second];
	}

	for (auto& pair : _desaturatedInstancedStaticMeshes)
	{
		desaturatedCount.FindOrAdd(pair.second);
	}

	// resize the visible ISMCs
	for (auto& pair : visibleCount)
	{
		UInstancedStaticMeshComponent * ismc = pair.first;

		int32 count = pair.second;

		// size down
		if (ismc->GetInstanceCount() > count)
		{
			do 
			{
				ismc->RemoveInstance(ismc->GetInstanceCount() - 1);
			} while (ismc->GetInstanceCount() > count);
		}
		// size up
		else if (ismc->GetInstanceCount() < count)
		{
			do
			{
				ismc->AddInstance(FTransform::Identity);
			} while (ismc->GetInstanceCount() < count);
		}

		_instanceToNodeHandle[ismc].resize(count);
	}

	// resize the desaturated ISMCs
	for (auto& pair : desaturatedCount)
	{
		UInstancedStaticMeshComponent * ismc = pair.Key;

		int32 count = pair.Value;

		// size down
		if (ismc->GetInstanceCount() > count)
		{
			do
			{
				ismc->RemoveInstance(ismc->GetInstanceCount() - 1);
			} while (ismc->GetInstanceCount() > count);
		}
		// size up
		else if (ismc->GetInstanceCount() < count)
		{
			do
			{
				ismc->AddInstance(FTransform::Identity);
			} while (ismc->GetInstanceCount() < count);
		}
	}

	// update positions
	TMap<UInstancedStaticMeshComponent*, int> visibleIndexTable;
	TMap<UInstancedStaticMeshComponent*, int> desaturatedIndexTable;

	for(FGraphMesh& mesh : meshes)
	{
		if( !mesh.visible || mesh.staticMesh == nullptr || !mesh.isValid())
			continue;

		FGraphNode& node = graph->node( mesh.nodeIndex );

		FTransform transform( node.orientation, node.position, FVector( node.scale ) );

		if (!mesh.desaturated)
		{
			UInstancedStaticMeshComponent * ismc = mesh._transientInstanceMesh;

			int& index = visibleIndexTable.FindOrAdd(ismc);

			ismc->UpdateInstanceTransform(index, transform, false, false, true);

			_instanceToNodeHandle[ismc][index] = mesh.nodeHandle();

			index++;
		}

		if (mesh.desaturated)
		{
			UInstancedStaticMeshComponent * ismc = _desaturatedInstancedStaticMeshes[mesh];

			int& index = desaturatedIndexTable.FindOrAdd(ismc);

			ismc->UpdateInstanceTransform(index, transform, false, false, true);

			index++;
		}
	}

	for(auto& pair : _instancedStaticeMeshes)
		pair.second->MarkRenderStateDirty();

	for (auto& pair : _desaturatedInstancedStaticMeshes)
		pair.second->MarkRenderStateDirty();
}

void UMeshSimulation::componentAdded( FGraphNodeHandle nodeHandle, ComponentType type )
{
	if (type != componentType<FGraphMesh>())
		return;

	FGraphMesh * mesh = graph->componentStorage<FGraphMesh>().componentPtrForNode(nodeHandle);

	mesh->_transientInstanceMesh = nullptr;
}

void UMeshSimulation::componentRemoved( FGraphNodeHandle nodeHandle, ComponentType type )
{

}

void UMeshSimulation::componentUpdated( FGraphNodeHandle nodeHandle, ComponentType type )
{

}

TArray<UInstancedStaticMeshComponent*> UMeshSimulation::instancedStaticMeshes()
{
	TArray<UInstancedStaticMeshComponent*> meshes;

	for(auto& pair : _instancedStaticeMeshes)
		meshes.Add( pair.second );

	return meshes;
}

FGraphNodeHandle UMeshSimulation::nodeForInstance(UInstancedStaticMeshComponent * mesh, int32 instance)
{
	auto& handles = _instanceToNodeHandle[mesh];

	if (handles.empty())
		return FGraphNodeHandle::null;

	return handles[instance];
}

FBox UMeshSimulation::boundsBoxForNode(FGraphNodeHandle handle)
{
	FBox box;

	if (!handle) return box;

	auto& storage = graph->componentStorage<FGraphMesh>();

	FGraphMesh * mesh = storage.componentPtrForNode(handle);

	if (!mesh) return box;

	FGraphNode& node = graph->node(handle);

	UInstancedStaticMeshComponent * ismc = mesh->_transientInstanceMesh;

	FTransform localTransform(node.orientation, node.position, FVector(node.scale));

	FTransform worldTransform = localTransform * ismc->GetComponentTransform();

	if (!ismc) return box;

	box = ismc->GetStaticMesh()->GetBoundingBox();

	return box.TransformBy(worldTransform);
}

UInstancedStaticMeshComponent * UMeshSimulation::_createISMCs(FGraphMesh& graphMesh)
{
	UInstancedStaticMeshComponent * instancedMesh = NewObject<UInstancedStaticMeshComponent>(actor, NAME_None, RF_Transient);

#if WITH_EDITOR
	if (GetWorld()->WorldType == EWorldType::Editor)
	{
		//instancedMesh->bIsEditorOnly = true;
		//instancedMesh->bHiddenInGame = true;
	}

#endif
	instancedMesh->SetStaticMesh(graphMesh.staticMesh);

	// we only have one material that we can set, but we need to set it for all the sections
	if (graphMesh.staticMesh)
	{
		int mati = 0;
		for (auto& mat : graphMesh.staticMesh->StaticMaterials)
			instancedMesh->SetMaterial(mati, graphMesh.material);
	}
	instancedMesh->SetMobility(EComponentMobility::Movable);
	instancedMesh->SetCollisionEnabled(ECollisionEnabled::NoCollision);
	instancedMesh->UseDynamicInstanceBuffer = true;
	instancedMesh->KeepInstanceBufferCPUAccess = true;
	instancedMesh->SetCollisionProfileName(TEXT("NoCollision"));
	instancedMesh->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
	instancedMesh->RegisterComponent();

	instancedMesh->SetRelativeScale3D(instancedMesh->GetRelativeTransform().GetScale3D());

	return instancedMesh;
}

UInstancedStaticMeshComponent * UMeshSimulation::getOrCreateInstancedMesh( FGraphMesh& graphMesh )
{
	UInstancedStaticMeshComponent * instancedMesh = nullptr;

	auto found = _instancedStaticeMeshes.find( graphMesh );

	if(found == _instancedStaticeMeshes.end())
	{
		instancedMesh = _createISMCs(graphMesh);

		UInstancedStaticMeshComponent * desaturatedMesh = _createISMCs(graphMesh);

		desaturatedMesh->SetRenderCustomDepth(true);

		_instancedStaticeMeshes[graphMesh] = instancedMesh;
		_desaturatedInstancedStaticMeshes[graphMesh] = desaturatedMesh;
	}
	else
		instancedMesh = found->second;

	return instancedMesh;
}

void UMeshSimulation::_clearTransients()
{
	// clear the transient mesh property, this could have been set to some other simulation
	auto& meshes = graph->componentStorage<FGraphMesh>();

	for (FGraphMesh& mesh : meshes)
	{
		mesh._transientInstanceMesh = nullptr;
	}
}