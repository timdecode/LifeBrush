// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"

#include "ParticleSlots.generated.h"	

USTRUCT()
struct LIFEBRUSH_API FParticleSlot : public FGraphObject
{
	GENERATED_BODY()

protected:
	UPROPERTY() FGraphNodeHandle _pigeon = FGraphNodeHandle::null;

public:
	UPROPERTY() bool seek = false;
	
	bool occupied() { return _pigeon != FGraphNodeHandle::null; }

	// Sets the pigeon if it doesn't have an FBoundObject attached
	void bind(FGraphNodeHandle target, FGraph& graph);

	// Detached the pigeon and removes its FBoundObject
	FGraphNodeHandle unbind(FGraph& graph);

	FGraphNodeHandle pigeon() { return _pigeon; }

public:
	// TSlotType is a FDNASimulationSlot and TWork is a lambda of type bool(TSlotType& slot, FGraphNodeHandle overlapHandle).
	// The queryFunction returns true if the element is valid for the slot.
	template<class TSlotType, typename TWork>
	static void tickSlots(float dt, TWork queryFunction, UObjectSimulation * simulation, float searchRadius)
	{
		typedef UMLElementSimulation::BVH_t BVH_T;

		auto doNothing = [](TSlotType& slot, FGraphNodeHandle) {};

		tickSlots_withCallback<TSlotType>(dt, queryFunction, simulation, searchRadius, doNothing);
	}

	template<class TSlotType, typename TWork, typename TBindCallback>
	static void tickSlots_withCallback(float dt, TWork queryFunction, UObjectSimulation * simulation, float searchRadius, TBindCallback bindCallback)
	{
		typedef UMLElementSimulation::BVH_t BVH_T;

		BVH_T& bvh = simulation->simulationManager->simulation<UMLElementSimulation>()->elementBVH;

		FGraph& graph = *simulation->graph;

		auto& elements = graph.componentStorage<FMLElement>();
		auto& slots = graph.componentStorage<TSlotType>();
		auto& velocities = graph.componentStorage<FVelocityGraphObject>();
		auto& boundObjects = graph.componentStorage<FBoundToParticleSlot>();
		auto& rigids = graph.componentStorage<FFlexRigidBodyObject>();

		for (TSlotType& slot : slots)
		{
			// attract the pigeon
			if (slot.occupied())
			{
				FGraphNodeHandle pigeon = slot.pigeon();
				FVector slotPosition = graph.node(slot.nodeHandle()).position;
				FQuat slotOrientation = graph.node(slot.nodeHandle()).orientation;

				// teleport
				FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(graph, pigeon);

				if (rigidHandle)
				{
					rigids.componentPtrForNode(rigidHandle)->setRotationPosition(slotOrientation, slotPosition, graph);
				}
				else
				{
					graph.node(pigeon).position = slotPosition;
					graph.node(pigeon).orientation = slotOrientation;
				}
			}
			// find a pigeon
			else if (!slot.occupied() && slot.seek)
			{
				FGraphNode node = graph.node(slot.nodeHandle());

				FGraphNodeHandle nearestHandle = FGraphNodeHandle::null;
				float minDistance = std::numeric_limits<float>::max();

				unrealAABB::AABB query(node.position, searchRadius);

				bvh.query(query, [&](unsigned int elementIndex) {
					FGraphNodeHandle handle(elementIndex);

					const FVector position = graph.node(handle).position;

					// we can't bind to something that someone else owns
					if (boundObjects.componentPtrForNode(handle)) return true;

					bool validElement = queryFunction(slot, handle);

					// this element doesn't match our query, keep searching
					if (!validElement) return true;

					float distance = FVector::Dist(position, node.position);

					if (distance < minDistance)
					{
						minDistance = distance;
						nearestHandle = handle;
					}

					return true; // keep querying
				});

				slot.bind(nearestHandle, graph);

				if( nearestHandle )
					bindCallback(slot, nearestHandle);
			}
		}
	}
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FBoundToParticleSlot : public FGraphObject
{
	GENERATED_BODY()

public:

};