// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/ObjectSimulation.h"

#include "PhysXSimulation.generated.h"

// Allows a node to be a member of a flex rigid object, but without a particle. The position of the node is updated
// relative to the position and orientation of the connected FFlexRigidBodyObject. Connect with a FFlexRigidBodyConnection.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FPXSphere : public FGraphObject
{
	GENERATED_BODY()
public:
	// These are transients, computed at load or component-add
	float radius;
	
	float mass;
};

UCLASS()
class LIFEBRUSH_API UPXSimulation : public UObjectSimulation
{
	GENERATED_BODY()

protected:
	virtual void attach();

public:
	virtual void tick(float deltaT) override;

	std::shared_ptr<tcodsMeshInterface> meshInterface;

	FTransform toWorld;
};