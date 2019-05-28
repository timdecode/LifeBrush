// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "Visualization/Timeline.h"

#include "VRTool.h"
#include "Tools/RegionGrowingToolInit.h"


#include "EventVisualizationTool.generated.h"

class UFlexSimulationComponent;
class UCameraComponent;
struct FFlexSimulation;


UCLASS(BlueprintType)
class LIFEBRUSH_API UEventVisualizationTool : public UBrushTool
{
	GENERATED_BODY()

public:
	void init(FUToolInitProperties& initProperties,
		FFlexSimulation * flexSimulation_in,
		UCameraComponent * camera_in)
	{
		UTool::init(initProperties);

		_flexSimulation = flexSimulation_in;
		_camera = camera_in;
	}

	virtual void gainFocus() override;
	virtual void loseFocus() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tick(float dt);

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void faceDown_released() override;
	virtual void faceUp_released(USceneComponent * interactionPoint /* = nullptr */) override;

protected:
	void _selectEvent(float dt, UPrimitiveComponent * hand, FTransform lastTransform);

	void _traceSelection(TSet<USEGraphEvent*>& selection);

protected:
	FFlexSimulation * _flexSimulation;
	UCameraComponent * _camera;

	bool physicalInteraction = false;

	UPROPERTY()
	TSet<USEGraphEvent*> _selection;
};





UCLASS(BlueprintType)
class LIFEBRUSH_API UAgentPathlineTool : public UBrushTool
{
	GENERATED_BODY()

public:
	void init(FUToolInitProperties& initProperties,
		FFlexSimulation * flexSimulation_in,
		UCameraComponent * camera_in)
	{
		UTool::init(initProperties);

		_flexSimulation = flexSimulation_in;
		_camera = camera_in;
	}

	virtual void gainFocus() override;
	virtual void loseFocus() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tick(float dt);

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void faceDown_released() override;
	virtual void faceUp_released(USceneComponent * interactionPoint /* = nullptr */) override;

protected:
	void _selectAgent(float dt, UPrimitiveComponent * hand, FTransform lastTransform);

	void _pathlinesForSelection(TSet<FGraphNodeHandle>& selection);

protected:
	FFlexSimulation * _flexSimulation;
	class UCameraComponent * _camera;

	bool physicalInteraction = false;

	UPROPERTY()
	TSet<FGraphNodeHandle> _selection;
};

UENUM(BlueprintType)
enum class EPhysicalInteractionType : uint8
{
	Punch UMETA(DisplayName = "Punch"),
	Grab UMETA(DisplayName = "Grab"),
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UPhysicalInteractionTool : public UBrushTool
{
	GENERATED_BODY()

public:
	void init(FUToolInitProperties& initProperties,
		FFlexSimulation * flexSimulation_in,
		UCameraComponent * camera_in)
	{
		UTool::init(initProperties);

		_flexSimulation = flexSimulation_in;
		_camera = camera_in;
	}

	virtual void loseFocus() override;

	virtual void faceDown_released() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual bool shouldShowBrush() override;

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	EPhysicalInteractionType interactionMode = EPhysicalInteractionType::Grab;

protected:
	FFlexSimulation * _flexSimulation;
	UCameraComponent * _camera;

	FGraphNodeHandle _grabbedNode;

	FMatrix _grabbedIdentityMatrix;
	FTransform _grabbedTransform;
	FTransform _startHand;

	TMap<FGraphNodeHandle, FVector> _cachedCalculatedVelocity;
};