// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "VRTool.h"

#include <functional>

#include "SnapController.h"

#include "GrabTool.generated.h"

UINTERFACE(Blueprintable)
class UGrabDelegate : public UInterface
{
	GENERATED_BODY()
};

class IGrabDelegate
{
	GENERATED_BODY()

public:
	UFUNCTION(BlueprintCallable, BlueprintNativeEvent, Category = "VRUI")
	void didCancel(UGrabTool * grabTool, AActor * draggingActor);

	UFUNCTION(BlueprintCallable, BlueprintNativeEvent, Category = "VRUI")
	void didPlace(UGrabTool * grabTool, AActor * draggingActor);
};

UCLASS(Blueprintable)
class LIFEBRUSH_API UGrabTool : public UTool
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	FSnapController snapController;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float placePastDistance = 20.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "VRUI", meta = (MustImplement = "GrabDelegate"))
	UObject * delegate;

	// return true if it should be selected
	std::function<bool(AActor*)> filter_isSelectable;


public:
	virtual void oneHandStart(UPrimitiveComponent * hand);
	virtual void oneHandEnd(UPrimitiveComponent * hand);

	virtual void twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB);
	virtual void twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB);

	// Called when the tool is released. Any action should be aborted.
	virtual void loseFocus();

	// Called when the tool gains control.
	virtual void focused();

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform);

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform transformA,
		FTransform transformB
	);

	virtual void tick(float dt);

	virtual void faceDown_released();
	virtual void faceUp_released(USceneComponent * interactionPoint /* = nullptr */);

	void selectActor(AActor * object, UPrimitiveComponent * selectionPoint);

	void finishPlacement();
	void cancelPlacement();

protected:

	FTransform _transform();

	void _triggerEnd();

protected:
	UPROPERTY()
	AActor * _selectedActor = nullptr;

	FVector _objectPosition;
	FQuat _objectRotation = FQuat::Identity;

	FQuat _lastObjectRotation = FQuat::Identity;


	float _previewScale = 1.0f;
	float _animationTime = 0.0f;
	float _endScale = 1.0f;
	float _startScale = 1.0f;
	float _animationDuration = 0.18f;

	FTransform _startTransform;
	FTransform _endTransform;

	float _twoHandStartDistance = 1.0f;
	float _twoHandStartScale = 0.1f;

	bool _placingObject = false;

	bool _needUpdate = false;

	enum class RotationType
	{
		None,
		OneHand,
		TwoHand
	};

	RotationType _rotationType = RotationType::None;
};