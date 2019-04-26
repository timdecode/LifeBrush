// Copyright 2018 Timothy Davison, all rights reserved.

#include "LifeBrush.h"

#include "GrabTool.h"

void UGrabTool::oneHandStart(UPrimitiveComponent * hand)
{
	_rotationType = RotationType::OneHand;

	if (hand->IsOverlappingActor(_selectedActor) || !filter_isSelectable)
		return;


	FVector closestPoint;
	float distance = _selectedActor->ActorGetDistanceToCollision(hand->GetComponentLocation(), ECollisionChannel::ECC_WorldDynamic, closestPoint);

	// too far, place it
	if (distance >= placePastDistance)
	{
		finishPlacement();
		return;
	}

	// we need to find a new overlapping actor
	//TSet<AActor*> overlaps;

	//hand->GetOverlappingActors(overlaps);

	//for (AActor * actor : overlaps)
	//{
	//	if (filter_isSelectable(actor))
	//		selectActor(actor, hand);
	//}
}

void UGrabTool::oneHandEnd(UPrimitiveComponent * hand)
{
	_triggerEnd();
}

void UGrabTool::twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
	_twoHandStartDistance = FVector::Dist(handA->GetComponentLocation(), handB->GetComponentLocation());
	_twoHandStartScale = _previewScale;// / _endScale;

	_rotationType = RotationType::TwoHand;
}

void UGrabTool::twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
}

void UGrabTool::loseFocus()
{
	_rotationType = RotationType::None;
	_placingObject = false;
}

void UGrabTool::gainFocus()
{

}

void UGrabTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform)
{
	_needUpdate = true;

	FVector position = _selectionA->GetComponentLocation();
	FQuat rotation = _selectionA->GetComponentQuat() * lastToWorldTransform.GetRotation().Inverse(); // same idea as (point - _lastPoint) to find delta position

	_objectPosition = rotation.RotateVector(_objectPosition - lastToWorldTransform.GetLocation()) + position;
	_objectRotation = rotation * _objectRotation;

	FVector endPosition = _endTransform.GetLocation();
	_endTransform.SetLocation(rotation.RotateVector(endPosition - lastToWorldTransform.GetLocation()) + position);
	_endTransform.SetRotation(rotation * _endTransform.GetRotation());
}

void UGrabTool::tickTwoHand(float dt, UPrimitiveComponent * handA, UPrimitiveComponent * handB, FTransform lastTransformA, FTransform lastTransformB)
{
	_needUpdate = true;

	FVector aPosition = handA->GetComponentLocation();
	FVector bPosition = handB->GetComponentLocation();

	FVector dirLast = lastTransformB.GetLocation() - lastTransformA.GetLocation();
	FVector dirCur = bPosition - aPosition;

	FQuat rotation = FQuat::FindBetweenVectors(dirLast, dirCur);

	float fracD = dirCur.Size() / dirLast.Size();

	_previewScale = (dirCur.Size() / _twoHandStartDistance) * _twoHandStartScale;

	_objectRotation = rotation * _objectRotation;

	_objectPosition = rotation.RotateVector(fracD * (_objectPosition - lastTransformA.GetLocation())) + aPosition;
}

void UGrabTool::tick(float dt)
{
	if (!_selectedActor)
		return;

	if (_animationTime < _animationDuration)
	{
		float alpha = FMath::InterpEaseInOut(0.0f, 1.0f, _animationTime / _animationDuration, 1.2f);

		_animationTime += dt;

		_needUpdate = false;

		FTransform transform = _startTransform;

		transform.LerpTranslationScale3D(_startTransform, _endTransform, ScalarRegister(alpha));
		transform.SetRotation(_endTransform.GetRotation()); // no slerp needed

		_objectRotation = transform.GetRotation();
		_objectPosition = transform.GetLocation();
		_previewScale = transform.GetScale3D().GetMax();

		_selectedActor->SetActorTransform(transform, false, nullptr, ETeleportType::TeleportPhysics);
	}

	if (_needUpdate)
	{
		FTransform transform = _transform();

		_selectedActor->SetActorTransform(transform, false, nullptr, ETeleportType::TeleportPhysics);

		_needUpdate = false;
	}
}

void UGrabTool::faceUp_released(USceneComponent * interactionPoint /* = nullptr */)
{
	cancelPlacement();
}

void UGrabTool::faceDown_released()
{
	finishPlacement();
}

void UGrabTool::selectActor( AActor * object, UPrimitiveComponent * selectionPoint )
{
	_updateLast();

	_placingObject = false;

	_selectedActor = object;

	if (!_selectedActor)
		return;

	_placingObject = true;

	FTransform selectedTransform = _selectedActor->GetRootComponent()->GetComponentTransform();

	_rotationType = RotationType::OneHand;

	_objectRotation = selectedTransform.GetRotation();
	_objectPosition = selectedTransform.GetLocation();

	_startTransform = selectedTransform;

	_endScale = selectedTransform.GetMaximumAxisScale();


	{
		_endTransform = _startTransform;
		_endTransform.SetScale3D(FVector(_endScale));

		FVector selectionLocation = selectionPoint->GetComponentLocation();
		FVector selectionRelativeLocation = _startTransform.InverseTransformPosition(selectionLocation);

		FVector selectionLocation_next = _endTransform.TransformPosition(selectionRelativeLocation);

		_endTransform.AddToTranslation(selectionLocation - selectionLocation_next);
	}

	_previewScale = selectedTransform.GetMaximumAxisScale();
	_startScale = _previewScale;

	_animationTime = 0.0f;

	_twoHandStartScale = 1.0f;
}

void UGrabTool::finishPlacement()
{
	_rotationType = RotationType::None;
	_placingObject = false;

	AActor * selected = _selectedActor;

	_selectedActor = nullptr;


	if (delegate && delegate->Implements<UGrabDelegate>())
	{
		IGrabDelegate::Execute_didPlace(delegate, this, selected);
	}
}

void UGrabTool::cancelPlacement()
{
	_rotationType = RotationType::None;
	_placingObject = false;

	AActor * selected = _selectedActor;

	_selectedActor = nullptr;

	if (delegate && delegate->Implements<UGrabDelegate>())
	{
		IGrabDelegate::Execute_didCancel(delegate, this, selected);
	}
}

FTransform UGrabTool::_transform()
{
	FTransform transform(
		_objectRotation,
		_objectPosition,
		FVector(_previewScale)
	);

	return transform;
}

void UGrabTool::_triggerEnd()
{
	if (!_selectedActor)
		return;

	_objectPosition = snapController.snappedPosition(_objectPosition);
	_objectRotation = snapController.snappedRotation(_objectRotation);

	_needUpdate = true;
}