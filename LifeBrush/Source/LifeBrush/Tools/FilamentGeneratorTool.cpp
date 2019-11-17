//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementEditor/FilamentGenerator.h"

#include "FilamentGeneratorTool.h"

void UFilamentGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	_generator = elementEditor->generator<UFilamentGenerator>();
}

void UFilamentGeneratorTool::oneHandStart(UPrimitiveComponent * hand)
{
	auto nearby = overlappingPrototypes(hand, brushMaxRadius * 1.5f);

	if (nearby.Num() == 0)
		_mode = Mode::Painting;
	else
		_mode = Mode::Selecting;

	if( _mode == Mode::Painting )
		Super::oneHandStart(hand);
}

void UFilamentGeneratorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_mode == Mode::Painting)
		Super::oneHandEnd(hand);
	else if (_mode == Mode::Selecting)
	{
		if (_selection)
			filamentGenerator()->setFilamentPrototype(_selection->filamentPrototype);
		else
			filamentGenerator()->clearFilamentPrototype();
	}
}

void UFilamentGeneratorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (_mode == Mode::Painting)
		Super::tickOneHand(dt, hand, lastTransform);
	else
	{
		_tickSelection(dt, hand, lastTransform);
	}
}

void UFilamentGeneratorTool::faceDown_released()
{
	filamentGenerator()->togglePlay();
}

void UFilamentGeneratorTool::_tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	// find the nearest filament actor
	auto filaments = overlappingPrototypes(hand, _brushRadius());

	AFilamentPrototypeActor * nearest = nullptr;
	float nearestDistance = std::numeric_limits<float>::max();

	FVector handPosition = hand->GetComponentLocation();

	for (auto actor : filaments)
	{
		float distance = FVector::Dist(handPosition, actor->GetActorLocation());

		if (distance < nearestDistance)
		{
			nearestDistance = distance;
			nearest = actor;
		}
	}

	if (_selection)
	{
		UStaticMeshComponent * mesh = _selection->FindComponentByClass<UStaticMeshComponent>();

		// we visualize selection through a post-process effect on the custom depth
		if (mesh)
		{
			mesh->SetRenderCustomDepth(false);
		}
	}

	if (nearest)
	{
		_selection = nearest;

		UStaticMeshComponent * mesh = _selection->FindComponentByClass<UStaticMeshComponent>();

		// we visualize selection through a post-process effect on the custom depth
		if (mesh)
		{
			mesh->SetRenderCustomDepth(true);
		}
	}
	else
		_selection = nullptr;
}

TArray<AFilamentPrototypeActor*> UFilamentGeneratorTool::overlappingPrototypes(UPrimitiveComponent * hand, float radius)
{
	// if we are near an AFilamentPrototypeActor we are in selection mode
	TArray<FOverlapResult> overlaps;

	FCollisionShape collisionShape = FCollisionShape::MakeSphere(radius);

	FCollisionQueryParams params;

	params.AddIgnoredActor(elementEditor->GetOwner());
	params.AddIgnoredActor(hand->GetOwner());

	bool overlap = elementEditor->GetWorld()->OverlapMultiByChannel(
		overlaps,
		hand->GetComponentLocation(),
		hand->GetComponentRotation().Quaternion(),
		ECollisionChannel::ECC_WorldDynamic,
		collisionShape,
		params
	);

	TArray<AFilamentPrototypeActor*> filamentActors;


	for (auto overlap : overlaps)
	{
		AFilamentPrototypeActor * filamentActor = Cast<AFilamentPrototypeActor>(overlap.Actor.Get());

		if (!filamentActor)
			continue;

		filamentActors.Add(filamentActor);
	}

	return filamentActors;
}

UFilamentGenerator * UFilamentGeneratorTool::filamentGenerator()
{
	return Cast<UFilamentGenerator>(_generator);
}


