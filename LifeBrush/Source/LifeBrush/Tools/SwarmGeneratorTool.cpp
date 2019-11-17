//
//  Created by Timothy Davison on 2019-10-25.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementEditor/SwarmGenerator.h"

#include "SwarmGeneratorTool.h"

void USwarmGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	_generator = elementEditor->generator<USwarmGenerator>();
}

USwarmGenerator * USwarmGeneratorTool::swarmGenerator()
{
	return Cast<USwarmGenerator>(_generator);
}

void USwarmGeneratorTool::oneHandStart(UPrimitiveComponent * hand)
{
	auto nearby = _overlappingPrototypes(hand, brushMaxRadius * 1.5f);

	if (nearby.Num() == 0)
		_mode = Mode::Painting;
	else
		_mode = Mode::Selecting;

	if (_mode == Mode::Painting)
		Super::oneHandStart(hand);
}

void USwarmGeneratorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_mode == Mode::Painting)
		Super::oneHandEnd(hand);
	else if (_mode == Mode::Selecting)
	{
		if (_selection)
			swarmGenerator()->setPrototype(_selection->handleInRuleGraph);
		else
			swarmGenerator()->setPrototype(FGraphNodeHandle::null);
	}
}

void USwarmGeneratorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (_mode == Mode::Painting)
		Super::tickOneHand(dt, hand, lastTransform);
	else
	{
		_tickSelection(dt, hand, lastTransform);
	}
}

TArray<ASGPrototypeActor*> USwarmGeneratorTool::_overlappingPrototypes(UPrimitiveComponent * hand, float radius)
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

	TArray<ASGPrototypeActor*> filamentActors;


	for (auto overlap : overlaps)
	{
		ASGPrototypeActor * filamentActor = Cast<ASGPrototypeActor>(overlap.Actor.Get());

		if (!filamentActor)
			continue;

		filamentActors.Add(filamentActor);
	}

	return filamentActors;
}

void USwarmGeneratorTool::_tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	// find the nearest filament actor
	auto filaments = _overlappingPrototypes(hand, _brushRadius());

	ASGPrototypeActor * nearest = nullptr;
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
