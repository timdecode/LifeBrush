//
//  Created by Timothy Davison on 2019-10-14
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementEditor/MolecularLegoGenerator.h"

#include "MolecularLegoGeneratorTool.h"

void UMolecularLegoGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	molecularGenerator = elementEditor->generator<UMolecularLegoGenerator>();

	if (!molecularGenerator)
		molecularGenerator = NewObject<UMolecularLegoGenerator>(this, TEXT("MolecularLegoGenerator"));

	_generator = molecularGenerator;
}

void UMolecularLegoGeneratorTool::faceUp_released(USceneComponent * interactionPoint /*= nullptr*/)
{
	if (molecularGenerator->mode == EMolecularLegotGeneratorMode::AlongPath)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::Grow;
	else if (molecularGenerator->mode == EMolecularLegotGeneratorMode::Grow)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::PlaceOne;
	else if (molecularGenerator->mode == EMolecularLegotGeneratorMode::PlaceOne)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::AlongPath;
}

void UMolecularLegoGeneratorTool::oneHandStart(UPrimitiveComponent * hand)
{
	auto nearby = _overlappingPrototypes(hand, brushMaxRadius * 1.5f);

	if (nearby.Num() == 0)
		_mode = Mode::Painting;
	else
		_mode = Mode::Selecting;

	if (_mode == Mode::Painting)
		Super::oneHandStart(hand);
}

void UMolecularLegoGeneratorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_mode == Mode::Painting)
		Super::oneHandEnd(hand);
	else if (_mode == Mode::Selecting)
	{
		TArray<AElementActor*> actors;

		if( _selection ) actors.Add(_selection);

		if (_selection)
			molecularGenerator->setSelection(actors);
		else
			molecularGenerator->setSelection(actors);
	}
}

void UMolecularLegoGeneratorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (_mode == Mode::Painting)
		Super::tickOneHand(dt, hand, lastTransform);
	else
	{
		_tickSelection(dt, hand, lastTransform);
	}
}

void UMolecularLegoGeneratorTool::loseFocus()
{
	Super::loseFocus();

	_hideSelection();

	_selection = nullptr;
}

TArray<AElementActor*> UMolecularLegoGeneratorTool::_overlappingPrototypes(UPrimitiveComponent * hand, float radius)
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

	TArray<AElementActor*> elementActors;


	for (auto overlap : overlaps)
	{
		AElementActor * elementActor = Cast<AElementActor>(overlap.Actor.Get());

		if (!elementActor)
			continue;

		elementActors.Add(elementActor);
	}

	return elementActors;
}

void UMolecularLegoGeneratorTool::_tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	// find the nearest element actor
	auto filaments = _overlappingPrototypes(hand, _brushRadius());

	AElementActor * nearest = nullptr;
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

	_hideSelection();

	if (nearest)
	{
		_selection = nearest;

		_showSelection();
	}
	else
		_selection = nullptr;
}

void UMolecularLegoGeneratorTool::_hideSelection()
{
	if (_selection)
		_selection->hideSelectionOutline();
}

void UMolecularLegoGeneratorTool::_showSelection()
{
	if (_selection)
		_selection->showSelectionOutline();
}



