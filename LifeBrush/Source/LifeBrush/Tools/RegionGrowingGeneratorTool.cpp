// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "RegionGrowingGeneratorTool.h"

// -------------------------------------------------------------------------------
// URegionGrowingElementGenerator_BrushTool
// -------------------------------------------------------------------------------

UGeneratorTool::~UGeneratorTool()
{
}

void UGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	UTool::init(initProperties);


}

void UGeneratorTool::gainFocus()
{
	Super::gainFocus();

	if(elementEditor)
		elementEditor->setCurrentGenerator(_generator);
}

void UGeneratorTool::loseFocus()
{
	Super::loseFocus();

	if(elementEditor)
		elementEditor->setCurrentGenerator(nullptr);
}

void UGeneratorTool::oneHandStart(UPrimitiveComponent * hand)
{
	Super::oneHandStart(hand);

	if (!_generator || !elementEditor)
		return;
	
	auto meshInterface = elementEditor->context()->meshInterface.get();

	FVector location = hand->GetComponentLocation();

	float radius = _brushRadius();

	FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface;

	if (_drawMode == EGenerativeDrawMode::Surface)
	{
		auto nearest = meshInterface->nearestPointOnMesh(location);

		surfaceIndex = nearest.surfaceIndex;
		location = nearest.point;
	}

	if (_tickMode == EGenerativeTickMode::Generating)
		_generator->beginBrushPath(location, radius, surfaceIndex);
	else if (_tickMode == EGenerativeTickMode::Erasing)
		_generator->beginEraserPath(location, radius, surfaceIndex);

	// choose the draw-mode
	if (automaticDrawMode == EAutomaticDrawMode::ByNearestSurface)
	{
		FVector location = hand->GetComponentLocation();

		auto meshInterface = elementEditor->context()->meshInterface.get();

		auto nearest = meshInterface->nearestPointOnMesh(location);

		if (FVector::DistSquared(nearest.point, location) < minSurfaceDistanceThreshold)
		{
			setDrawMode(EGenerativeDrawMode::Surface);
		}
		else
		{
			setDrawMode(EGenerativeDrawMode::Volumetric);
		}
	}
}

void UGeneratorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	Super::oneHandEnd(hand);

	if (!_generator)
		return;

	_generator->endBrushPath();
	_generator->endEraserPath();
}

void UGeneratorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	Super::tickOneHand(dt, hand, lastTransform);

	if (!_generator || !elementEditor)
		return;
	
	if (_tickMode == EGenerativeTickMode::Generating)
		_tickOneHand_generate(dt, hand, lastTransform);
	else
		_tickOneHand_erase(dt, hand, lastTransform);
}

void UGeneratorTool::_tickOneHand_generate(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	auto meshInterface = elementEditor->context()->meshInterface.get();

	FVector location = hand->GetComponentLocation();

	float radius = _brushRadius();

	if (_drawMode == EGenerativeDrawMode::Volumetric)
	{
		_generator->addBrushPoint(location, radius);
	}
	else if (_drawMode == EGenerativeDrawMode::Surface)
	{
		auto nearest = meshInterface->nearestPointOnMesh(location);

		_generator->addBrushPoint(nearest.point, radius, nearest.surfaceIndex);
	}
}

void UGeneratorTool::_tickOneHand_erase(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	auto meshInterface = elementEditor->context()->meshInterface.get();

	FVector location = hand->GetComponentLocation();

	float radius = _brushRadius();

	if (_drawMode == EGenerativeDrawMode::Volumetric)
	{
		_generator->eraseInRadiusAt(location, radius);
	}
	else if (_drawMode == EGenerativeDrawMode::Surface)
	{
		auto nearest = meshInterface->nearestPointOnMesh(location);

		_generator->eraseInRadiusAt(nearest.point, radius, nearest.surfaceIndex);
	}
}


void URegionGrowingGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	URegionGrowingGenerator * generator = elementEditor->generator<URegionGrowingGenerator>();

	if (!generator)
		generator = NewObject<URegionGrowingGenerator>(this, TEXT("generator"));

	generator->exemplarActor = exemplarActor;

	_generator = generator;
}

void URegionGrowingGeneratorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (!generator()) return;

	if (_isSelecting)
	{
		_updateBrushMesh();
		_tickSelection(dt, hand, lastTransform);
	}
	else
	{
		generator()->brushSize = _brushRadius();

		Super::tickOneHand(dt, hand, lastTransform);
	}
}

void URegionGrowingGeneratorTool::oneHandStart(UPrimitiveComponent * hand)
{
	const float radius = 20.0f;

	auto overlaps = _overlappingElementActors(hand, radius);

	if (overlaps.Num() > 0)
	{
		_isSelecting = true;

		_clearSelection();
	}
	else
	{
		_isSelecting = false;

		Super::oneHandStart(hand);
	}
}

void URegionGrowingGeneratorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_isSelecting)
	{
		_isSelecting = false;

		URegionGrowingGenerator * rgGenerator = generator();

		auto asVector = _toVector(_selection);

		rgGenerator->setExampleSelection(asVector);
	}
	else
	{
		Super::oneHandEnd(hand);
	}
}

void URegionGrowingGeneratorTool::gainFocus()
{
	if (auto rgGenerator = generator())
	{
		std::vector<AElementActor*> nothing;
		rgGenerator->setExampleSelection(nothing);
	}

	_clearSelection();

	_isSelecting = false;

	Super::gainFocus();
}

void URegionGrowingGeneratorTool::loseFocus()
{
	if (auto rgGenerator = generator())
	{
		std::vector<AElementActor*> nothing;
		rgGenerator->setExampleSelection(nothing);
	}

	_clearSelection();

	_isSelecting = false;

	Super::loseFocus();
}

void URegionGrowingGeneratorTool::_tickOneHand_generate(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (_drawMode == EGenerativeDrawMode::Volumetric)
	{
		generator()->setGenerationMode(EGenerationMode::SpacePainting);
	}
	else if (_drawMode == EGenerativeDrawMode::Surface)
	{
		generator()->setGenerationMode(EGenerationMode::SurfacePainting);
	}

	Super::_tickOneHand_generate(dt, hand, lastTransform);
}

void URegionGrowingGeneratorTool::_tickOneHand_erase(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	Super::_tickOneHand_erase(dt, hand, lastTransform);
}

void URegionGrowingGeneratorTool::_tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	// find the nearest element actor
	auto elementActors = _overlappingElementActors(hand, _brushRadius());

	for (auto actor : elementActors)
	{
		_selection.Add(actor);
	}

	_showSelection();
}

void URegionGrowingGeneratorTool::_clearSelection()
{
	for (AElementActor * elementActor : _selection)
	{
		elementActor->hideSelectionOutline();
	}

	_selection.Empty();
}

void URegionGrowingGeneratorTool::_showSelection()
{
	for (AElementActor * elementActor : _selection)
	{
		elementActor->showSelectionOutline();
	}
}

std::vector<AElementActor*> URegionGrowingGeneratorTool::_toVector(TSet<AElementActor*> aSet)
{
	std::vector<AElementActor*> result;

	for (AElementActor * actor : aSet)
		result.push_back(actor);

	return result;
}

TArray<AElementActor*> URegionGrowingGeneratorTool::_overlappingElementActors(UPrimitiveComponent * hand, float radius)
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
