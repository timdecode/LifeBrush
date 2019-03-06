// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include <cmath>

#include "CollectionSpace.h"


UInteractionSpace::UInteractionSpace()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	// ...
}

UCollectionSpace::UCollectionSpace()
{
	PrimaryComponentTick.bCanEverTick = true;

	bAutoActivate = true;
	bTickInEditor = false;
}


void UCollectionSpace::BeginPlay()
{
	Super::BeginPlay();

	_boundsMesh = NewObject<UStaticMeshComponent>(this);

	_boundsMesh->AttachToComponent(this, FAttachmentTransformRules::KeepRelativeTransform);

	_boundsMesh->SetStaticMesh(backgroundMesh);
	_boundsMesh->SetMaterial(0, backgroundMaterial);
	_boundsMesh->SetCollisionEnabled(ECollisionEnabled::NoCollision);

	_boundsMesh->RegisterComponent();
}


// Called every frame
void UCollectionSpace::TickComponent(float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	_velocity -= _velocity * damping * DeltaTime; 

	// apply velocity
	this->AddRelativeLocation(_velocity * DeltaTime, false, nullptr, ETeleportType::TeleportPhysics);

	FVector location = GetRelativeTransform().GetLocation();

	FBox boundsInRoot = _bounds.TransformBy(GetRelativeTransform());

	const float halfWidth = width * 0.5f;

	if( boundsInRoot.Min[Forward] > -halfWidth )
	{
		_velocity = FVector::ZeroVector;

		location[Forward] += (-halfWidth - boundsInRoot.Min[Forward]);

		this->SetRelativeLocation(location, false, nullptr, ETeleportType::TeleportPhysics);
	}
	else if( boundsInRoot.Max[Forward] < halfWidth )
	{
		_velocity = FVector::ZeroVector;

		location[Forward] += (halfWidth - boundsInRoot.Max[Forward]);

		this->SetRelativeLocation(location, false, nullptr, ETeleportType::TeleportPhysics);
	}
}

void UCollectionSpace::insertItemsAt(TArray<uint32> indices)
{

}

void UCollectionSpace::reloadData()
{
	_cells.Empty();

	if(!dataSource->Implements<UCollectionSpaceDataSource>())
		return;

	int32 n = ICollectionSpaceDataSource::Execute_count(dataSource);

	for (int32 i = 0; i < n; ++i)
	{
		UPrimitiveComponent * cell = ICollectionSpaceDataSource::Execute_primitiveCellAt(dataSource, i);

		_cells.Add(cell);
	}

	_layout();
}

void UCollectionSpace::begin_oneHand(UPrimitiveComponent * interactionPoint)
{
	FVector localPoint = this->GetComponentTransform().InverseTransformPosition(interactionPoint->GetComponentLocation());

	_interactionMode = InteractionMode::Pan;

	// check if we overlap one of the items
	for (int i = 0; i < _localBounds.Num(); ++i)
	{
		FBox& bounds = _localBounds[i];

		if (bounds.IsInside(localPoint) && delegate)
		{
			 ICollectionSpaceDelegate::Execute_didGrab(delegate, i, interactionPoint->GetComponentTransform(), _cells[i], _cells[i]->GetComponentTransform(), bounds);

			_interactionMode = InteractionMode::Grab;

			return;
		}
	}
}

void UCollectionSpace::update_oneHand(float dt, UPrimitiveComponent * interactionPoint, FTransform lastTransform)
{
	if (_interactionMode == InteractionMode::Pan)
	{
		_update_pan(dt, interactionPoint, lastTransform);
	}
}

void UCollectionSpace::_update_pan(float dt, UPrimitiveComponent * interactionPoint, FTransform lastTransform)
{
	FVector localPoint = this->GetComponentTransform().InverseTransformPosition(interactionPoint->GetComponentLocation());

	if (_bounds.IsInside(localPoint))
	{
		FTransform parentTransform = this->GetAttachParent()->GetComponentTransform();

		FVector pointInParent = parentTransform.InverseTransformPosition(interactionPoint->GetComponentLocation());
		FVector lastInParent = parentTransform.InverseTransformPosition(lastTransform.GetLocation());

		FVector delta = pointInParent - lastInParent;

		// only move in forward direction
		delta.X = 0.0f;
		delta.Z = 0.0f;

		_velocity = delta / dt;
	}
}

void UCollectionSpace::end_oneHand(UPrimitiveComponent * interactionPoint)
{

}

void UCollectionSpace::grab(UPrimitiveComponent * interactionPoint)
{

}

void UCollectionSpace::query(UPrimitiveComponent * interactionPoint)
{

}


TSet<int32> UCollectionSpace::selection()
{
	return _selection;
}

void UCollectionSpace::_layout()
{
	_layoutCells();
	_layoutBackground();
}

void UCollectionSpace::_layoutCells()
{
	const float totalWidth = _totalWidth();
	const float stepSize = _stepSize();

	FVector _startP = FVector::ZeroVector;
	_startP[Forward] = -totalWidth / 2.0f;
	_startP.Z = -(cellExtents * rows + cellSpacing * (rows - 1)) * .5f;

	const FVector startP = _startP;

	_bounds = FBox(EForceInit::ForceInitToZero);

	int n = _cells.Num();

	_localBounds.SetNumZeroed(_cells.Num(), true);

	FTransform toComponent = this->GetComponentTransform().Inverse();

	const int32 nColumns = std::ceil(float(n) / float(rows));

	for (int32 i = 0; i < n; ++i)
	{
		UPrimitiveComponent * cell = _cells[i];

		// cell->SetCollisionEnabled(ECollisionEnabled::QueryOnly);

		checkf(cell->GetAttachParent() == this, TEXT("The AttachParent of the cell must be the UCollectionSpace"));


		const float scale = _scaleForCell(*cell);
		const FVector offset = _offsetForCell(toComponent, *cell);

		FVector p = startP;;
		{
			const int32 rowIndex = std::floor(i / nColumns);
			const int32 colIndex = i % nColumns;

			p[Forward] += stepSize * colIndex;
			p.Z += stepSize * rowIndex;
		}

		FTransform transform(FQuat::Identity, p + offset, FVector(scale));

		cell->SetRelativeTransform(transform, false, nullptr, ETeleportType::TeleportPhysics);




		//p[Forward] += stepSize;

		FBoxSphereBounds localBounds = cell->Bounds.TransformBy(toComponent);

		_localBounds[i] = localBounds.GetBox();

		_bounds += localBounds.GetBox();
	}

	_bounds.Min.Z -= cellSpacing;
	_bounds.Max.Z += cellSpacing;

	_bounds.Min.Z -= cellSpacing;
	_bounds.Max.Z += cellSpacing;
}

void UCollectionSpace::_layoutBackground()
{
	FBoxSphereBounds bounds = _boundsMesh->GetStaticMesh()->GetBounds();

	const float yScale = _bounds.GetSize()[Forward] / bounds.GetBox().GetSize()[Forward];
	const float zScale = _bounds.GetSize().Z / bounds.GetBox().GetSize().Z;

	FVector location = _bounds.GetCenter();
	location.X = _bounds.Min.X;

	_boundsMesh->SetRelativeLocation(location);

	FVector scaleVector(0.02f, yScale, zScale);

	_boundsMesh->SetRelativeScale3D(scaleVector);
}

float UCollectionSpace::_totalWidth()
{
	return (float(_cells.Num()) * cellExtents) - (float(_cells.Num() - 1) * cellSpacing);

}

float UCollectionSpace::_stepSize()
{
	return cellExtents + cellSpacing;
}

float UCollectionSpace::_scaleForCell(UPrimitiveComponent& component)
{
	return (cellExtents * 0.5f) / component.Bounds.SphereRadius;
}

FVector UCollectionSpace::_offsetForCell(const FTransform& toComponent, UPrimitiveComponent& component)
{
	FBoxSphereBounds bounds = component.Bounds;

	// get it into the coordinate frame of this
	bounds = bounds.TransformBy(toComponent);

	FVector offset = bounds.Origin;
	offset[Forward] = -bounds.GetBox().Max[Forward];

	return -offset * _scaleForCell(component);
}
