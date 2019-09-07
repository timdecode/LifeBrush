// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "RegionGrowingComponent.h"
#include "WidgetComponent.h"
#include "ElementEditor/DiscreteElementEditorComponent.h"
#include "RegionGrowingGeneratorTool.h"
#include "SelectionTool.h"


USelectionTool::~USelectionTool()
{
	_destroySelectionMeshComponent();
}

void USelectionTool::postInit()
{

}

void USelectionTool::oneHandStart( UPrimitiveComponent * hand )
{
	_createSelectionMeshComponent( hand );
}

void USelectionTool::oneHandEnd( UPrimitiveComponent * hand )
{
	_destroySelectionMeshComponent();
}

float USelectionTool::_brushRadius()
{
	return 10.0f * selectionATriggerValue();
}

void USelectionTool::tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform )
{
	if(_selectionMeshComponent)
	{
		_selectionMeshComponent->SetRelativeScale3D( FVector( _brushRadius() * selectionMeshScaleFactor ) );
	}

	if (useGroupSelection)
		_tickGroupSelection(dt, hand);
	else
		_tickSelection(dt, hand);
}

void USelectionTool::_tickGroupSelection(float dt, UPrimitiveComponent * hand)
{
	// perform an overlap test at the hand position to find actors that overlap
	TArray<FOverlapResult> overlaps;

	FCollisionShape collisionShape = FCollisionShape::MakeSphere(_brushRadius());

	FCollisionQueryParams params;

	UActorComponent * root;
	if (discreteEditorComponent)
		root = discreteEditorComponent;
	else
		root = regionGrowingComponent;


	params.AddIgnoredActor(root->GetOwner());
	params.AddIgnoredActor(hand->GetOwner());

	bool overlap = root->GetWorld()->OverlapMultiByChannel(
		overlaps,
		hand->GetComponentLocation(),
		hand->GetComponentRotation().Quaternion(),
		ECollisionChannel::ECC_WorldDynamic,
		collisionShape,
		params
	);

	// bail if there is nothing
	if (overlaps.Num() == 0)
		return;

	// otherwise, set the selection based on the groupName of the first overlapping object
	TSet<AElementActor*> newSelection;

	AElementActor * elementActor = nullptr;

	for (FOverlapResult overlap : overlaps)
	{
		elementActor = Cast<AElementActor>(overlap.Actor.Get());

		if (!elementActor)
			continue;

		if (!elementActor->GetRootComponent()->IsVisible())
			continue;

		break;
	}

	if( elementActor )
	{
		USceneComponent * parent = elementActor->GetRootComponent()->GetAttachParent();

		TArray<USceneComponent*> children;
		parent->GetChildrenComponents(false, children);

		for (USceneComponent * child : children)
		{
			AElementActor * childElementActor = Cast<AElementActor>(child->GetOwner());

			if (!childElementActor)
				continue;

			if (childElementActor->groupName != elementActor->groupName)
				continue;

			newSelection.Add(childElementActor);
		}

		clearSelection();

		*_selection = newSelection;

		_showSelection(*_selection);
	}
}

void USelectionTool::_tickSelection(float dt, UPrimitiveComponent * hand)
{
	// perform an overlap test at the hand position to find actors that overlap
	TArray<FOverlapResult> overlaps;

	FCollisionShape collisionShape = FCollisionShape::MakeSphere(_brushRadius());

	FCollisionQueryParams params;

	UActorComponent * root;
	if (discreteEditorComponent)
		root = discreteEditorComponent;
	else
		root = regionGrowingComponent;

	params.AddIgnoredActor(root->GetOwner());
	params.AddIgnoredActor(hand->GetOwner());

	bool overlap = root->GetWorld()->OverlapMultiByChannel(
		overlaps,
		hand->GetComponentLocation(),
		hand->GetComponentRotation().Quaternion(),
		ECollisionChannel::ECC_WorldDynamic,
		collisionShape,
		params
	);

	// update the selection
	for (FOverlapResult overlap : overlaps)
	{
		AElementActor * elementActor = Cast<AElementActor>(overlap.Actor.Get());

		if (!elementActor)
			continue;

		if (!elementActor->GetRootComponent()->IsVisible())
			continue;

		if (selectionMode == ESelectionToolMode::Adding)
		{
			bool inSet = false;
			_selection->Add(elementActor, &inSet);


			if (!inSet)
				elementActor->showSelectionOutline();
		}
		else if (selectionMode == ESelectionToolMode::Removing)
		{
			bool wasInSet = _selection->Contains(elementActor);

			_selection->Remove(elementActor);

			if (wasInSet)
				elementActor->hideSelectionOutline();
		}
	}
}

void USelectionTool::faceDown_released()
{
	auto mode = selectionMode == ESelectionToolMode::Adding ? ESelectionToolMode::Removing : ESelectionToolMode::Adding;

	setSelectionMode( mode );
}

void USelectionTool::faceDown_touchStart()
{

}

void USelectionTool::faceDown_touchEnd()
{

}

void USelectionTool::faceUp_released(USceneComponent * interactionPoint /*= nullptr*/)
{

}

void USelectionTool::setSelectionMode( ESelectionToolMode mode )
{
	auto oldMode = selectionMode;

	selectionMode = mode;

	if(widgetComponent)
	{
		UUserWidget * widget = widgetComponent->GetUserWidgetObject();

		if(!widget || !widget->Implements<USelectionToolDelegate>())
			return;

		ISelectionToolDelegate::Execute_didChangeSelectionMode( widget, mode, oldMode );
	}

	if(selectionPointWidgetComponent)
	{
		UUserWidget * widget = selectionPointWidgetComponent->GetUserWidgetObject();

		if(!widget || !widget->Implements<USelectionToolDelegate>())
			return;

		ISelectionToolDelegate::Execute_didChangeSelectionMode( widget, mode, oldMode );
	}
}

void USelectionTool::gainFocus()
{
	if( widgetComponent && widgetClass)
	{
		widgetComponent->SetWidgetClass( widgetClass );

		UUserWidget * widget = widgetComponent->GetUserWidgetObject();

		if(!widget->Implements<USelectionToolDelegate>())
			return;

		ISelectionToolDelegate::Execute_didChangeSelectionMode( widget, selectionMode, selectionMode );
	}

	if(selectionPointWidgetComponent && selectionPointWidgetClass)
	{
		selectionPointWidgetComponent->SetWidgetClass( selectionPointWidgetClass );

		UUserWidget * widget = selectionPointWidgetComponent->GetUserWidgetObject();

		if(!widget->Implements<USelectionToolDelegate>())
			return;

		ISelectionToolDelegate::Execute_didChangeSelectionMode( widget, selectionMode, selectionMode );
	}

	if (!useGroupSelection)
	{
		clearSelection();
	}
}

std::vector<AElementActor*> USelectionTool::_toVector( TSet<AElementActor*> aSet )
{
	std::vector<AElementActor*> result;

	for(AElementActor * actor : aSet)
		result.push_back( actor );

	return result;
}

void USelectionTool::loseFocus()
{
	std::vector<AElementActor*> vecA = _toVector( selectionA );

	if(regionGrowingComponent)
		regionGrowingComponent->updateExampleSelection( vecA, weightA );

	if (generatorBrushTool && generatorBrushTool->generator())
		generatorBrushTool->generator()->setSelection(vecA);

	// hide the brush mesh component
	if(_selectionMeshComponent)
		_selectionMeshComponent->SetRelativeScale3D( FVector::ZeroVector );
}

void USelectionTool::_createSelectionMeshComponent( UPrimitiveComponent * selectionPoint )
{
	if(_selectionMeshComponent)
		return;

	_selectionMeshComponent = NewObject<UStaticMeshComponent>( selectionPoint );
	_selectionMeshComponent->AttachToComponent( selectionPoint, FAttachmentTransformRules::KeepRelativeTransform );
	_selectionMeshComponent->SetStaticMesh( selectionMesh );
	_selectionMeshComponent->SetMaterial( 0, selectionMeshMaterial );
	_selectionMeshComponent->SetRelativeScale3D( FVector( _brushRadius() * selectionMeshScaleFactor ) );


	_selectionMeshComponent->RegisterComponent();
}

void USelectionTool::_destroySelectionMeshComponent()
{
	if(_selectionMeshComponent == nullptr || !_selectionMeshComponent->IsValidLowLevel())
		return;

	_selectionMeshComponent->DestroyComponent();

	_selectionMeshComponent = nullptr;
}

void USelectionTool::makeSelectionAActive()
{
	_hideSelection( selectionB );

	_selection = &selectionA;

	_showSelection( selectionA );
}

void USelectionTool::makeSelectionBActive()
{
	_hideSelection( selectionA );

	_selection = &selectionB;

	_showSelection( selectionB );
}

void USelectionTool::_hideSelection( TSet<AElementActor*>& selected )
{
	for(AElementActor * elementActor : selected)
	{
		elementActor->hideSelectionOutline();
	}
}

void USelectionTool::_showSelection( TSet<AElementActor*>& selected )
{
	for(AElementActor * elementActor : selected)
	{
		elementActor->showSelectionOutline();
	}
}

void USelectionTool::clearSelection()
{
	_hideSelection( *_selection );

	_selection->Empty();
}