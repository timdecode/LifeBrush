//
//  Created by Timothy Davison on 2018-03-06, based on my code for Ship Editor.
//  Copyright (c) Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "RegionGrowingComponent.h"

#include "MotionControllerComponent.h"
#include "WidgetComponent.h"

#include "VRTool.h"

void UTool::loseFocus()
{
	_hideWidgets();
}

void UTool::_hideWidgets()
{
	if(widgetComponent)
	{
		widgetComponent->SetWidgetClass( nullptr );
		widgetComponent->SetVisibility( false );
	}

	if(selectionPointWidgetComponent)
	{
		selectionPointWidgetComponent->SetWidgetClass( nullptr );
		selectionPointWidgetComponent->SetVisibility( false );
	}
}

void UTool::_loadWidgets()
{
	if( widgetComponent && getWidgetClass() )
	{
		widgetComponent->SetWidgetClass( getWidgetClass() );
		widgetComponent->SetVisibility( true );
	}

	if( selectionPointWidgetComponent && getSelectionWidgetClass() )
	{
		selectionPointWidgetComponent->SetWidgetClass( getSelectionWidgetClass() );
		selectionPointWidgetComponent->SetVisibility( true );
	}
}



// ----------------------------------------------------------
// UBrushTool
// ----------------------------------------------------------

UBrushTool::~UBrushTool()
{
	if (_brushMeshComponent && _brushMeshComponent->IsValidLowLevel())
	{
		_brushMeshComponent->DestroyComponent();
		_brushMeshComponent = nullptr;
	}
}

void UBrushTool::loseFocus()
{
	if (_brushMeshComponent && _brushMeshComponent->IsValidLowLevel())
	{
		_brushMeshComponent->DestroyComponent();
		_brushMeshComponent = nullptr;
	}
}

void UBrushTool::oneHandStart(UPrimitiveComponent * hand)
{
	_createBrushMeshComponent(hand);

	_brushMeshComponent->SetVisibility(shouldShowBrush());
}

void UBrushTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_brushMeshComponent)
		_brushMeshComponent->SetVisibility(false);
}

void UBrushTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform)
{
	_updateBrushMesh();
}

void UBrushTool::_createBrushMeshComponent(UPrimitiveComponent * selectionPoint)
{
	if (_brushMeshComponent)
		return;

	AActor * actor = selectionPoint->GetOwner();

	_brushMeshComponent = NewObject<UStaticMeshComponent>(actor);

	_brushMeshComponent->AttachToComponent(selectionPoint, FAttachmentTransformRules::KeepRelativeTransform);
	_brushMeshComponent->SetStaticMesh(brushMesh);
	_brushMeshComponent->SetMaterial(0, brushMeshMaterial);
	_brushMeshComponent->SetVisibility(true);
	_brushMeshComponent->SetRelativeScale3D(FVector(_brushRadius() * brushMeshScaleFactor));

	_brushMeshComponent->RegisterComponent();
}

float UBrushTool::_brushRadius()
{
	float delta = brushMaxRadius - brushMinRadius;

	return brushMinRadius + delta * selectionATriggerValue();
}

void UBrushTool::_updateBrushMesh()
{
	if (_brushMeshComponent)
	{
		float radius = _brushRadius();

		_brushMeshComponent->SetVisibility(shouldShowBrush());
		_brushMeshComponent->SetRelativeScale3D(FVector(radius * brushMeshScaleFactor));
	}
}
