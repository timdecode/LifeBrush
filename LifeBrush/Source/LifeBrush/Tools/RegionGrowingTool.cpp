// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "WidgetComponent.h"
#include "RegionGrowingComponent.h"

#include "RegionGrowingTool.h"

UGenerativeBrushTool::~UGenerativeBrushTool()
{
	_destroyBrushMeshComponent();
}

void UGenerativeBrushTool::focused()
{
	if(widgetComponent && widgetClass)
	{
		widgetComponent->SetWidgetClass( widgetClass );

		UUserWidget * widget = widgetComponent->GetUserWidgetObject();

		if(!widget->Implements<UGenerativeBrushToolDelegate>())
			return;

		IGenerativeBrushToolDelegate::Execute_didChangeTickMode( widget, _tickMode, _tickMode );
		IGenerativeBrushToolDelegate::Execute_didChangeDrawMode( widget, _drawMode, _drawMode );
	}

	if(selectionPointWidgetComponent && selectionPointWidgetClass)
	{
		selectionPointWidgetComponent->SetWidgetClass( selectionPointWidgetClass );

		UUserWidget * widget = selectionPointWidgetComponent->GetUserWidgetObject();

		if(!widget->Implements<UGenerativeBrushToolDelegate>())
			return;

		IGenerativeBrushToolDelegate::Execute_didChangeTickMode( widget, _tickMode, _tickMode );
		IGenerativeBrushToolDelegate::Execute_didChangeDrawMode( widget, _drawMode, _drawMode );
	}
}

void UGenerativeBrushTool::loseFocus()
{
	if( _brushMeshComponent )
		_brushMeshComponent->SetRelativeScale3D( FVector::ZeroVector );
}

void UGenerativeBrushTool::oneHandStart( UPrimitiveComponent * hand )
{
	_createBrushMeshComponent( hand );
}

void UGenerativeBrushTool::oneHandEnd( UPrimitiveComponent * hand )
{
}

float UGenerativeBrushTool::_brushRadius()
{
	return 2.0f + 6.0 * selectionATriggerValue();
}

void UGenerativeBrushTool::tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform )
{
	if(!_brushMeshComponent)
		return;

	float radius = _brushRadius();

	_brushMeshComponent->SetRelativeScale3D( FVector( radius * brushMeshScaleFactor ) );
}

void UGenerativeBrushTool::faceDown_released()
{
	setTickMode( _tickMode == EGenerativeTickMode::Generating ? EGenerativeTickMode::Erasing : EGenerativeTickMode::Generating );
}

void UGenerativeBrushTool::faceUp_released(USceneComponent * interactionPoint /*= nullptr*/)
{
	setDrawMode( _drawMode == EGenerativeDrawMode::Surface ? EGenerativeDrawMode::Volumetric : EGenerativeDrawMode::Surface );
}

void UGenerativeBrushTool::setDrawMode( EGenerativeDrawMode mode )
{
	auto oldMode = _drawMode;
	
	_drawMode = mode;

	if (!widgetComponent)
		return;

	UUserWidget * widget = widgetComponent->GetUserWidgetObject();

	if(!widget || !widget->Implements<UGenerativeBrushToolDelegate>())
		return;

	IGenerativeBrushToolDelegate::Execute_didChangeDrawMode( widget, mode, oldMode );
}

void UGenerativeBrushTool::setTickMode( EGenerativeTickMode mode )
{
	auto oldMode = _tickMode;

	// toggle the tick mode
	_tickMode = mode;

	if( widgetComponent)
	{
		UUserWidget * widget = widgetComponent->GetUserWidgetObject();

		if(!widget || !widget->Implements<UGenerativeBrushToolDelegate>())
			return;

		IGenerativeBrushToolDelegate::Execute_didChangeTickMode( widget, mode, oldMode );
	}

	if(selectionPointWidgetComponent)
	{
		UUserWidget * widget = selectionPointWidgetComponent->GetUserWidgetObject();

		if(!widget || !widget->Implements<UGenerativeBrushToolDelegate>())
			return;

		IGenerativeBrushToolDelegate::Execute_didChangeTickMode( widget, mode, oldMode );
	}
}



void UGenerativeBrushTool::_createBrushMeshComponent( UPrimitiveComponent * selectionPoint )
{
	if(_brushMeshComponent)
		return;

	AActor * actor = selectionPoint->GetOwner();

	_brushMeshComponent = NewObject<UStaticMeshComponent>( actor );
	
	_brushMeshComponent->AttachToComponent( selectionPoint, FAttachmentTransformRules::KeepRelativeTransform );
	_brushMeshComponent->SetStaticMesh( brushMesh );
	_brushMeshComponent->SetMaterial( 0, brushMeshMaterial );
	_brushMeshComponent->SetVisibility( true );
	_brushMeshComponent->SetRelativeScale3D( FVector( _brushRadius() * brushMeshScaleFactor ) );

	_brushMeshComponent->RegisterComponent();
}

void UGenerativeBrushTool::_destroyBrushMeshComponent()
{
	if(_brushMeshComponent == nullptr || !_brushMeshComponent->IsValidLowLevel())
		return;

	_brushMeshComponent->DestroyComponent();

	_brushMeshComponent = nullptr;
}























// -------------------------------------------------------------------------------
// URegionGrowing_GenerativeBrushTool
// -------------------------------------------------------------------------------

URegionGrowing_GenerativeBrushTool::~URegionGrowing_GenerativeBrushTool()
{
	_destroyBrushMeshComponent();
}

void URegionGrowing_GenerativeBrushTool::focused()
{
	Super::focused();
}

void URegionGrowing_GenerativeBrushTool::loseFocus()
{
	Super::loseFocus();
}

void URegionGrowing_GenerativeBrushTool::oneHandStart( UPrimitiveComponent * hand )
{
	Super::oneHandStart(hand);

	if (_tickMode == EGenerativeTickMode::Generating)
		regionGrowingComponent->startPaint();
	else if( _tickMode == EGenerativeTickMode::Erasing )
		regionGrowingComponent->startErase();

	// choose the draw-mode
	if (automaticDrawMode == EAutomaticDrawMode::ByNearestSurface)
	{
		FVector location = hand->GetComponentLocation();

		auto meshInterface = regionGrowingComponent->meshInterface();

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

void URegionGrowing_GenerativeBrushTool::oneHandEnd( UPrimitiveComponent * hand )
{
	Super::oneHandEnd(hand);

	regionGrowingComponent->endPaint();
	regionGrowingComponent->endErase();
}

void URegionGrowing_GenerativeBrushTool::tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform )
{
	Super::tickOneHand(dt, hand, lastTransform);

	if (!regionGrowingComponent)
		return;

	regionGrowingComponent->brushSize = _brushRadius();

	if(_tickMode == EGenerativeTickMode::Generating)
		_tickOneHand_generate( dt, hand, lastTransform );
	else
		_tickOneHand_erase( dt, hand, lastTransform );
}

void URegionGrowing_GenerativeBrushTool::_tickOneHand_generate( float dt, UPrimitiveComponent * hand, FTransform lastTransform )
{
	auto meshInterface = regionGrowingComponent->meshInterface();

	FVector location = hand->GetComponentLocation();

	if(_drawMode == EGenerativeDrawMode::Volumetric )
	{
		regionGrowingComponent->generationMode = EGenerationMode::SpacePainting;
		regionGrowingComponent->addBrushPoint( location );
	}
	else if( _drawMode == EGenerativeDrawMode::Surface )
	{
		auto nearest = meshInterface->nearestPointOnMesh( location );

		regionGrowingComponent->generationMode = EGenerationMode::SurfacePainting;
		regionGrowingComponent->addBrushPoint( nearest.point, nearest.surfaceIndex );
	}
}

void URegionGrowing_GenerativeBrushTool::_tickOneHand_erase( float dt, UPrimitiveComponent * hand, FTransform lastTransform )
{
	auto meshInterface = regionGrowingComponent->meshInterface();

	FVector location = hand->GetComponentLocation();

	float radius = _brushRadius();

	if (_drawMode == EGenerativeDrawMode::Volumetric)
	{
		regionGrowingComponent->eraseInRadiusAt( location, radius );
	}
	else if (_drawMode == EGenerativeDrawMode::Surface)
	{
		auto nearest = meshInterface->nearestPointOnMesh( location );

		regionGrowingComponent->eraseInRadiusAt( nearest.point, radius, nearest.surfaceIndex );
	}
}

