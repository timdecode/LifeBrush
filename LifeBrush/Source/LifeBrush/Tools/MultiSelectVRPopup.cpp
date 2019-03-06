// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "WidgetComponent.h"

#include "MultiSelectVRPopup.h"

void AMultiSelectVRPopup::setSelectionTool( USelectionTool * tool )
{
	selectionTool = tool;

	widgetComponent->SetWidgetClass( widgetClass );

	UUserWidget * widget = widgetComponent->GetUserWidgetObject();

	if(!widget || !widget->Implements<UMultiSelectPopupDelegate>())
		return;

	IMultiSelectPopupDelegate::Execute_setSelectionTool( widget, selectionTool );
}

AMultiSelectVRPopup::AMultiSelectVRPopup()
{
	RootComponent = CreateDefaultSubobject<USceneComponent>( TEXT( "Root" ) );

	widgetComponent = CreateDefaultSubobject<UWidgetComponent>( TEXT( "WidgetComponent" ) );

	widgetComponent->AttachToComponent( RootComponent, FAttachmentTransformRules::KeepRelativeTransform );
}

void AMultiSelectVRPopup::BeginPlay()
{
	Super::BeginPlay();
}

