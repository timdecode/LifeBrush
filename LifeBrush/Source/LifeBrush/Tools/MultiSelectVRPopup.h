// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "WidgetComponent.h"

#include "MultiSelectVRPopup.generated.h"

UINTERFACE( BlueprintType )
class UMultiSelectPopupDelegate : public UInterface
{
	GENERATED_BODY()
};

class IMultiSelectPopupDelegate
{
	GENERATED_BODY()

public:
	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void setSelectionTool( USelectionTool * tool );
};

UCLASS( DefaultToInstanced )
class LIFEBRUSH_API AMultiSelectVRPopup : public AActor
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor", meta = (MustImplement = "MultiSelectPopupDelegate") )
	TSubclassOf<class UUserWidget> widgetClass;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor" )
	UWidgetComponent * widgetComponent;

	void setSelectionTool( USelectionTool * tool );



public:
	AMultiSelectVRPopup();

	virtual void BeginPlay() override;

protected:
	UPROPERTY()
	USelectionTool * selectionTool;
};
