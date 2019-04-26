// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "VRTool.h"

#include "RegionGrowingToolInit.h"


#include "RegionGrowingTool.generated.h"

UENUM( BlueprintType )
enum class EGenerativeTickMode : uint8
{
	Generating,
	Erasing,
};

UENUM( BlueprintType )
enum class EGenerativeDrawMode : uint8
{
	Volumetric,
	Surface
};

UENUM(BlueprintType)
enum class EAutomaticDrawMode : uint8
{
	Disabled,
	ByNearestSurface,
	ByElementHint
};

UINTERFACE( BlueprintType )
class UGenerativeBrushToolDelegate : public UInterface
{
	GENERATED_BODY()
};

class IGenerativeBrushToolDelegate
{
	GENERATED_BODY()

public:
	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void setGenerativeBrush( class UGenerativeBrushTool * brush );

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void didChangeDrawMode(EGenerativeDrawMode newMode, EGenerativeDrawMode oldMode );

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void didChangeTickMode( EGenerativeTickMode newMode, EGenerativeTickMode oldMode );
};

UCLASS( Blueprintable )
class UGenerativeBrushTool : public UBrushTool
{
	GENERATED_BODY()

public:
	virtual ~UGenerativeBrushTool();

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Automatic Draw Mode")
	EAutomaticDrawMode automaticDrawMode = EAutomaticDrawMode::Disabled;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Automatic Draw Mode")
	float minSurfaceDistanceThreshold = 10.0f; // If we start painting within this minimum distance, we will start painting in surface mode.

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement="GenerativeBrushToolDelegate") )
	TSubclassOf<class UUserWidget> widgetClass;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement = "GenerativeBrushToolDelegate") )
	TSubclassOf<class UUserWidget> selectionPointWidgetClass;

	virtual TSubclassOf<class UUserWidget> getWidgetClass() { return widgetClass; }

	virtual TSubclassOf<class UUserWidget> getSelectionWidgetClass() { return selectionPointWidgetClass; }

public:
	UFUNCTION( BlueprintCallable, Category = Generation )
	void setDrawMode(EGenerativeDrawMode mode );

	UFUNCTION( BlueprintCallable, Category = Generation )
	EGenerativeDrawMode getDrawMode() { return _drawMode; }

	UFUNCTION( BlueprintCallable, Category = Generation )
	void setTickMode( EGenerativeTickMode mode );

	UFUNCTION( BlueprintCallable, Category = Generation )
	EGenerativeTickMode getTickMode() { return _tickMode; }

public:
	void init(FRGC_UToolInitProperties& initProperties)
	{
		UTool::init(initProperties);
	}

	virtual void gainFocus() override;

	virtual void twoHandStart( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) override {}
	virtual void twoHandEnd( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) override {}

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform transformA,
		FTransform transformB
	) override {}

	virtual void faceDown_released() override;
	virtual void faceUp_released(USceneComponent * interactionPoint = nullptr) override;

protected:
	EGenerativeDrawMode _drawMode = EGenerativeDrawMode::Volumetric;
	EGenerativeTickMode _tickMode = EGenerativeTickMode::Generating;
};

UCLASS(Blueprintable)
class URegionGrowing_GenerativeBrushTool : public UGenerativeBrushTool
{
	GENERATED_BODY()

public:
	virtual ~URegionGrowing_GenerativeBrushTool();

	URegionGrowingComponent * regionGrowingComponent = nullptr;

public:
	void init(FRGC_UToolInitProperties& initProperties)
	{
		UTool::init(initProperties);

		regionGrowingComponent = initProperties.regionGrowingComponent;
	}

	virtual void gainFocus() override;
	virtual void loseFocus() override;

	virtual void oneHandStart( UPrimitiveComponent * hand ) override;
	virtual void oneHandEnd( UPrimitiveComponent * hand ) override;

	virtual void twoHandStart( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) override {}
	virtual void twoHandEnd( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) override {}

	virtual void tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform ) override;

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform transformA,
		FTransform transformB
	) override {}


protected:
	void _tickOneHand_generate( float dt, UPrimitiveComponent * hand, FTransform lastTransform );
	void _tickOneHand_erase( float dt, UPrimitiveComponent * hand, FTransform lastTransform );
};
