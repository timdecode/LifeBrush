// Copyright 2017 Timothy Davison, all rights reserved.

#pragma once

#include "VRTool.h"
#include "RegionGrowingToolInit.h"

#include "VolumeComponent.h"

#include "DigTool.generated.h"

UENUM( BlueprintType )
enum class EDigMode : uint8
{
	Adding UMETA( DisplayName = "Adding" ),
	Removing UMETA( DisplayName = "Removing" ),
	Smoothing UMETA( DisplayName = "Smoothing" )
};

UINTERFACE( BlueprintType )
class UDigToolDelegate : public UInterface
{
	GENERATED_BODY()
};

class IDigToolDelegate
{
	GENERATED_BODY()

public:

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void didChangeDigMode( EDigMode newMode, EDigMode oldMode );

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void setDigTool( class UDigTool * digTool );
};

/**
 * 
 */
UCLASS( Blueprintable )
class LIFEBRUSH_API UDigTool : public UTool
{
	GENERATED_BODY()

public:

	void init(FRGC_UToolInitProperties& initProperties);

	virtual void oneHandStart( UPrimitiveComponent * hand ) override;
	virtual void oneHandEnd( UPrimitiveComponent * hand ) override;

	virtual void tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform ) override;

	virtual void focused() override;
	virtual void loseFocus() override;

	virtual void faceDown_released() override;
	virtual void faceUp_released( USceneComponent * interactionPoint = nullptr ) override;
	virtual void faceLeft_released() override;
	virtual void faceRight_released() override;

	void setDigMode( EDigMode mode );
	EDigMode getDigMode() { return _digMode; }

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float radius = 5.0f; // cell units

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float fillRate = 1.0f; 

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float maxValue = 100.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UStaticMesh * selectionMesh;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UMaterialInterface * selectionMeshMaterial;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float selectionMeshScaleFactor = 2.0f;

	UPROPERTY()
	URegionGrowingComponent * regionGrowingComponent = nullptr;

	UPROPERTY()
	UDiscreteElementEditorComponent * editorComponent = nullptr;

	// ------------------------------
	// Widget Stuff
	// ------------------------------
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor", meta = (MustImplement = "DigToolDelegate") )
	TSubclassOf<class UUserWidget> widgetClass;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor", meta = (MustImplement = "DigToolDelegate") )
	TSubclassOf<class UUserWidget> selectionPointWidgetClass;

	virtual TSubclassOf<class UUserWidget> getWidgetClass() { return widgetClass; }

	virtual TSubclassOf<class UUserWidget> getSelectionWidgetClass() { return selectionPointWidgetClass; }

protected:
	void _dig( class UChunkedVolumeComponent * volume, FVector volumePosition, float radius, float dt);
	void _smooth(UChunkedVolumeComponent * volume, FVector volumePosition, float volumeRadius );

	void _createSelectionMeshComponent( UPrimitiveComponent * selectionPoint );

	void _destroySelectionMeshComponent();

	float _brushRadius();

	void _convolve(ChunkGrid<float>& grid, FIntVector& index );

	lb::BasicGrid<float> _gaussianKernel( size_t n );

	void _notifyDidModeDelegates( EDigMode newMode, EDigMode oldMode );

protected:
	UPROPERTY()
	EDigMode _digMode = EDigMode::Adding;

	UPROPERTY()
	UStaticMeshComponent * _selectionMeshComponent;

	lb::BasicGrid<float> _kernel;

	UChunkedVolumeComponent * _lastVolume = nullptr;

private:
	void _notifyDigToolDelegates();
};
