// Copyright 2017 Timothy Davison, all rights reserved.

#pragma once

#include "VRTool.h"
#include "RegionGrowingToolInit.h"

#include "VolumeComponent.h"

#include "DigTool.generated.h"

class UFlexSimulationComponent;

UCLASS(Blueprintable)
class LIFEBRUSH_API ADigToolPopupActor : public AActor
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	class UDigTool * digTool = nullptr;
};

UENUM(BlueprintType)
enum class EDigShape : uint8
{
	Sphere UMETA(DisplayName = "Sphere"),
	Capsule UMETA(DisplayName = "Capsule"),
};

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

	void init(FRGC_UToolInitProperties& initProperties, UCameraComponent * camera);

	virtual void oneHandStart( UPrimitiveComponent * hand ) override;
	virtual void oneHandEnd( UPrimitiveComponent * hand ) override;

	virtual void tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastTransform ) override;
	virtual void tickTwoHand(float dt, UPrimitiveComponent * handA, UPrimitiveComponent * handB, FTransform lastTransformA, FTransform lastTransformB) override;

	virtual void gainFocus() override;
	virtual void loseFocus() override;

	virtual void faceDown_released() override;
	virtual void faceUp_released( USceneComponent * interactionPoint = nullptr ) override;
	virtual void faceLeft_released() override;
	virtual void faceRight_released() override;

	void setDigMode( EDigMode mode );
	EDigMode getDigMode() { return _digMode; }


	virtual bool consume_rightShoulder_pressed() override;



public:
	UFUNCTION(BlueprintCallable, Category = Generation)
	void SaveMeshToCollection();

	UFUNCTION(BlueprintCallable, Category = Generation)
	void SaveMeshToScene();

	UFUNCTION(BlueprintCallable, Category = Generation)
	void ExitSession();

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float radius = 5.0f; // cell units

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float fillRate = 1.0f; 

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxValue = 100.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	EDigShape digShape = EDigShape::Sphere;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float capsuleLength = 10.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	UStaticMesh * selectionMesh;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	UMaterialInterface * selectionMeshMaterial;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float selectionMeshScaleFactor = 2.0f;

	
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TSubclassOf<class ADigToolPopupActor> popupActorClass;

	UPROPERTY()
	ADigToolPopupActor * popupActor = nullptr;

	UPROPERTY()
	UFlexSimulationComponent * flexComponent = nullptr;

	// ------------------------------
	// Widget Stuff
	// ------------------------------
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement = "DigToolDelegate") )
	TSubclassOf<class UUserWidget> widgetClass;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement = "DigToolDelegate") )
	TSubclassOf<class UUserWidget> selectionPointWidgetClass;

	virtual TSubclassOf<class UUserWidget> getWidgetClass() { return widgetClass; }

	virtual TSubclassOf<class UUserWidget> getSelectionWidgetClass() { return selectionPointWidgetClass; }

protected:
	void _dig( class UChunkedVolumeComponent * volume, FVector volumePosition, float radius, float dt);
	void _capsuleDig(UChunkedVolumeComponent * volume, FVector volumeStart, FVector volumeEnd, float volumeRadius, float dt);
	void _smooth(UChunkedVolumeComponent * volume, FVector volumePosition, float volumeRadius);

	void _createSelectionMeshComponent( UPrimitiveComponent * selectionPoint );

	void _destroySelectionMeshComponent();

	float _brushRadius();

	void _convolve(ChunkGrid<float>& grid, FIntVector& index );

	lb::BasicGrid<float> _gaussianKernel( size_t n );

	void _notifyDidModeDelegates( EDigMode newMode, EDigMode oldMode );

	FString _actorLabelByDate(FString baseName);

	void _setRMCVisibility(bool visibile);

protected:
	UPROPERTY()
	EDigMode _digMode = EDigMode::Adding;

	UPROPERTY()
	UStaticMeshComponent * _selectionMeshComponent;

	lb::BasicGrid<float> _kernel;

	UCameraComponent * _camera = nullptr;

	UChunkedVolumeComponent * _lastVolume = nullptr;

private:
	void _notifyDigToolDelegates();
};
