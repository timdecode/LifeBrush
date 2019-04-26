// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "VRTool.h"
#include "RegionGrowingToolInit.h"

#include "Algorithm/Algorithm.h"

#include "SelectionTool.generated.h"

class URegionGrowingComponent;
class AElementActor;
class UDiscreteElementEditorComponent;
class URegionGrowingGeneratorTool;

UENUM( BlueprintType )
enum class ESelectionToolMode : uint8
{
	Adding UMETA( DisplayName = "Adding" ),
	Removing UMETA( DisplayName = "Removing" ),
};

UINTERFACE( BlueprintType )
class USelectionToolDelegate : public UInterface
{
	GENERATED_BODY()
};

class ISelectionToolDelegate
{
	GENERATED_BODY()

public:
	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void setSelectionTool( USelectionTool * tool );

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void didChangeSelectionMode( ESelectionToolMode newMode, ESelectionToolMode oldMode );

	UFUNCTION( BlueprintCallable, BlueprintImplementableEvent, Category = "LifeBrush" )
	void didChangeModifier( ESelectionToolMode newMode, ESelectionToolMode oldMode );
};

UCLASS( Blueprintable )
class USelectionTool : public UTool
{
	GENERATED_BODY()

public:
	virtual ~USelectionTool();

public:
	UFUNCTION( BlueprintCallable, Category = "LifeBrush" )
	void makeSelectionAActive();

	UFUNCTION( BlueprintCallable, Category = "LifeBrush" )
	void makeSelectionBActive();

	UFUNCTION( BlueprintCallable, Category = "LifeBrush" )
	void clearSelection();

public:

	UPROPERTY()
	URegionGrowingComponent * regionGrowingComponent = nullptr;

	UPROPERTY()
	UDiscreteElementEditorComponent * discreteEditorComponent = nullptr;

	UPROPERTY()
	URegionGrowingGeneratorTool * generatorBrushTool = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool useGroupSelection = false;

	UPROPERTY()
	TSet<AElementActor*> selectionA;

	UPROPERTY()
	TSet<AElementActor*> selectionB;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float weightA = 1.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float weightB = 0.0f;


	TSet<AElementActor*> * _selection = &selectionA;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	UStaticMesh * selectionMesh;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	UMaterialInterface * selectionMeshMaterial;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float selectionMeshScaleFactor = 0.01f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement = "SelectionToolDelegate") )
	TSubclassOf<class UUserWidget> widgetClass;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (MustImplement = "SelectionToolDelegate") )
	TSubclassOf<class UUserWidget> selectionPointWidgetClass;

	virtual TSubclassOf<class UUserWidget> getWidgetClass() { return widgetClass; }

	virtual TSubclassOf<class UUserWidget> getSelectionWidgetClass() { return selectionPointWidgetClass; }

protected:
	UPROPERTY()
	UStaticMeshComponent * _selectionMeshComponent;

	UPROPERTY()
	ESelectionToolMode selectionMode = ESelectionToolMode::Adding;



public:
	void init(FRGC_UToolInitProperties& initProperties)
	{
		UTool::init(initProperties);

		regionGrowingComponent = initProperties.regionGrowingComponent;

		makeSelectionAActive();
	}

	// we have this weird postInit method because we need to init this component with a region growing component, but when
	// we send this guy an init, the RGC isn't ready yet. Therefore, we have this postInit, that is called on the VRSketchyPawn tick.
	// This is an ugly dirty hack.
	void postInit();
	
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

	virtual void faceDown_released() override;
	virtual void faceUp_released(USceneComponent * interactionPoint = nullptr) override;

	virtual void faceDown_touchStart() override;
	virtual void faceDown_touchEnd() override;

	virtual void gainFocus() override;

	virtual void loseFocus() override;

	UFUNCTION( BlueprintCallable, Category = LifeBrush )
	void setSelectionMode( ESelectionToolMode mode );

	UFUNCTION( BlueprintCallable, Category = LifeBrush )
	ESelectionToolMode getSelectionMode() { return selectionMode; }

protected:
	float _brushRadius();

	void _createSelectionMeshComponent( UPrimitiveComponent * selectionPoint );
	void _destroySelectionMeshComponent();

	void _tickGroupSelection(float dt, UPrimitiveComponent * hand);
	void _tickSelection(float dt, UPrimitiveComponent * hand);


	void _hideSelection( TSet<AElementActor*>& selected );
	void _showSelection( TSet<AElementActor*>& selected );

	std::vector<AElementActor*> _toVector( TSet<AElementActor*> aSet );
};