// Copyright 2019, Timothy Davison. All rights reserved.

#pragma once

#include "GameFramework/Pawn.h"

#include "VRTool.h"
#include "ShipEditorSimulation/GraphSnapshot.h"

#include "VRSketchyPawn.generated.h"


class UMotionControllerComponent;
class URegionGrowingComponent;
class UFlexSimulationComponent;
class USelectionTool;
class UGenerativeBrushTool;
class UDigTool;
class UMeshCollectionTool;
class UEventVisualizationTool;
class UAgentPathlineTool;
class UPhysicalInteractionTool;
class UDiscreteElementEditorComponent;
class URegionGrowingGeneratorTool;
class UStringGeneratorTool;
class UCollagenGeneratorTool;
class USwarmGeneratorTool;
class URegionGrowing_GenerativeBrushTool;
class UExemplarInspectorTool;
struct FFlexSimulation;

UENUM(BlueprintType)
enum class ESketchyInteractionMode : uint8
{
	Painting UMETA(DisplayName = "Painting"),
	Simulating UMETA(DisplayName = "Simulating"),
};

UCLASS( DefaultToInstanced )
class LIFEBRUSH_API AVRSketchyPawn : public APawn, public UToolDelegate
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	bool developerMode = true;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	USceneComponent * vrScene;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UCameraComponent * camera;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UMotionControllerComponent * leftController;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UMotionControllerComponent * rightController;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UStaticMeshComponent * rightInteractionPoint;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UStaticMeshComponent * leftInteractionPoint;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UWidgetComponent * rightPadWidget;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UWidgetComponent * rightSelectionPointWidget;

	UPROPERTY( EditInstanceOnly, BlueprintReadWrite, Category = "ShipEditor" )
	AActor * regionGrowingActor = nullptr;

	UPROPERTY( EditInstanceOnly, BlueprintReadWrite, Category = "ShipEditor" )
	AActor * flexSimulationActor = nullptr;

	UPROPERTY(EditInstanceOnly, BlueprintReadWrite, Category = "ShipEditor")
	AActor * discreteElementEditorActor = nullptr;

	UPROPERTY(EditInstanceOnly, BlueprintReadWrite, Category = "ShipEditor")
	AActor * exemplarActor = nullptr;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor" )
	URegionGrowing_GenerativeBrushTool * generativeBrushTool;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor" )
	URegionGrowingGeneratorTool * regionGrowingElementGeneratorTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UStringGeneratorTool * stringGeneratorTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UCollagenGeneratorTool * collagenGeneratorTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	USwarmGeneratorTool * swarmGeneratorTool;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor" )
	USelectionTool * selectionTool;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor" )
	UDigTool * digTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UMeshCollectionTool * meshCollectionTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UEventVisualizationTool * eventVisualizationTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UAgentPathlineTool * agentPathlineTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UPhysicalInteractionTool * physicalInteractionTool;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "ShipEditor")
	UExemplarInspectorTool * exemplarInspectorTool;

public:
	UFUNCTION(BlueprintCallable, Category = Generation)
	ESketchyInteractionMode interactionMode() { return _interactionMode; }

	UFUNCTION(BlueprintCallable, Category = Generation)
	void setInteractionMode( ESketchyInteractionMode mode );

	UFUNCTION( BlueprintCallable, Category = Generation )
	bool isSimulating();

	UFUNCTION( BlueprintCallable, Category = Generation )
	void setSimulating( bool simulating );

	UFUNCTION( BlueprintCallable, Category = Generation )
	void leftTrigger( float value );

	UFUNCTION( BlueprintCallable, Category = Generation )
	void rightTrigger( float value );

	UFUNCTION( BlueprintCallable, Category = Generation )
	void setCurrentTool(UTool * newCurrentTool);

	UFUNCTION( BlueprintCallable, Category = Generation )
	UTool * getCurrentTool() { return _currentTool; }

	// --------------------
	// Simulation snapshots
	// --------------------
	// Experimental: Snapshot the simulation state, to later be restored with restoreSimulation. This stores
	// the snapshot in memory, not on disk.
	UFUNCTION(BlueprintCallable, Category = Generation)
	void snapshotSimulation();

	// Restore a previously snapshotted simulation.
	UFUNCTION(BlueprintCallable, Category = Generation)
	void restoreSimulation();

	UFUNCTION(BlueprintNativeEvent)
	void ShowToolSelectPopup();
	virtual void ShowToolSelectPopup_Implementation();

	// Snapshots the complete state of the simulation back to the editor world.
	UFUNCTION(BlueprintCallable, Category = Generation)
	void snapshotSimulationStateToEditor();

	// Snapshots the runtime-mesh component of the simulation to the passed actor. This is
	// useful for capturing a chunked-volume mesh created with the UDigTool.
	void snapshotMeshInterfaceRuntimeMeshComponentToActor(AActor* actor);

	// Snapshots just the graphics of a simulation, such as RMCs and the ISMCs, to the editor
	// world. This is useful for creating fancy screen shots later on. The simulation state is
	// not capture.
	UFUNCTION(BlueprintCallable, Category = Generation)
	void snapshotSimulationGraphicsToEditor();

protected:
	FGraphSnapshot _snapshot;

protected:
	URegionGrowingComponent * regionGrowingComponent;
	UFlexSimulationComponent * flexComponent;
	FFlexSimulation * flexSimulation;
	UDiscreteElementEditorComponent * editorComponent;

	bool _didInit = false;
	bool _didInitSimulation = false;

	UPROPERTY(EditInstanceOnly, BlueprintReadWrite, Category = "ShipEditor")
	ESketchyInteractionMode	_interactionMode = ESketchyInteractionMode::Painting;

	UTool * _currentTool;

	bool _leftTriggerDown = false;
	bool _rightTriggerDown = false;

	float _lastLeftTriggerValue = 0.0f;
	float _lastRightTriggerValue = 0.0f;

	bool _wantsLeftSpiderMan = false;
	bool _wantsRightSpiderMan = false;


	FVector _lastLeft;
	FVector _lastRight;


	int _simulationSnapShotCount = 0;

	bool _leftTouchActive = false;
	bool _rightTouchActive = false;

public:
	// Sets default values for this pawn's properties
	AVRSketchyPawn();

	virtual void PreInitializeComponents() override;

	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void Tick( float DeltaSeconds ) override;

	// Called to bind functionality to input
	virtual void SetupPlayerInputComponent(class UInputComponent* InputComponent) override;

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManLeftStart();
	
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManLeftEnd();
	
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManLeftUpdate();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManRightStart();
	
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManRightEnd();
	
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void spiderManRightUpdate();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void takeGraphicalSnapshotAndHighResShot();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void takeHighReshShot();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void snapshotElementDomain();

	void leftController_touchStart();
	void leftController_touchUpdated();
	void leftController_touchEnd();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void leftController_upFace_released();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void leftController_downFace_released();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void leftController_leftFace_released();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void leftController_rightFace_released();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_touchStart();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_touchUpdated();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_touchEnd();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceUp_pressed();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceUp_released();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceDown_pressed();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceDown_released();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceLeft_pressed();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceLeft_released();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceRight_pressed();
	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_faceRight_released();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void rightController_shoulder_released();

protected:
	FVector2D _getLeftTouchPoint();
	FVector2D _getRightTouchPoint();

	void _initSimulation_oneTime();
	void _initSimulationBounds();
	void _initTools();

	FString _actorLabelByDate(FString baseName);



	void _startSimulation();
	void _endSimulation();

	void _startCurrentToolLeftTrigger();
	void _startCurrentToolRightTrigger();


public: // UTool delegate
	virtual void cedeFocus(UTool * tool) override;

};
