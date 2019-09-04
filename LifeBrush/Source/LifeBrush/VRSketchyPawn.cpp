// Copyright 2016, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "MotionControllerComponent.h"
#include "RegionGrowingComponent.h"
#include "WidgetComponent.h"

#include "Simulation/FlexElements.h"

#include "Tools/RegionGrowingGeneratorTool.h"
#include "Tools/StringGeneratorTool.h"
#include "Tools/SelectionTool.h"
#include "Tools/RegionGrowingTool.h"
#include "Tools/DigTool.h"
#include "Tools/MeshCollectionTool.h"
#include "Tools/ExemplarInspectorTool.h"

#include "Visualization/EventVisualizationTool.h"

#include "ElementEditor/DiscreteElementEditorComponent.h"


#include "VRSketchyPawn.h"

#if WITH_EDITOR
#include "SimulationSnapshotActor.h"
#endif

// Sets default values
AVRSketchyPawn::AVRSketchyPawn()
{
 	// Set this pawn to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	RootComponent = CreateDefaultSubobject<USceneComponent>( TEXT( "Root" ) );
	RootComponent->SetMobility( EComponentMobility::Movable );

	// copied from ShipEditor 2016-11-16 Copyright Timothy Davison of Code Monkey Castle, all rights reserved.
	vrScene = CreateDefaultSubobject<USceneComponent>( TEXT( "VRScene" ) );
	vrScene->SetMobility( EComponentMobility::Movable );
	vrScene->AttachToComponent( RootComponent, FAttachmentTransformRules::KeepRelativeTransform );

	leftController = CreateDefaultSubobject<UMotionControllerComponent>( TEXT( "LeftController" ) );
	leftController->Hand = EControllerHand::Left;
	leftController->SetCollisionProfileName( TEXT( "BlockAll" ) );
	leftController->AttachToComponent( vrScene, FAttachmentTransformRules::KeepRelativeTransform );

	leftInteractionPoint = CreateDefaultSubobject<UStaticMeshComponent>( TEXT( "LeftSelectionPoint" ) );
	leftInteractionPoint->AttachToComponent( leftController, FAttachmentTransformRules::KeepRelativeTransform );

	rightController = CreateDefaultSubobject<UMotionControllerComponent>( TEXT( "RightController" ) );
	rightController->Hand = EControllerHand::Right;
	rightController->SetCollisionProfileName( TEXT( "BlockAll" ) );
	rightController->AttachToComponent( vrScene, FAttachmentTransformRules::KeepRelativeTransform );

	rightInteractionPoint = CreateDefaultSubobject<UStaticMeshComponent>( TEXT( "RightSelectionPoint" ) );
	rightInteractionPoint->AttachToComponent( rightController, FAttachmentTransformRules::KeepRelativeTransform );

	rightPadWidget = CreateDefaultSubobject<UWidgetComponent>( TEXT( "RightWidgetComponent" ) );
	rightPadWidget->AttachToComponent( rightController, FAttachmentTransformRules::KeepRelativeTransform );

	rightSelectionPointWidget = CreateDefaultSubobject<UWidgetComponent>( TEXT( "RightSelectionPointWidgetComponent" ) );
	rightSelectionPointWidget->AttachToComponent( rightInteractionPoint, FAttachmentTransformRules::KeepRelativeTransform );

	camera = CreateDefaultSubobject<UCameraComponent>( TEXT( "Camera" ) );
	camera->AttachToComponent( vrScene, FAttachmentTransformRules::KeepRelativeTransform );

	generativeBrushTool = CreateDefaultSubobject<URegionGrowing_GenerativeBrushTool>( TEXT( "RegionGrowingGenerativeBrushTool" ) );
	regionGrowingElementGeneratorTool = CreateDefaultSubobject<URegionGrowingGeneratorTool>(TEXT("RegionGrowingElementGenerator_BrushTool"));
	stringGeneratorTool = CreateDefaultSubobject<UStringGeneratorTool>(TEXT("StringGeneratorTool"));
	collagenGeneratorTool = CreateDefaultSubobject<UCollagenGeneratorTool>(TEXT("CollagenGeneratorTool"));
	swarmGeneratorTool = CreateDefaultSubobject<USwarmGeneratorTool>(TEXT("SwarmGeneratorTool"));
	selectionTool = CreateDefaultSubobject<USelectionTool>( TEXT( "SelectionTool" ) );
	digTool = CreateDefaultSubobject<UDigTool>( TEXT( "DigTool" ) );
	meshCollectionTool = CreateDefaultSubobject<UMeshCollectionTool>( TEXT("MeshCollectionTool") );
	eventVisualizationTool = CreateDefaultSubobject<UEventVisualizationTool>( TEXT("EventVisualizationTool") );
	agentPathlineTool = CreateDefaultSubobject<UAgentPathlineTool>(TEXT("AgentPathlineTool"));
	physicalInteractionTool = CreateDefaultSubobject<UPhysicalInteractionTool>(TEXT("PhysicalInteractionTool"));
	exemplarInspectorTool = CreateDefaultSubobject<UExemplarInspectorTool>(TEXT("ExemplarInspectorTool"));

	// necessary for VR 
	BaseEyeHeight = 0.0f;
}

void AVRSketchyPawn::PreInitializeComponents()
{
	Super::PreInitializeComponents();

	if (regionGrowingActor)
		regionGrowingComponent = regionGrowingActor->FindComponentByClass<URegionGrowingComponent>();

	if (discreteElementEditorActor)
	{
		editorComponent = discreteElementEditorActor->FindComponentByClass<UDiscreteElementEditorComponent>();
		editorComponent->camera = camera;
	}

	// tick after the RGC
	AddTickPrerequisiteComponent(regionGrowingComponent);
	AddTickPrerequisiteComponent(editorComponent);
	AddTickPrerequisiteComponent(flexComponent);
}

// Called when the game starts or when spawned
void AVRSketchyPawn::BeginPlay()
{
	Super::BeginPlay();
	
	if (flexSimulationActor)
	{
		flexComponent = flexSimulationActor->FindComponentByClass<UFlexSimulationComponent>();
	
		if( flexComponent ) flexSimulation = flexComponent->flexSimulation();
	}

	_initTools();

	_initSimulation_oneTime();

	if (_interactionMode == ESketchyInteractionMode::Painting && editorComponent)
		editorComponent->start();
}

// Called every frame
void AVRSketchyPawn::Tick( float DeltaTime )
{
	Super::Tick( DeltaTime );

	if(!_didInit)
	{
		selectionTool->postInit();

		_didInit = true;
	}

	if(_wantsLeftSpiderMan)
		spiderManLeftUpdate();

	if(_wantsRightSpiderMan)
		spiderManRightUpdate();

	if( _rightTouchActive )
		rightController_touchUpdated();

	if( _leftTouchActive ) 
		leftController_touchUpdated();

	if(_currentTool)
	{
		_currentTool->doTick( DeltaTime );
	}
}

// Called to bind functionality to input
void AVRSketchyPawn::SetupPlayerInputComponent(class UInputComponent* InputComponent)
{
	Super::SetupPlayerInputComponent(InputComponent);

	InputComponent->BindKey( EKeys::MotionController_Left_FaceButton1, IE_Released, this, &AVRSketchyPawn::leftController_upFace_released );
	InputComponent->BindKey( EKeys::MotionController_Left_FaceButton3, IE_Released, this, &AVRSketchyPawn::leftController_downFace_released );
	InputComponent->BindKey( EKeys::MotionController_Left_FaceButton4, IE_Released, this, &AVRSketchyPawn::leftController_leftFace_released );
	InputComponent->BindKey( EKeys::MotionController_Left_FaceButton2, IE_Released, this, &AVRSketchyPawn::leftController_rightFace_released );

	InputComponent->BindKey( EKeys::MotionController_Right_FaceButton1, IE_Released, this, &AVRSketchyPawn::rightController_faceUp_released );
	InputComponent->BindKey( EKeys::MotionController_Right_FaceButton3, IE_Released, this, &AVRSketchyPawn::rightController_faceDown_released );
	InputComponent->BindKey( EKeys::MotionController_Right_FaceButton4, IE_Released, this, &AVRSketchyPawn::rightController_faceLeft_released );
	InputComponent->BindKey( EKeys::MotionController_Right_FaceButton2, IE_Released, this, &AVRSketchyPawn::rightController_faceRight_released );

	InputComponent->BindKey(EKeys::MotionController_Right_FaceButton1, IE_Pressed, this, &AVRSketchyPawn::rightController_faceUp_pressed);
	InputComponent->BindKey(EKeys::MotionController_Right_FaceButton3, IE_Pressed, this, &AVRSketchyPawn::rightController_faceDown_pressed);
	InputComponent->BindKey(EKeys::MotionController_Right_FaceButton4, IE_Pressed, this, &AVRSketchyPawn::rightController_faceLeft_pressed);
	InputComponent->BindKey(EKeys::MotionController_Right_FaceButton2, IE_Pressed, this, &AVRSketchyPawn::rightController_faceRight_pressed);

	InputComponent->BindKey( EKeys::Steam_Touch_0, IE_Pressed, this, &AVRSketchyPawn::leftController_touchStart );
	InputComponent->BindKey( EKeys::Steam_Touch_0, IE_Released, this, &AVRSketchyPawn::leftController_touchEnd );

	InputComponent->BindKey( EKeys::Steam_Touch_1, IE_Pressed, this, &AVRSketchyPawn::rightController_touchStart );
	InputComponent->BindKey( EKeys::Steam_Touch_1, IE_Released, this, &AVRSketchyPawn::rightController_touchEnd );

	InputComponent->BindAxisKey( EKeys::MotionController_Left_Thumbstick_X );
	InputComponent->BindAxisKey( EKeys::MotionController_Left_Thumbstick_Y );

	InputComponent->BindAxisKey( EKeys::MotionController_Right_Thumbstick_X );
	InputComponent->BindAxisKey( EKeys::MotionController_Right_Thumbstick_Y );

	InputComponent->BindAxisKey( EKeys::MotionController_Left_TriggerAxis, this, &AVRSketchyPawn::leftTrigger );
	InputComponent->BindAxisKey( EKeys::MotionController_Right_TriggerAxis, this, &AVRSketchyPawn::rightTrigger );

	InputComponent->BindKey( EKeys::MotionController_Left_Grip1, IE_Pressed, this, &AVRSketchyPawn::spiderManLeftStart );
	InputComponent->BindKey( EKeys::MotionController_Left_Grip1, IE_Released, this, &AVRSketchyPawn::spiderManLeftEnd );

	InputComponent->BindKey( EKeys::MotionController_Right_Grip1, IE_Pressed, this, &AVRSketchyPawn::spiderManRightStart );
	InputComponent->BindKey( EKeys::MotionController_Right_Grip1, IE_Released, this, &AVRSketchyPawn::spiderManRightEnd );

	InputComponent->BindKey( EKeys::MotionController_Left_Shoulder, IE_Released, this, &AVRSketchyPawn::takeGraphicalSnapshotAndHighResShot );
	InputComponent->BindKey( EKeys::MotionController_Right_Shoulder, IE_Released, this, &AVRSketchyPawn::rightController_shoulder_released );
}



void AVRSketchyPawn::_initTools()
{
	auto regionGrowingRoot = discreteElementEditorActor ? 
		discreteElementEditorActor->GetRootComponent() :
		regionGrowingActor->GetRootComponent();

	// set the tool widgets
	generativeBrushTool->widgetComponent = rightPadWidget;
	generativeBrushTool->selectionPointWidgetComponent = rightSelectionPointWidget;

	regionGrowingElementGeneratorTool->widgetComponent = rightPadWidget;
	regionGrowingElementGeneratorTool->selectionPointWidgetComponent = rightSelectionPointWidget;

	stringGeneratorTool->widgetComponent = rightPadWidget;
	stringGeneratorTool->selectionPointWidgetComponent = rightSelectionPointWidget;

	collagenGeneratorTool->widgetComponent = rightPadWidget;
	collagenGeneratorTool->selectionPointWidgetComponent = rightSelectionPointWidget;

	swarmGeneratorTool->widgetComponent = rightPadWidget;
	swarmGeneratorTool->selectionPointWidgetComponent = rightSelectionPointWidget;
	
	selectionTool->widgetComponent = rightPadWidget;
	selectionTool->selectionPointWidgetComponent = rightSelectionPointWidget;
	selectionTool->discreteEditorComponent = editorComponent;
	selectionTool->generatorBrushTool = regionGrowingElementGeneratorTool;

	digTool->widgetComponent = rightPadWidget;
	digTool->selectionPointWidgetComponent = rightSelectionPointWidget;

	// init the tools
	FRGC_UToolInitProperties initProperties;
	initProperties.editor = editorComponent;
	initProperties.flexSimulation = flexComponent;
	initProperties.leftSelectionPoint = leftInteractionPoint;
	initProperties.rightSelectionPoint = rightInteractionPoint;
	initProperties.toolDelegate = this;
	initProperties.targetComponent = regionGrowingRoot;
	initProperties.regionGrowingComponent = regionGrowingComponent;
	initProperties.developerMode = developerMode;

	generativeBrushTool->init(initProperties);

	regionGrowingElementGeneratorTool->exemplarActor = exemplarActor;
	regionGrowingElementGeneratorTool->elementEditor = editorComponent;
	regionGrowingElementGeneratorTool->init(initProperties);

	stringGeneratorTool->exemplarActor = exemplarActor;
	stringGeneratorTool->elementEditor = editorComponent;
	stringGeneratorTool->init(initProperties);

	collagenGeneratorTool->exemplarActor = exemplarActor;
	collagenGeneratorTool->elementEditor = editorComponent;
	collagenGeneratorTool->init(initProperties);

	swarmGeneratorTool->exemplarActor = exemplarActor;
	swarmGeneratorTool->elementEditor = editorComponent;
	swarmGeneratorTool->init(initProperties);
	
	selectionTool->init(initProperties);

	digTool->init(initProperties, camera);

	meshCollectionTool->init(initProperties, camera);

	exemplarInspectorTool->init(initProperties, editorComponent,  camera);

	eventVisualizationTool->init(initProperties, flexSimulation, camera);

	agentPathlineTool->init(initProperties, flexSimulation, camera);

	physicalInteractionTool->init(initProperties, flexSimulation, camera);

	if (_interactionMode == ESketchyInteractionMode::Simulating)
		_currentTool = eventVisualizationTool;
	else 
		_currentTool = regionGrowingElementGeneratorTool;

	_currentTool->gainControl();
}

FString AVRSketchyPawn::_actorLabelByDate(FString baseName)
{
	FDateTime now = FDateTime::Now();

	return baseName += now.ToString();
}

void AVRSketchyPawn::takeGraphicalSnapshotAndHighResShot()
{
	if (!developerMode)
		return;

	if (flexComponent && _interactionMode == ESketchyInteractionMode::Simulating )
	{
		takeHighReshShot();
		snapshotSimulationGraphicsToEditor();
	}
	else if( _interactionMode == ESketchyInteractionMode::Painting)
	{
		takeHighReshShot();
		snapshotElementDomain();
	}
}

void AVRSketchyPawn::takeHighReshShot()
{
	APlayerController * player = UGameplayStatics::GetPlayerController( GetWorld(), 0 );

	if(player)
	{
		player->ConsoleCommand( TEXT( "HighResShot 2" ), true );
	}
}

void AVRSketchyPawn::snapshotElementDomain()
{
#if WITH_EDITOR
	UWorld * world = GEditor->EditorWorld;

	ASimulationSnapshotActor * ismcSnapshot = regionGrowingComponent->createGraphicalSnapshotActor(world);

	FString baseName = regionGrowingComponent->GetOwner()->GetName() + "_elementSnapshot_";

	FString name = _actorLabelByDate(baseName);

	ismcSnapshot->SetActorLabel(name);

	ismcSnapshot->GetRootComponent()->SetVisibility(false, true);
	ismcSnapshot->SetActorEnableCollision(false);

	_simulationSnapShotCount++;
#endif // WITH_EDITOR
}

void AVRSketchyPawn::snapshotSimulationStateToEditor()
{
#if WITH_EDITOR
	if (!flexComponent)
		return;

	UWorld * world = GEditor->EditorWorld;

	AActor * actor = flexComponent->snapshotSimulationStateToWorld(world);

	FString baseName = flexSimulationActor->GetName() + "_simulationStateSnapshot_";

	FString name = _actorLabelByDate(baseName);

	actor->SetActorLabel(name);
#endif
}

void AVRSketchyPawn::snapshotMeshInterfaceRuntimeMeshComponentToActor(AActor* actor)
{
	if (!regionGrowingComponent)
		return;

	auto rmc = regionGrowingComponent->meshInterfaceRuntimeMesh;

	Utility::duplicateRuntimeMeshComponentToActor(rmc, actor);
}

void AVRSketchyPawn::snapshotSimulationGraphicsToEditor()
{
#if WITH_EDITOR
	if (!flexComponent)
		return;

	UWorld * world = GEditor->EditorWorld;

	ASimulationSnapshotActor * ismcSnapshot = flexComponent->createGraphicalSnapshotActor(world);

	ismcSnapshot->GetRootComponent()->SetVisibility(false, true);
	ismcSnapshot->SetActorEnableCollision(false);

	snapshotMeshInterfaceRuntimeMeshComponentToActor(ismcSnapshot);

	FString baseName = flexSimulationActor->GetName() + "_graphicalSnapshot_";

	// name it
	FString name = _actorLabelByDate(baseName);

	ismcSnapshot->SetActorLabel(name);
#endif // WITH_EDITOR
}

// ------------------------------------------------------------
// Left touch
// ------------------------------------------------------------

FVector2D AVRSketchyPawn::_getRightTouchPoint()
{
	FVector2D p;

	p.X = InputComponent->GetAxisKeyValue( EKeys::MotionController_Right_Thumbstick_X );
	p.Y = InputComponent->GetAxisKeyValue( EKeys::MotionController_Right_Thumbstick_Y );

	return p;
}

void AVRSketchyPawn::rightController_touchStart()
{
	_rightTouchActive = true;

	FVector2D p = _getRightTouchPoint();

	_currentTool->rightTouchStart( p );
}

void AVRSketchyPawn::rightController_touchUpdated()
{
	if(!_rightTouchActive)
		return;
	
	FVector2D p = _getRightTouchPoint();

	_currentTool->rightTouchUpdated( p );
}

void AVRSketchyPawn::rightController_touchEnd()
{
	_rightTouchActive = false;

	FVector2D p = _getRightTouchPoint();

	_currentTool->rightTouchEnd();
}

// ------------------------------------------------------------
// Right touch
// ------------------------------------------------------------

FVector2D AVRSketchyPawn::_getLeftTouchPoint()
{
	FVector2D p;

	p.X = InputComponent->GetAxisKeyValue( EKeys::MotionController_Left_Thumbstick_X );
	p.Y = InputComponent->GetAxisKeyValue( EKeys::MotionController_Left_Thumbstick_Y );

	return p;
}

void AVRSketchyPawn::leftController_touchStart()
{
	_leftTouchActive = true;
}

void AVRSketchyPawn::leftController_touchUpdated()
{

}

void AVRSketchyPawn::leftController_touchEnd()
{
	_leftTouchActive = false;
}

// ------------------------------------------------------------

void AVRSketchyPawn::leftController_upFace_released()
{
	if (!developerMode)
		return;

	restoreSimulation();
}

void AVRSketchyPawn::leftController_downFace_released()
{
	if (!developerMode)
		return;

	snapshotSimulation();
}

void AVRSketchyPawn::leftController_leftFace_released()
{

}

void AVRSketchyPawn::leftController_rightFace_released()
{	

}


void AVRSketchyPawn::rightController_faceUp_pressed()
{
	if (!_currentTool)
		return;

	_currentTool->faceUp_pressed();
}

void AVRSketchyPawn::rightController_faceUp_released()
{
	if(!_currentTool)
		return;

	_currentTool->faceUp_released(rightInteractionPoint);
}

void AVRSketchyPawn::rightController_shoulder_released()
{
	if (!_currentTool || !_currentTool->consume_rightShoulder_pressed())
		ShowToolSelectPopup();
}






void AVRSketchyPawn::rightController_faceDown_pressed()
{
	if (!_currentTool)
		return;

	_currentTool->faceDown_pressed();
}

void AVRSketchyPawn::rightController_faceDown_released()
{
	if (!_currentTool)
		return;

	_currentTool->faceDown_released();
}




void AVRSketchyPawn::rightController_faceLeft_pressed()
{
	if (!_currentTool)
		return;

	_currentTool->faceLeft_pressed();
}

void AVRSketchyPawn::rightController_faceLeft_released()
{
	if (!_currentTool)
		return;

	_currentTool->faceLeft_released();
}




void AVRSketchyPawn::rightController_faceRight_pressed()
{
	if(!_currentTool)
		return;

	_currentTool->faceRight_pressed();
}

void AVRSketchyPawn::rightController_faceRight_released()
{
	if (!_currentTool)
		return;

	_currentTool->faceRight_released();
}



void AVRSketchyPawn::_initSimulation_oneTime()
{
	if (_didInitSimulation)
		return;

	if (flexComponent && editorComponent)
	{
		_initSimulationBounds();

		// it's already in world space
		FTransform meshInterfaceToWorld = FTransform::Identity;

		flexComponent->init(editorComponent->context->meshInterface, meshInterfaceToWorld, camera);
	}

	_didInitSimulation = true;
}

void AVRSketchyPawn::_initSimulationBounds()
{
	FTransform transform = flexSimulationActor->GetRootComponent()->GetComponentTransform();

	FBox limits = editorComponent->context->limits;

	limits.Min = transform.InverseTransformPosition(limits.Min);
	limits.Max = transform.InverseTransformPosition(limits.Max);

	flexSimulation->setInstanceManagerBounds(limits);
}

void AVRSketchyPawn::_startSimulation()
{
	if (!flexSimulation || !editorComponent)
		return;

	editorComponent->stop();

	_initSimulation_oneTime();

	auto exportedDomain = editorComponent->exportOutputDomain();
	flexSimulation->loadExportedDomainInfo(exportedDomain);

	flexSimulation->begin();

	// keep it paused
	flexSimulation->pause();
}

void AVRSketchyPawn::_endSimulation()
{
	if(!flexSimulation || !editorComponent)
		return;

	auto elementDomain = flexSimulation->exportElementDomain();

	editorComponent->loadElementDomain(elementDomain);
	editorComponent->start();

	flexSimulation->pause();
	flexSimulation->clear();
}

void AVRSketchyPawn::leftTrigger( float value )
{
	if(!_currentTool)
		return;

	_currentTool->setLeftTriggerValue( value );
	_lastLeftTriggerValue = value;

	if(_leftTriggerDown && value <= 0.1f)
	{
		_leftTriggerDown = false;
		_currentTool->leftEnd();
	}

	if(!_leftTriggerDown && value > 0.1f)
	{
		_leftTriggerDown = true;
		_currentTool->leftStart();
	}
}

void AVRSketchyPawn::rightTrigger( float value )
{
	if(!_currentTool)
		return;

	_currentTool->setRightTriggerValue( value );
	_lastRightTriggerValue = value;

	if(_rightTriggerDown && value <= 0.1f)
	{
		_rightTriggerDown = false;
		_currentTool->rightEnd();
	}

	if(!_rightTriggerDown && value > 0.1f)
	{
		_rightTriggerDown = true;
		_currentTool->rightStart();
	}
}

void AVRSketchyPawn::_startCurrentToolLeftTrigger()
{
	_currentTool->setLeftTriggerValue( _lastLeftTriggerValue );

	if(_leftTriggerDown && _lastLeftTriggerValue > 0.1f)
	{
		_currentTool->leftStart();
	}
}

void AVRSketchyPawn::_startCurrentToolRightTrigger()
{
	_currentTool->setRightTriggerValue( _lastRightTriggerValue );

	if(_rightTriggerDown && _lastRightTriggerValue > 0.1f)
	{
		_currentTool->rightStart();
	}
}

void AVRSketchyPawn::cedeFocus(UTool * tool)
{
	if (tool == digTool)
	{
		_currentTool->loseFocus();
		_currentTool = nullptr;

		ShowToolSelectPopup();
	}
}

void AVRSketchyPawn::spiderManLeftStart()
{
	_wantsLeftSpiderMan = true;

	_lastLeft = leftController->GetComponentLocation();
}

void AVRSketchyPawn::spiderManLeftEnd()
{
	_wantsLeftSpiderMan = false;
}

void AVRSketchyPawn::spiderManLeftUpdate()
{
	FVector left = leftController->GetComponentLocation();

	FVector delta = (left - _lastLeft);

	this->AddActorWorldOffset( -delta, false, nullptr, ETeleportType::TeleportPhysics );

	_lastLeft = leftController->GetComponentLocation();
}

void AVRSketchyPawn::spiderManRightStart()
{
	_wantsRightSpiderMan = true;

	_lastRight = rightController->GetComponentLocation();
}

void AVRSketchyPawn::spiderManRightEnd()
{
	_wantsRightSpiderMan = false;
}

void AVRSketchyPawn::spiderManRightUpdate()
{
	FVector right = rightController->GetComponentLocation();

	FVector delta = (right - _lastRight);

	this->AddActorWorldOffset( -delta, false, nullptr, ETeleportType::TeleportPhysics );

	_lastRight = rightController->GetComponentLocation();
}

void AVRSketchyPawn::setCurrentTool( UTool * newCurrentTool )
{
	if(_currentTool == newCurrentTool)
		return;

	if(_currentTool)
		_currentTool->releaseControl();

	_currentTool = newCurrentTool;

	if(_currentTool)
	{
		_currentTool->gainControl();

		// fire off previous trigger values
		_startCurrentToolLeftTrigger();
		_startCurrentToolRightTrigger();

		if(_rightTouchActive)
			rightController_touchStart();
	}

	if (_currentTool == generativeBrushTool)
	{
		if (generativeBrushTool->automaticDrawMode == EAutomaticDrawMode::ByElementHint && selectionTool->_selection )
		{
			for (auto e : *selectionTool->_selection)
			{
				if (e->spaceModeHint == ESpaceMode::Surface)
					generativeBrushTool->setDrawMode(EGenerativeDrawMode::Surface);
				else if (e->spaceModeHint == ESpaceMode::Volume)
					generativeBrushTool->setDrawMode(EGenerativeDrawMode::Volumetric);
			}
		}
	}
}

void AVRSketchyPawn::setInteractionMode(ESketchyInteractionMode newMode)
{
	if (newMode == _interactionMode)
		return;

	if (_interactionMode == ESketchyInteractionMode::Simulating && newMode == ESketchyInteractionMode::Painting)
	{
		_endSimulation();

		if (editorComponent)
			editorComponent->start();
	}
	else if (_interactionMode == ESketchyInteractionMode::Painting && newMode == ESketchyInteractionMode::Simulating)
	{
		if (editorComponent)
			editorComponent->stop();

		_startSimulation();
	}

	_interactionMode = newMode;
}

bool AVRSketchyPawn::isSimulating()
{
	if (!flexSimulation)
		return false;

	return flexSimulation->isPlaying();
}

void AVRSketchyPawn::setSimulating( bool simulating )
{
	if (!flexSimulation)
		return;

	if (_interactionMode != ESketchyInteractionMode::Simulating)
		return;

	if (simulating)
	{
		_initSimulation_oneTime();
		flexSimulation->play();
	}
	else
		flexSimulation->pause();
}

void AVRSketchyPawn::snapshotSimulation()
{
	if (!flexSimulation)
		return;

	_snapshot.snapshot(flexSimulation->graphSimulation);
}

void AVRSketchyPawn::restoreSimulation()
{
	if (!flexSimulation)
		return;

	_snapshot.restore(flexSimulation->graphSimulation);
}

void AVRSketchyPawn::ShowToolSelectPopup_Implementation()
{

}
