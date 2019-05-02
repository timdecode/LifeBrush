//
//  VRTool.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2018-03-06, based on my code for Ship Editor.
//  Copyright (c) Timothy Davison. All rights reserved.
//

#pragma once

#include "GameFramework/Pawn.h"
#include "VRTool.generated.h"

class UMotionControllerComponent;
class URegionGrowingComponent;
class UWidgetComponent;
class UDiscreteElementEditorComponent;
class UFlexSimulationComponent;

class UToolDelegate
{
public:
	virtual void cedeFocus(class UTool * tool) {};
};

struct FUToolInitProperties
{
public:
	UDiscreteElementEditorComponent * editor;
	UFlexSimulationComponent * flexSimulation;
	UPrimitiveComponent * leftSelectionPoint;
	UPrimitiveComponent * rightSelectionPoint;
	USceneComponent * targetComponent;
	bool developerMode;
	UToolDelegate * toolDelegate;
};



/**
*
*/
UCLASS( BlueprintType )
class UTool : public UObject
{
	GENERATED_BODY()

public:
	void init
	( 
		FUToolInitProperties& initProperties
	)
	{
		developerMode = initProperties.developerMode;

		targetComponent = initProperties.targetComponent;

		toolDelegate = initProperties.toolDelegate;

		_leftSelectionPoint = initProperties.leftSelectionPoint;
		_rightSelectionPoint = initProperties.rightSelectionPoint;

		_selectionA = initProperties.leftSelectionPoint;
		_selectionB = initProperties.rightSelectionPoint;

		_handMode = HandMode::None;


		_updateLast();
	}

	virtual void tick(float dt) {}

	virtual void oneHandStart( UPrimitiveComponent * hand ) {}
	virtual void oneHandEnd( UPrimitiveComponent * hand ) {}

	virtual void twoHandStart( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) {}
	virtual void twoHandEnd( UPrimitiveComponent * handA, UPrimitiveComponent * handB ) {}

	// Called when the tool is released. Any action should be aborted.
	virtual void loseFocus();

	// Called when the tool gains control.
	virtual void gainFocus() {}

	virtual void tickOneHand( float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform ) {}

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform lastTransformA,
		FTransform lastTransformB
	) {}

	// Face Down
	// -----------------------------
	virtual void faceDown_released()
	{

	}

	virtual void faceDown_touchStart()
	{

	}

	virtual void faceDown_touchEnd()
	{

	}

	virtual void faceDown_pressed()
	{

	}

	// Face Up
	// -----------------------------
	virtual void faceUp_released(USceneComponent * interactionPoint = nullptr)
	{

	}

	virtual void faceUp_touchStart()
	{

	}

	virtual void faceUp_touchEnd()
	{

	}

	virtual void faceUp_pressed()
	{

	}

	// Face Right
	// -----------------------------
	virtual void faceRight_released()
	{

	}

	virtual void faceRight_pressed()
	{

	}

	// Face Left
	// -----------------------------
	virtual void faceLeft_released()
	{

	}

	virtual void faceLeft_pressed()
	{

	}

	// Shoulder Right
	// -----------------------------
	// The right shoulder button was pressed. If the event is not consumed, other tools can process the event.
	// \return Whether the shoulder button event pressed event is consumed. True, for consumed, false if it wasn't consumed by this tool.
	virtual bool consume_rightShoulder_pressed()
	{
		return false;
	}
	
	virtual void warmupOtherTool(UTool * other)
	{
		other->_selectionA = _selectionA;
		other->_selectionB = _selectionB;

		other->_lastA = _lastA;
		other->_lastB = _lastB;

		other->_triggerA = _triggerA;
		other->_triggerB = _triggerB;

		other->_rightTouchActive = _rightTouchActive;

		other->_rightTouchStartDirection = _rightTouchStartDirection;
	}

	void doTick( float dt )
	{
		if(_handMode == HandMode::OneHand)
			tickOneHand( dt, _selectionA, _lastA );
		else if(_handMode == HandMode::TwoHand)
			tickTwoHand( dt, _selectionA, _selectionB, _lastA, _lastB );

		tick(dt);

		_updateLast();
	}

	void rightStart()
	{
		if(_handMode == HandMode::None)
		{
			_selectionA = _rightSelectionPoint;
			_selectionB = _leftSelectionPoint;

			_updateLast();

			_handMode = HandMode::OneHand;

			oneHandStart( _selectionA );
		}
		else if(_handMode == HandMode::OneHand)
			_twoHandStart();
	}

	void leftStart()
	{
		if(_handMode == HandMode::None)
		{

			_selectionA = _leftSelectionPoint;
			_selectionB = _rightSelectionPoint;

			_updateLast();

			_handMode = HandMode::OneHand;

			oneHandStart( _selectionA );
		}
		else if(_handMode == HandMode::OneHand)
			_twoHandStart();
	}

	void rightEnd()
	{
		if(_handMode == HandMode::OneHand)
		{
			_handMode = HandMode::None;

			oneHandEnd( _selectionA );
		}
		else if(_handMode == HandMode::TwoHand)
		{
			_handMode = HandMode::OneHand;

			_selectionA = _leftSelectionPoint;
			_selectionB = _rightSelectionPoint;

			twoHandEnd( _selectionA, _selectionB );
		}

		_updateLast();
	}

	void leftEnd()
	{
		if(_handMode == HandMode::OneHand)
		{
			_handMode = HandMode::None;

			oneHandEnd( _selectionA );
		}
		else if(_handMode == HandMode::TwoHand)
		{
			_handMode = HandMode::OneHand;

			_selectionA = _rightSelectionPoint;
			_selectionB = _leftSelectionPoint;

			twoHandEnd( _selectionA, _selectionB );
		}

		_updateLast();
	}

	void setLeftTriggerValue( float value )
	{
		if(_selectionA == _leftSelectionPoint)
			_triggerA = value;
		else if(_selectionB == _leftSelectionPoint)
			_triggerB = value;

		if(_handMode == HandMode::OneHand)
			_triggerB = 0.0f;
		else if(_handMode == HandMode::None)
		{
			_triggerA = 0.0f;
			_triggerB = 0.0f;
		}
	}

	void setRightTriggerValue( float value )
	{
		if(_selectionA == _rightSelectionPoint)
			_triggerA = value;
		else if(_selectionB == _rightSelectionPoint)
			_triggerB = value;

		if(_handMode == HandMode::OneHand)
			_triggerB = 0.0f;
		else if(_handMode == HandMode::None)
		{
			_triggerA = 0.0f;
			_triggerB = 0.0f;
		}
	}

	void rightTouchStart( FVector2D p )
	{
		_rightTouchActive = true;

		if( p.Y >= 0.0f )
		{
			_rightTouchStartDirection = TouchDirection::Down;
			faceDown_touchStart();
		}
		else
		{
			_rightTouchStartDirection = TouchDirection::Up;
			faceUp_touchStart();
		}
	}

	void rightTouchUpdated( FVector2D p )
	{

	}

	void rightTouchEnd()
	{
		if( _rightTouchStartDirection == TouchDirection::Down )
			faceDown_touchEnd();
		else if( _rightTouchStartDirection == TouchDirection::Up )
			faceUp_touchEnd();

		_rightTouchActive = false;
	}

	float selectionATriggerValue()
	{
		return _triggerA;
	}

	float selectionBTriggerValue()
	{
		return _triggerB;
	}

	void gainControl()
	{
		_loadWidgets();
		
		gainFocus();
	}

	void releaseControl()
	{
		loseFocus();

		_hideWidgets();

		_handMode = HandMode::None;

		_selectionA = _leftSelectionPoint;
		_selectionB = _rightSelectionPoint;

		_lastA = FTransform::Identity;
		_lastB = FTransform::Identity;

		_rightTouchActive = false;
	}

	// Call this if both triggers should be considered released. It will trigger a oneHandEnd, or a twoHandEnd followed by a oneHandEnd.
	void endAction()
	{
		if(_handMode == HandMode::OneHand)
		{
			_handMode = HandMode::None;

			oneHandEnd( _selectionA );
		}
		else if(_handMode == HandMode::TwoHand)
		{
			_handMode = HandMode::None;

			twoHandEnd( _selectionA, _selectionB );
			oneHandEnd( _selectionA );

			_selectionA = _leftSelectionPoint;
			_selectionB = _rightSelectionPoint;
		}
	}

public:
	FTransform targetLocalTransform(FTransform toWorld)
	{
		FTransform inverseGraph = targetComponent->GetOwner()->GetRootComponent()->GetComponentTransform().Inverse();

		return toWorld * inverseGraph;
	}

protected:
	void _twoHandStart()
	{
		_updateLast();

		_handMode = HandMode::TwoHand;

		twoHandStart( _selectionA, _selectionB );
	}

	void _updateLast()
	{
		_lastA = _selectionA->GetComponentTransform();
		_lastB = _selectionB->GetComponentTransform();
	}



	void _hideWidgets();

	void _loadWidgets();

public:
	USceneComponent * targetComponent = nullptr;

	UWidgetComponent * widgetComponent = nullptr;
	UWidgetComponent * selectionPointWidgetComponent = nullptr;

	UToolDelegate * toolDelegate = nullptr;

	bool developerMode = true;

	virtual TSubclassOf<class UUserWidget> getWidgetClass() { return nullptr; }
	virtual TSubclassOf<class UUserWidget> getSelectionWidgetClass() { return nullptr; }

	enum class HandMode
	{
		None,
		OneHand,
		TwoHand
	};

protected:

	HandMode _handMode = HandMode::None;

	UPrimitiveComponent * _leftSelectionPoint;
	UPrimitiveComponent * _rightSelectionPoint;

	UPrimitiveComponent * _selectionA;
	UPrimitiveComponent * _selectionB;

	FTransform _lastA;
	FTransform _lastB;

	float _triggerA = 0.0f;
	float _triggerB = 0.0f;

	bool _rightTouchActive = false;
	enum class TouchDirection
	{
		Up,
		Down
	};

	TouchDirection _rightTouchStartDirection;
};

/**
* A tool that draws a brush sphere controlled by the trigger.
*/
UCLASS(BlueprintType)
class UBrushTool : public UTool
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * brushMesh;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * brushMeshMaterial;

	// The brushMesh will be scaled by the trigger value and the scale factor. The size of the mesh component should
	// match the brush radius.
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float brushMeshScaleFactor = 2.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float brushMinRadius = 2.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	float brushMaxRadius = 6.0f;


public:
	~UBrushTool();

	virtual void loseFocus() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform) override;

	virtual bool shouldShowBrush() { return true; }

protected:
	void _createBrushMeshComponent(UPrimitiveComponent * selectionPoint);
	float _brushRadius();

protected:
	UStaticMeshComponent * _brushMeshComponent = nullptr;
};