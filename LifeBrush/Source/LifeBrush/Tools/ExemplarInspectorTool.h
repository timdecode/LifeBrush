// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "VRTool.h"

#include "RegionGrowingToolInit.h"
#include "GrabTool.h"

#include "WidgetComponent.h"
#include "IStructureDetailsView.h"
#include "MeshCollectionTool.h"

#include "ExemplarInspectorTool.generated.h"

class UExemplarInspectorTool;
class UDiscreteElementEditorComponent;
class AMeshCollectionSpaceActor;
class UMeshCollectionTool;

UCLASS(DefaultToInstanced)
class LIFEBRUSH_API AExemplarInspectorActor : public AActor
{
	GENERATED_BODY()

public:
	AExemplarInspectorActor();

	void setSelection(UObject * object);

public:
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Instanced, Category = "LifeBrush")
	UWidgetComponent * widgetComponent;

protected:
	UPROPERTY()
	UObject * _selection;

	TSharedPtr<class IStructureDetailsView> _structDetailsView;

};



UCLASS( Blueprintable )
class UExemplarInspectorTool : public UTool, public IGrabDelegate, public MeshCollectionSpaceActorDelegate
{
	GENERATED_BODY()

public:
	virtual ~UExemplarInspectorTool();

	void init(FUToolInitProperties& initProperties, UDiscreteElementEditorComponent * elementEditor, UCameraComponent * camera);

	// UTool overrids
	virtual void focused() override;
	virtual void loseFocus() override;

	virtual void tick(float dt) override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;

	auto nearestElementActor(UPrimitiveComponent * hand) -> std::pair<AElementActor*, float>;

	virtual void oneHandEnd(UPrimitiveComponent * hand) override;
	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform) override;

	virtual void twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override;
	virtual void twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override;
	virtual void tickTwoHand(float dt, UPrimitiveComponent * handA, UPrimitiveComponent * handB, FTransform lastTransformA, FTransform lastTransformB) override;


	virtual void faceDown_released() override;


	virtual void faceDown_pressed() override;

public: // IGrabDelegate overrides
	virtual void didCancel_Implementation(UGrabTool * grabTool, AActor * draggingActor) override;
	virtual void didPlace_Implementation(UGrabTool * grabTool, AActor * draggingActor) override;

public:
	virtual void didGrabItem(UCollectionSpace * collectionSpace, int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds) override;

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TSubclassOf<class AExemplarInspectorActor> inspectorActorClass;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	TSubclassOf<class AMeshCollectionSpaceActor> meshCollectionActorClass;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float overlapRadius = 10.0f;

public:
	UPROPERTY()
	AExemplarInspectorActor * inspectorActor;

	UCameraComponent * _camera;
	UDiscreteElementEditorComponent * _elementEditor;

	AElementActor * _selection = nullptr;

	UPROPERTY()
	UGrabTool * _grabTool = nullptr;

	UPROPERTY()
	AMeshCollectionSpaceActor * _meshCollectionActor = nullptr;

	enum class Mode {
		MeshCollection,
		Selection,
		DuplicateDrag,
		None
	};

	Mode mode = Mode::None;

protected:
	void _spawnMeshCollectionActor();
	void _spawnInspectorActor();
	float _brushRadius();
	void _initGrabTool();

	void _select(AElementActor * actor);

	FTransform _inspectorSpawnTransform();
};

