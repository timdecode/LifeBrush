// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "WidgetComponent.h"
#include "VRInterface/CollectionSpace.h"
#include "RegionGrowingComponent.h"
#include "ElementEditor/DiscreteElementEditorComponent.h"
#include "GrabTool.h"
#include "ElementActor.h"
#include "Kismet/KismetMathLibrary.h"

#include "MeshCollectionTool.h"



AMeshCollectionSpaceActor::AMeshCollectionSpaceActor()
{
	RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));
	RootComponent->SetMobility(EComponentMobility::Movable);

	collectionSpace = CreateDefaultSubobject<UCollectionSpace>(TEXT("CollectionSpace"));

	collectionSpace->AttachToComponent(RootComponent, FAttachmentTransformRules::KeepRelativeTransform);

	PrimaryActorTick.bCanEverTick = false;
}

bool AMeshCollectionSpaceActor::ShouldTickIfViewportsOnly() const
{
	return false;
}

void AMeshCollectionSpaceActor::Tick(float DeltaSeconds)
{
	if (!_didBegingPlayInEditor)
	{
		collectionSpace->dataSource = this;
		collectionSpace->delegate = this;

		collectionSpace->reloadData();

		_didBegingPlayInEditor = true;
	}
}

int32 AMeshCollectionSpaceActor::count_Implementation()
{
	return staticMeshes.Num();
}

UPrimitiveComponent * AMeshCollectionSpaceActor::primitiveCellAt_Implementation(int32 index)
{
	UStaticMesh * mesh = staticMeshes[index];

	UStaticMeshComponent * cell = NewObject<UStaticMeshComponent>(collectionSpace);

	cell->AttachToComponent(collectionSpace, FAttachmentTransformRules::KeepRelativeTransform);
	cell->SetStaticMesh(mesh);
	cell->SetMobility(EComponentMobility::Movable);
	cell->SetCollisionProfileName(TEXT("NoCollision"));

	for (int mi = 0; mi < mesh->GetNumSections(0); ++mi)
	{
		cell->SetMaterial(mi, mesh->GetMaterial(mi));
	}

	cell->RegisterComponent();

	return cell;
}

void AMeshCollectionSpaceActor::didGrab_Implementation(int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds)
{
	if (!delegate)
		return;

	delegate->didGrabItem(collectionSpace, itemAtIndex, grabTransform, grabbedCell, cellTransform, cellBounds);
}

void AMeshCollectionSpaceActor::BeginPlay()
{
	Super::BeginPlay();

	if (!_didBegingPlayInEditor)
	{
		collectionSpace->dataSource = this;
		collectionSpace->delegate = this;

		collectionSpace->reloadData();

		_didBegingPlayInEditor = true;
	}
}


void UMeshCollectionTool::init(FRGC_UToolInitProperties& initProperties, UCameraComponent * camera)
{
	UTool::init(initProperties);

	_elementEditor = initProperties.editor;
	_camera = camera;
}

void UMeshCollectionTool::gainFocus()
{
	if (_collectionSpaceActor || !meshCollectionClass || !_elementEditor)
		return;

	spawnCollectionSpaceActor();
}

void UMeshCollectionTool::spawnCollectionSpaceActor()
{
	if (_collectionSpaceActor || !meshCollectionClass || !_camera)
		return;

	UWorld * world = _camera->GetWorld();

	const FVector cameraLocation = _camera->GetComponentLocation();
	const FVector handLocation = _rightSelectionPoint->GetComponentLocation();

	FRotator rotator = UKismetMathLibrary::FindLookAtRotation(handLocation, cameraLocation);

	FTransform transform(rotator.Quaternion(), handLocation, FVector(1.0f));

	_collectionSpaceActor = world->SpawnActor<AMeshCollectionSpaceActor>(meshCollectionClass, transform);
	_collectionSpaceActor->delegate = this;

	_collectionSpace = _collectionSpaceActor->FindComponentByClass<UCollectionSpace>();
}

void UMeshCollectionTool::loseFocus()
{
}

void UMeshCollectionTool::oneHandStart(UPrimitiveComponent * hand)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->oneHandStart(hand);
	
		return;
	}


	if (_collectionSpace)
		_collectionSpace->begin_oneHand(hand);
}

void UMeshCollectionTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->oneHandEnd(hand);

		return;
	}

	if (_collectionSpace)
		_collectionSpace->end_oneHand(hand);
}

void UMeshCollectionTool::twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->twoHandStart(handA, handB);

		return;
	}
}

void UMeshCollectionTool::twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->twoHandEnd(handA, handB);

		return;
	}
}

void UMeshCollectionTool::tick(float dt)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tick(dt);

		return;
	}
}

void UMeshCollectionTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tickOneHand(dt, hand, lastTransform);

		return;
	}

	if (_collectionSpace)
		_collectionSpace->update_oneHand(dt, hand, lastTransform);
}

void UMeshCollectionTool::tickTwoHand(float dt, UPrimitiveComponent * handA, UPrimitiveComponent * handB, FTransform transformA, FTransform transformB)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tickTwoHand(dt, handA, handB, transformA, transformB);

		return;
	}
}

void UMeshCollectionTool::faceDown_released()
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->faceDown_released();

		return;
	}
}

void UMeshCollectionTool::faceUp_released(USceneComponent * interactionPoint /* = nullptr */)
{
	if (_grabTool)
	{
		warmupOtherTool(_grabTool);
		_grabTool->faceUp_released(interactionPoint);

		return;
	}
}


void UMeshCollectionTool::didGrabItem(UCollectionSpace * collectionSpace, int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds)
{
	if (_draggingCell)
	{
		_draggingCell->DestroyComponent();
		_draggingCell = nullptr;
	}

	if (_draggingActor)
	{
		_draggingActor->Destroy();
		_draggingActor = nullptr;
	}

	UStaticMeshComponent * asMesh = Cast<UStaticMeshComponent>(grabbedCell);

	UWorld * world = _elementEditor->GetWorld();


	FVector location = _selectionA->GetComponentLocation();
	FQuat rotation = FQuat::Identity;

	_draggingActor = world->SpawnActorAbsolute<AActor>(AActor::StaticClass(), cellTransform);

	_draggingCell = NewObject<UStaticMeshComponent>(_draggingActor);

	_draggingCell->SetStaticMesh(asMesh->GetStaticMesh());
	for (int i = 0; i < asMesh->GetNumMaterials(); ++i)
	{
		_draggingCell->SetMaterial(i, asMesh->GetMaterial(i));
	}

	_draggingCell->AttachToComponent(_draggingActor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
	_draggingCell->SetMobility(EComponentMobility::Movable);

	_draggingActor->SetRootComponent(_draggingCell);

	_draggingCell->RegisterComponent();

	_draggingActor->SetActorTransform(cellTransform, false, nullptr, ETeleportType::TeleportPhysics);

	if (!_grabTool)
	{
		_grabTool = NewObject<UGrabTool>(this);

		FUToolInitProperties initProperties;
		initProperties.leftSelectionPoint = _leftSelectionPoint;
		initProperties.rightSelectionPoint = _rightSelectionPoint;
		initProperties.targetComponent = targetComponent;

		_grabTool->init(initProperties);

		_grabTool->delegate = this;
	}

	_grabTool->snapController.snapRotation = true;
	_grabTool->selectActor(_draggingActor, _selectionA);

	warmupOtherTool(_grabTool);
}


// -----------------------------------------------
// IGrabDelegate
// -----------------------------------------------

void UMeshCollectionTool::didPlace_Implementation(UGrabTool * grabTool, AActor * draggingActor)
{
	_grabTool = nullptr; // supposedly the GC will eat it later

	_draggingActor = nullptr;
	_draggingCell = nullptr;

	AActor * exemplar = _elementEditor->exemplarActor();

	if (!exemplar)
	{
		draggingActor->Destroy();
		return;
	}

	// create a new element actor based on the mesh
	UWorld * world = _elementEditor->GetWorld();

	FTransform transform = draggingActor->GetTransform();

	AElementActor * element = world->SpawnActorAbsolute<AElementActor>(AElementActor::StaticClass(), transform);

	UStaticMeshComponent * elementMesh = element->GetStaticMeshComponent();
	UStaticMeshComponent * draggedMesh = draggingActor->FindComponentByClass<UStaticMeshComponent>();

	if (!(elementMesh && draggedMesh))
		return;

	elementMesh->SetStaticMesh(draggedMesh->GetStaticMesh());

	int32 n = draggedMesh->GetNumMaterials();
	for (int i = 0; i < n; ++i)
	{
		elementMesh->SetMaterial(i, draggedMesh->GetMaterial(i));
	}

	// add it to the exemplar
	element->AttachToActor(exemplar, FAttachmentTransformRules::KeepWorldTransform);

	// we don't need the draggingActor anymore
	draggingActor->Destroy();
}

void UMeshCollectionTool::didCancel_Implementation(UGrabTool * grabTool, AActor * draggingActor)
{
	_grabTool = nullptr; // supposedly the GC will eat it later

	_draggingActor = nullptr;
	_draggingCell = nullptr;

	draggingActor->Destroy();
}
