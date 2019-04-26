// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Toolkits/ToolkitManager.h"
#include "IStructureDetailsView.h"

#include "ElementEditor/DiscreteElementEditorComponent.h"
#include "Tools/MeshCollectionTool.h"

#include "Kismet/KismetMathLibrary.h"

#include "ExemplarInspectorTool.h"

UExemplarInspectorTool::~UExemplarInspectorTool()
{

}

void UExemplarInspectorTool::init(FUToolInitProperties& initProperties, UDiscreteElementEditorComponent * elementEditor, UCameraComponent * camera)
{
	UTool::init(initProperties);

	_elementEditor = elementEditor;
	_camera = camera;
}

void UExemplarInspectorTool::gainFocus()
{
	_selection = nullptr;

	_initGrabTool();
	_spawnInspectorActor();
	_spawnMeshCollectionActor();
}

void UExemplarInspectorTool::loseFocus()
{
	// nuke the inspector
	if (inspectorActor) 
	{
		inspectorActor->Destroy();
		inspectorActor = nullptr;
	}

	if (_meshCollectionActor)
	{
		_meshCollectionActor->Destroy();
		_meshCollectionActor = nullptr;
	}

	if (_selection)
	{
		_selection->hideSelectionOutline();
		_selection = nullptr;
	}
}

float UExemplarInspectorTool::_brushRadius()
{
	return 10.0f * selectionATriggerValue();
}



auto UExemplarInspectorTool::nearestElementActor(UPrimitiveComponent * hand) -> std::pair<AElementActor*, float>
{
	// perform an overlap test at the hand position to find actors that overlap
	TArray<FOverlapResult> overlaps;
	FCollisionShape collisionShape = FCollisionShape::MakeSphere(overlapRadius);

	FCollisionQueryParams params;

	UActorComponent * root;
	if (_elementEditor)
		root = _elementEditor;
	else
		root = _elementEditor;

	params.AddIgnoredActor(root->GetOwner());
	params.AddIgnoredActor(hand->GetOwner());

	bool overlap = root->GetWorld()->OverlapMultiByChannel(
		overlaps,
		hand->GetComponentLocation(),
		hand->GetComponentRotation().Quaternion(),
		ECollisionChannel::ECC_WorldDynamic,
		collisionShape,
		params
	);

	float minDistance = std::numeric_limits<float>::max();
	AElementActor * minActor = nullptr;

	// update the selection
	for (FOverlapResult overlap : overlaps)
	{
		AElementActor * elementActor = Cast<AElementActor>(overlap.Actor.Get());

		if (!elementActor) continue;

		elementActor->hideSelectionOutline();

		if (!elementActor->GetRootComponent()->IsVisible())
			continue;

		FVector hitPoint = FVector::ZeroVector;
		float distance = elementActor->ActorGetDistanceToCollision(hand->GetComponentLocation(), ECollisionChannel::ECC_WorldDynamic, hitPoint);

		if (distance < minDistance)
		{
			minDistance = distance;
			minActor = elementActor;
		}
	}

	return std::make_pair(minActor, minDistance);
}

void UExemplarInspectorTool::_initGrabTool()
{
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
}

void UExemplarInspectorTool::_select(AElementActor * actor)
{
	_selection = actor;

	_selection->showSelectionOutline();

	// init the grab tool
	_grabTool->selectActor(_selection, _selectionA);

	warmupOtherTool(_grabTool);
	_grabTool->oneHandStart(_selectionA);

	// show the UI
	inspectorActor->setSelection(_selection);
}

void UExemplarInspectorTool::didGrabItem(UCollectionSpace * collectionSpace, int32 itemAtIndex, FTransform grabTransform, UPrimitiveComponent * grabbedCell, FTransform cellTransform, FBox cellBounds)
{
	// create a new element actor based on the mesh
	UWorld * world = _elementEditor->GetWorld();

	AElementActor * element = world->SpawnActorAbsolute<AElementActor>(AElementActor::StaticClass(), cellTransform);

	element->SetMobility(EComponentMobility::Movable);

	// eww, hacky
	element->overrideRadius = 1.0f;
	element->scaleOverrideRadius = false;
	element->generationParameters.radius = 12.0f;
	element->generationParameters.kNearest = 20;
	element->optimizationParameters.radius = 0.0f;
	element->minAssignmentDistance = 0.5f;
	element->freespaceRadius = 0.5f;
	element->generationInnerRadius = 10.0f;
	element->particleOrMember = EAggregateProxyParticleOrMember::Particle;

	UStaticMeshComponent * elementMesh = element->GetStaticMeshComponent();
	UStaticMeshComponent * draggedMesh = Cast<UStaticMeshComponent>(grabbedCell);

	if (!(elementMesh && draggedMesh)) return;

	elementMesh->SetStaticMesh(draggedMesh->GetStaticMesh());

	int32 n = draggedMesh->GetNumMaterials();
	for (int i = 0; i < n; ++i)
	{
		elementMesh->SetMaterial(i, draggedMesh->GetMaterial(i));
	}

	// add it to the exemplar
	AActor * exemplar = _elementEditor->exemplarActor();

	if (!exemplar) return;

	element->AttachToActor(exemplar, FAttachmentTransformRules::KeepWorldTransform);

	mode = Mode::Selection;

	_select(element);
}

void UExemplarInspectorTool::_spawnMeshCollectionActor()
{
	if (!_camera) return;
	if (!meshCollectionActorClass) return;

	if (!_meshCollectionActor)
	{
		UWorld * world = _camera->GetWorld();

		FTransform inspectorTransform = _inspectorSpawnTransform();

		// move down 100 units in the inspector space
		FVector offset(0.0f, 25.0f, 0.0f);
		offset = inspectorTransform.TransformPosition(offset);

		FQuat rotation = FQuat(FVector::ForwardVector, M_PI * 0.5f);
		rotation = inspectorTransform.GetRotation() * rotation;

		inspectorTransform.SetLocation(offset);
		inspectorTransform.SetRotation(rotation);

		_meshCollectionActor = world->SpawnActor<AMeshCollectionSpaceActor>(meshCollectionActorClass, inspectorTransform);
		_meshCollectionActor->delegate = this;
	}
}

FTransform UExemplarInspectorTool::_inspectorSpawnTransform()
{
	const FVector cameraLocation = _camera->GetComponentLocation();
	const FVector handLocation = _rightSelectionPoint->GetComponentLocation();

	FQuat rotation = UKismetMathLibrary::FindLookAtRotation(handLocation, cameraLocation).Quaternion();

	FTransform transform(rotation, handLocation, FVector(1.0f));

	return transform;
}

void UExemplarInspectorTool::_spawnInspectorActor()
{
	if (!_camera) return;
	if (!inspectorActorClass) return;

	UWorld * world = _camera->GetWorld();

	FTransform transform = _inspectorSpawnTransform();

	inspectorActor = world->SpawnActor<AExemplarInspectorActor>(inspectorActorClass, transform);
}


void UExemplarInspectorTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastToWorldTransform)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tickOneHand(dt, hand, lastToWorldTransform);
	}
	else if (mode == Mode::MeshCollection)
	{
		_meshCollectionActor->collectionSpace->update_oneHand(dt, hand, lastToWorldTransform);
	}
}

void UExemplarInspectorTool::oneHandStart(UPrimitiveComponent * hand)
{
	auto minActorPair = nearestElementActor(hand);

	AElementActor * minActor = minActorPair.first;

	FVector collectionNearest = _meshCollectionActor->collectionSpace->nearest(hand);

	float nearestElementDistance = minActorPair.second;
	float nearestCollectionDistance = FVector::Dist(collectionNearest, hand->GetComponentLocation());

	_selection = nullptr;

	mode = Mode::None;

	if (nearestCollectionDistance < nearestElementDistance)
	{
		mode = Mode::MeshCollection;

		_meshCollectionActor->collectionSpace->begin_oneHand(hand);
	}
	else if (nearestElementDistance < nearestCollectionDistance && minActor)
	{
		mode = Mode::Selection;

		_select(minActor);
	}
}

void UExemplarInspectorTool::oneHandEnd(UPrimitiveComponent * hand)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->oneHandEnd(hand);
	}
	else if (mode == Mode::MeshCollection)
	{
		_meshCollectionActor->collectionSpace->end_oneHand(hand);
	}
}

void UExemplarInspectorTool::twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->twoHandStart(handA, handB);
	}
}

void UExemplarInspectorTool::twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->twoHandEnd(handA, handB);
	}
}

void UExemplarInspectorTool::tickTwoHand(float dt, UPrimitiveComponent * handA, UPrimitiveComponent * handB, FTransform lastTransformA, FTransform lastTransformB)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tickTwoHand(dt, handA, handB, lastTransformA, lastTransformB);
	}
}

void UExemplarInspectorTool::faceDown_released()
{

}

void UExemplarInspectorTool::faceDown_pressed()
{
	if (mode == Mode::Selection && _selection)
	{
		FActorSpawnParameters params;

		params.Template = _selection;

		FTransform transform = _selection->GetTransform();

		// shift to the right a bit
		FBox bounds = _selection->GetComponentsBoundingBox();
		FVector right = _camera->GetComponentTransform().GetRotation().RotateVector(FVector::RightVector);
		transform.AddToTranslation(right * bounds.GetSize().GetMax() * 0.1f);

		// spawn it
		UWorld * world = _selection->GetWorld();

		AElementActor * element = world->SpawnActorAbsolute<AElementActor>(AElementActor::StaticClass(), transform, params);

		element->SetActorTransform(transform, false, nullptr, ETeleportType::TeleportPhysics);

		_select(element);

	}
}

void UExemplarInspectorTool::tick(float dt)
{
	if (mode == Mode::Selection)
	{
		warmupOtherTool(_grabTool);
		_grabTool->tick(dt);
	}
	else if (mode == Mode::DuplicateDrag)
	{
		// emulate a tickOneHand
		warmupOtherTool(_grabTool);
		_grabTool->tickOneHand(dt, _selectionA, _lastA);
	}
}

void UExemplarInspectorTool::didCancel_Implementation(UGrabTool * grabTool, AActor * draggingActor)
{
}

void UExemplarInspectorTool::didPlace_Implementation(UGrabTool * grabTool, AActor * draggingActor)
{
}




AExemplarInspectorActor::AExemplarInspectorActor()
{
	RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("Root"));

	widgetComponent = CreateDefaultSubobject<UWidgetComponent>(TEXT("WidgetComponent"));

	widgetComponent->AttachToComponent(RootComponent, FAttachmentTransformRules::KeepRelativeTransform);
}

void AExemplarInspectorActor::setSelection(UObject * object)
{
	_selection = object;

	FPropertyEditorModule& propertyModule = FModuleManager::GetModuleChecked<FPropertyEditorModule>("PropertyEditor");

	FDetailsViewArgs detailArgs;

	auto detailView = propertyModule.CreateDetailView(detailArgs);
	
	widgetComponent->SetSlateWidget(detailView);

	detailView->SetObject(object, true);







	/*FStructureDetailsViewArgs structDetailsArgs;

	auto graph = &tool->_elementEditor->graph;

	TSharedRef<FStructOnScope> structData(new FStructOnScope(FGraph::StaticStruct(), (uint8*)graph));

	auto structView = propertyModule.CreateStructureDetailView(detailArgs, structDetailsArgs, structData);

	_structDetailsView = structView;

	TSharedPtr<SWidget> slatePtr = structView->GetWidget();

	widgetComponent->SetSlateWidget(slatePtr);*/
}
