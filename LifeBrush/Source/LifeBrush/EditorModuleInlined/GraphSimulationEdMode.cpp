// Copyright 2019 Timothy Davison, Inc. All Rights Reserved.
// Based on my code for Ship Editor.

#include "LifeBrush.h"

#if WITH_EDITOR

#include "GraphSimulationEdMode.h"
#include "Toolkits/ToolkitManager.h"

#include "Simulation/FlexElements.h"
#include "ShipEditorSimulation/MeshSimulation.h"

#include "IStructureDetailsView.h"
#include "Simulation/Aggregates.h"

#define LOCTEXT_NAMESPACE "FGraphSimulationEdMode"

const FEditorModeID FGraphSimulationEdMode::EM_GraphSimulationEdModeId = TEXT("EM_GraphSimulationEdMode");

FGraphSimulationEdMode::FGraphSimulationEdMode()
{

}

FGraphSimulationEdMode::~FGraphSimulationEdMode()
{

}

void FGraphSimulationEdMode::Enter()
{
	FEdMode::Enter();

	if (!Toolkit.IsValid())
	{
		Toolkit = MakeShareable(new FGraphSimulationToolKit(this));
		Toolkit->Init(Owner->GetToolkitHost());

		graphToolKit = StaticCastSharedPtr<FGraphSimulationToolKit>(Toolkit);
	}
}

void FGraphSimulationEdMode::Exit()
{
	if (Toolkit.IsValid())
	{
		FToolkitManager::Get().CloseToolkit(Toolkit.ToSharedRef());
		Toolkit.Reset();
	}

	// Call base Exit method to ensure proper cleanup
	FEdMode::Exit();
}

/** FEdMode: Called when a key is pressed */
bool FGraphSimulationEdMode::InputKey(FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey, EInputEvent InEvent)
{
	const bool isLeftButtonDown = (InKey == EKeys::LeftMouseButton && InEvent != IE_Released) || InViewport->KeyState(EKeys::LeftMouseButton);
	const bool isCtrlDown = ((InKey == EKeys::LeftControl || InKey == EKeys::RightControl) && InEvent != IE_Released) || InViewport->KeyState(EKeys::LeftControl) || InViewport->KeyState(EKeys::RightControl);
	const bool isShiftDown = ((InKey == EKeys::LeftShift || InKey == EKeys::RightShift) && InEvent != IE_Released) || InViewport->KeyState(EKeys::LeftShift) || InViewport->KeyState(EKeys::RightShift);
	const bool isAltDown = ((InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt) && InEvent != IE_Released) || InViewport->KeyState(EKeys::LeftAlt) || InViewport->KeyState(EKeys::RightAlt);

	// this is a click
	_shouldProcessInputDelta = !isAltDown && isLeftButtonDown;

	if (InKey == EKeys::F)
	{
		if (focusSelection(InViewportClient)) return true;
	}

	return false;
}

bool FGraphSimulationEdMode::focusSelection(FEditorViewportClient* InViewportClient)
{
	if (!selectedGraphComponent || !lastSelectedNode)
		return false;

	FGraphNode& node = selectedGraphComponent->graph.node(lastSelectedNode);

	UMeshSimulation * meshSim = selectedGraphComponent->simulationManager->simulation<UMeshSimulation>();

	if (!meshSim)
		return false;

	FBox box = meshSim->boundsBoxForNode(lastSelectedNode);

	InViewportClient->FocusViewportOnBox(box, true);

	return true;
}

bool FGraphSimulationEdMode::InputDelta( FEditorViewportClient* InViewportClient, FViewport* InViewport, FVector& InDrag, FRotator& InRot, FVector& InScale )
{
	if (!_shouldProcessInputDelta)
		return false;

	if(!selectedGraphComponent || !lastSelectedNode )
		return false;

	FGraphNode& node = selectedGraphComponent->graph.node( lastSelectedNode );

	FTransform transform = selectedGraphComponent->GetOwner()->GetTransform();

	FVector drag = transform.InverseTransformVector( InDrag );
	FQuat rotation = InRot.Quaternion() * transform.GetRotation().Inverse();


	node.position += drag;
	node.orientation =  rotation * node.orientation;

	UMeshSimulation * meshSimulation = selectedGraphComponent->simulationManager->registerSimulation<UMeshSimulation>();

	meshSimulation->updateInstances();

	return true;
}

void FGraphSimulationEdMode::getRay( FEditorViewportClient* InViewportClient, FVector& start_out, FVector& dir_out, FVector& end_out )
{
	FViewport * viewport = InViewportClient->Viewport;

	const int32 HitX = viewport->GetMouseX();
	const int32 HitY = viewport->GetMouseY();

	FSceneViewFamilyContext viewFamily( FSceneViewFamily::ConstructionValues( InViewportClient->Viewport,
		InViewportClient->GetScene(),
		InViewportClient->EngineShowFlags )
		.SetRealtimeUpdate( InViewportClient->IsRealtime() ) );
	FSceneView* view = InViewportClient->CalcSceneView( &viewFamily );
	FViewportCursorLocation mouseViewportRay( view, (FEditorViewportClient*)viewport->GetClient(), HitX, HitY );

	// we have to translate into the components coordinate frame

	start_out = mouseViewportRay.GetOrigin();
	dir_out = mouseViewportRay.GetDirection();
	end_out = start_out + dir_out * HALF_WORLD_MAX;
}

bool FGraphSimulationEdMode::HandleClick( FEditorViewportClient* InViewportClient, HHitProxy *HitProxy, const FViewportClick &Click )
{
	auto viewport = InViewportClient->Viewport;

	bool selected = 
		[&]() {
		if (!selectedGraphComponent)
		{
			if (!HitProxy)
				return false;

			if (!HitProxy->IsA(HActor::StaticGetType()))
				return false;

			HActor* actorProxy = ((HActor*)HitProxy);

			if (!actorProxy->PrimComponent)
				return false;

			UInstancedStaticMeshComponent * mesh = const_cast<UInstancedStaticMeshComponent*>(Cast<UInstancedStaticMeshComponent>(actorProxy->PrimComponent));

			if (!mesh)
				return false;


			selectedGraphComponent = mesh->GetOwner()->FindComponentByClass<UGraphComponent>();

			if (!selectedGraphComponent)
				return false;
		}

		UMeshSimulation * meshSimulation = selectedGraphComponent->simulationManager->registerSimulation < UMeshSimulation>();

		TArray<UInstancedStaticMeshComponent*> meshes = meshSimulation->instancedStaticMeshes();

		FVector traceStart, traceDirection, traceEnd;
		getRay(InViewportClient, traceStart, traceDirection, traceEnd);

		float nearestDistance = std::numeric_limits<float>::max();
		UInstancedStaticMeshComponent * nearestISMC = nullptr;
		int32 nearestInstance = -1;

		for (auto mesh : meshes)
		{
			// needed in the editor
			if (!mesh->IsPhysicsStateCreated())
			{
				mesh->SetCollisionEnabled(ECollisionEnabled::QueryOnly);
				mesh->CreatePhysicsState();
			}

			FHitResult traceResult(1.0f);

			// doing a trace on the ISMC directly doesn't work in the editor for whatever reason
			// this is just a UI so brute-force it
			for (auto body : mesh->InstanceBodies)
			{
				if (body->LineTrace(traceResult, traceStart, traceEnd, true, false))
				{
					if (traceResult.Distance < nearestDistance)
					{
						nearestISMC = mesh;
						nearestInstance = traceResult.Item;
						nearestDistance = traceResult.Distance;
					}
				}
			}
		}

		// we have a hit! add it to our selection
		if (nearestISMC)
		{
			FGraphNodeHandle hitNode = meshSimulation->nodeForInstance(nearestISMC, nearestInstance);

			addOrRemoveToSelection(hitNode);

			return true;
		}
		else
			return false;
	}();

	if (!selected)
	{
		// we suck, no hit
		clearSelection();

		selectedGraphComponent = nullptr;

		return false;
	}
	else
		return true;
}

void FGraphSimulationEdMode::clearSelection()
{
	if (selectedGraphComponent)
		selectedGraphComponent->graph.unrealEditorSelection.Empty();

	selectedGraphComponent = nullptr;
	lastSelectedNode = FGraphNodeHandle::null;

	if (graphToolKit.IsValid())
		graphToolKit->selectionChanged();
}

void FGraphSimulationEdMode::addOrRemoveToSelection(FGraphNodeHandle handle)
{
	if (selectedGraphComponent && handle)
	{
		auto& selection = selectedGraphComponent->graph.unrealEditorSelection;

		// update the new graph's selection
		if (selection.Contains(handle))
		{
			selection.Remove(handle);

			if (selection.Num() > 0)
				lastSelectedNode = *selection.CreateConstIterator();
			else
				lastSelectedNode = FGraphNodeHandle::null;
		}
		else
		{
			selection.Add(handle);
			lastSelectedNode = handle;
		}
	}
	else
		lastSelectedNode = FGraphNodeHandle::null;

	if (graphToolKit.IsValid())
		graphToolKit->selectionChanged();
}

bool FGraphSimulationEdMode::IsSelectionAllowed( AActor* InActor, bool bInSelection ) const
{
	UGraphComponent * graph = InActor->FindComponentByClass<UGraphComponent>();

	if(graph)
		return true;

	return false;
}

void FGraphSimulationEdMode::SelectionChanged()
{
	clearSelection();
}

FVector FGraphSimulationEdMode::GetWidgetLocation() const
{
	if( !selectedGraphComponent || !lastSelectedNode )
		return FVector::ZeroVector;
	
	FGraphNode& node = selectedGraphComponent->graph.node( lastSelectedNode );

	AActor * actor = selectedGraphComponent->GetOwner();
	
	FTransform transform = actor->GetTransform();

	return transform.TransformPosition(node.position);
}


bool FGraphSimulationEdMode::UsesTransformWidget() const
{
	return true;
}

bool FGraphSimulationEdMode::UsesTransformWidget( FWidget::EWidgetMode CheckMode ) const
{
	return true;

}

bool FGraphSimulationEdMode::UsesToolkits() const
{
	return true;
}






FGraphSimulationToolKit::FGraphSimulationToolKit(class FGraphSimulationEdMode * InOwningMode)
	: editorMode(InOwningMode)
{

}

void FGraphSimulationToolKit::RegisterTabSpawners(const TSharedRef<class FTabManager>& TabManager)
{

}

void FGraphSimulationToolKit::UnregisterTabSpawners(const TSharedRef<class FTabManager>& TabManager)
{

}

void FGraphSimulationToolKit::Init(const TSharedPtr< class IToolkitHost >& InitToolkitHost)
{
	FPropertyEditorModule& propertyModule = FModuleManager::GetModuleChecked<FPropertyEditorModule>("PropertyEditor");

	FDetailsViewArgs detailArgs;
	//FStructureDetailsViewArgs structDetailArgs;

	//FGraph * graph = editorMode->graphComponent ? &editorMode->graphComponent->graph : nullptr;

	//TSharedRef<FStructOnScope> structData(new FStructOnScope(FGraph::StaticStruct(), (uint8*)graph));

	//detailView = propertyModule.CreateStructureDetailView(detailArgs, structDetailArgs, structData);

	detailView = propertyModule.CreateDetailView(detailArgs);

	contentView = constructContentView(detailView);

	FModeToolkit::Init(InitToolkitHost);
}

FName FGraphSimulationToolKit::GetToolkitFName() const
{
	return FName("GraphSimulationToolKit");
}

FText FGraphSimulationToolKit::GetBaseToolkitName() const
{
	return LOCTEXT("Graph Simulation Tool Kit", "Graph Simulation Tool Kit");
}

class FEdMode* FGraphSimulationToolKit::GetEditorMode() const
{
	return editorMode;
}

TSharedPtr<class SVerticalBox> FGraphSimulationToolKit::constructContentView(TSharedPtr<IDetailsView> details)
{
	TSharedPtr<class SVerticalBox> box;

	SAssignNew(box, SVerticalBox)
	+SVerticalBox::Slot()
	.AutoHeight()
	.HAlign(HAlign_Left)
	.VAlign(VAlign_Top)
	[
		SNew(SButton)
		.Text(LOCTEXT("Aggregate and rigidify selection", "Aggregate and rigidify selection"))
		.OnClicked_Lambda([&]() ->FReply {
			UGraphComponent * graphComponent = editorMode->selectedGraphComponent;

			if (!graphComponent) return FReply::Handled();

			FGraph& graph = graphComponent->graph;
	
			graph.beginTransaction();

			auto arraySelection = graph.unrealEditorSelection.Array();

			auto& rigid = FFlexRigidBodyObject::createRigidBody(graph, arraySelection);

			FMLAggregateNO::aggregateNodes(arraySelection, graph, rigid.nodeHandle());

			graph.endTransaction();

			selectionChanged();

			return FReply::Handled();
		})
	]
	+ SVerticalBox::Slot()
	.HAlign(HAlign_Fill)
	.VAlign(VAlign_Fill)
	[
		detailView.ToSharedRef()
	];

	return box;
}

TSharedPtr<class SWidget> FGraphSimulationToolKit::GetInlineContent() const
{
	return contentView;
	//return detailView->GetWidget();
}

void FGraphSimulationToolKit::selectionChanged()
{
	if (!detailView.IsValid())
		return;

	if (editorMode->selectedGraphComponent)
		detailView->SetObject(editorMode->selectedGraphComponent, true);
	else
		detailView->SetObject(nullptr, true);


	//FGraphNode& node = editorMode->graphComponent->graph.node(editorMode->selectedNode);
	//if (!node.isValid())
	//	return;

	//TSharedRef<FStructOnScope> nodeData(new FStructOnScope(FGraphNode::StaticStruct(), (uint8*)&node));

	//detailView->SetStructureData(nodeData);
}

#undef LOCTEXT_NAMESPACE

#endif