// Copyright 2019, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#if WITH_EDITOR

#include "GraphPropertyCustomization.h"
#include "IPropertyUtilities.h"

#define LOCTEXT_NAMESPACE "FGraphPropertyCustomization"

// ----------------------------------------------------------------------
// FGraph - customization
// ----------------------------------------------------------------------

TSharedRef<IPropertyTypeCustomization> FGraphPropteryCustomization::MakeInstance()
{
	return MakeShareable(new FGraphPropteryCustomization);
}

void FGraphPropteryCustomization::CustomizeHeader(
	TSharedRef<class IPropertyHandle> structPropertyHandle,
	class FDetailWidgetRow& headerRow, 
	IPropertyTypeCustomizationUtils& structCustomizationUtils)
{
	TArray<void*> rawData;

	structPropertyHandle->AccessRawData(rawData);

	FGraph * graph = reinterpret_cast<FGraph*>(rawData.Top());

	TSharedPtr<SHorizontalBox> HorizontalBox;

	headerRow.NameContent()
	[
		SNew(STextBlock).Text(LOCTEXT("Nodes", "Nodes"))
	]
	.ValueContent()
	[
		SAssignNew(HorizontalBox, SHorizontalBox)
	];


	if (!graph)
		return;

	TSharedPtr<IPropertyUtilities> propertyUtilities = structCustomizationUtils.GetPropertyUtilities();

	// remove node button
	HorizontalBox->AddSlot()
	[
		SNew(SButton)
		.Text(LOCTEXT("Add Node", "Add Node"))
		.OnClicked_Lambda([graph, propertyUtilities]() ->FReply {
			graph->beginTransaction();

			FGraphNodeHandle handle(graph->addNode(FVector::ZeroVector));

			graph->endTransaction();

			FGraphNode& node = graph->node(handle);

			TSharedRef<FStructOnScope> nodeData(new FStructOnScope(FGraphNode::StaticStruct(), (uint8*)&node));

			graph->unrealEditorSelection.Add(handle);

			propertyUtilities->ForceRefresh();

			return FReply::Handled();
		})
	];
}

void FGraphPropteryCustomization::CustomizeChildren(
	TSharedRef<class IPropertyHandle> structPropertyHandle,
	class IDetailChildrenBuilder& structBuilder,
	IPropertyTypeCustomizationUtils& structCustomizationUtils)
{
	TArray<void*> rawData;

	structPropertyHandle->AccessRawData(rawData);

	FGraph * graph = reinterpret_cast<FGraph*>(rawData.Top());

	for (FGraphNode& node : graph->allNodes)
	{
		if (!node.isValid())
			continue;

		// if we have a selection, only show the selection
		// if there is no selection, show everything
		if (graph->unrealEditorSelection.Num() > 0 && !graph->unrealEditorSelection.Contains(node.handle()))
			continue;

		TSharedRef<class FGraphNode_DetailsCustomNodeBuilder> nodeBuilder = MakeShareable(new FGraphNode_DetailsCustomNodeBuilder());
		nodeBuilder->_graph = graph;
		nodeBuilder->_nodeHandle = node.handle();

		structBuilder.AddCustomBuilder(nodeBuilder);
	}
}

void FGraphPropteryCustomization::OnChildValueChanged()
{

}


// ----------------------------------------------------------------------
// FGraphNode - customization - without components
// ----------------------------------------------------------------------

TSharedRef<IPropertyTypeCustomization> FGraphNodePropteryCustomization::MakeInstance()
{
	return MakeShareable(new FGraphNodePropteryCustomization);
}

void FGraphNodePropteryCustomization::CustomizeHeader(
	TSharedRef<class IPropertyHandle> structPropertyHandle, 
	class FDetailWidgetRow& headerRow, 
	IPropertyTypeCustomizationUtils& structCustomizationUtils)
{
	const bool bDisplayResetToDefault = false;
	const FText DisplayNameOverride = FText::GetEmpty();
	const FText DisplayToolTipOverride = FText::GetEmpty();

	headerRow.NameContent()
	[
		structPropertyHandle->CreatePropertyNameWidget(DisplayNameOverride, DisplayToolTipOverride, bDisplayResetToDefault)
	];
}

void FGraphNodePropteryCustomization::CustomizeChildren(
	TSharedRef<class IPropertyHandle> structPropertyHandle,
	class IDetailChildrenBuilder& structBuilder, 
	IPropertyTypeCustomizationUtils& structCustomizationUtils)
{
	// Add the regular properties, except for components
	uint32 n;
	structPropertyHandle->GetNumChildren(n);

	FName componentsPropertyName = GET_MEMBER_NAME_CHECKED(FGraphNode, components);

	for (uint32 i = 0; i < n; ++i)
	{
		auto childHandle = structPropertyHandle->GetChildHandle(i).ToSharedRef();

		UProperty * uProperty = childHandle->GetProperty();
		FName name = uProperty->GetFName();

		// don't show the components property array
		if (componentsPropertyName == name)
			continue;

		structBuilder.AddProperty(childHandle);
	}
}

FGraph* FGraphNodePropteryCustomization::getGraph(TSharedRef<class IPropertyHandle> structPropertyHandle)
{
	{
		auto parentHandle = structPropertyHandle->GetParentHandle();

		// make sure we are in an array (FGraph::allNodes)
		auto arrayHandle = parentHandle->AsArray();

		if (!arrayHandle.Get())
			return nullptr;

		// Now, we'll grab the graph
		auto graphHandle = parentHandle->GetParentHandle();

		if (!graphHandle->GetProperty()->IsA(UStructProperty::StaticClass()))
			return nullptr;

		// make sure we have a FGraph
		UStructProperty * parentStruct = Cast<UStructProperty>(graphHandle->GetProperty());

		if (!parentStruct)
			return nullptr;

		if (parentStruct->Struct != FGraph::StaticStruct())
			return nullptr;

		TArray<void*> rawData;

		parentHandle->AccessRawData(rawData);

		return reinterpret_cast<FGraph*>(rawData.Top());
	}
}

// ----------------------------------------------------------------------
// FGraphNode - details custom node builder - (can remove nodes from a graph)
// ----------------------------------------------------------------------

void FGraphNode_DetailsCustomNodeBuilder::SetOnRebuildChildren(FSimpleDelegate InOnRegenerateChildren)
{
	onRegenerateChildren = InOnRegenerateChildren;
}

void FGraphNode_DetailsCustomNodeBuilder::GenerateHeaderRowContent(FDetailWidgetRow& nodeRow)
{
	FText name = FText::FromName(FGraphNode::StaticStruct()->GetFName());

	TSharedPtr<SHorizontalBox> HorizontalBox;
	
	nodeRow.NameContent()
	[
		SNew(STextBlock).Text(name)
	]
	.ValueContent()
	[
		SAssignNew(HorizontalBox, SHorizontalBox)
	];


	if (!_graph)
		return;

	if (!_nodeHandle)
		return;

	if (!_graph->node(_nodeHandle).isValid())
		return;

	// remove node button
	HorizontalBox->AddSlot()
	[
		SNew(SButton)
		.Text(LOCTEXT("Remove Node", "Remove Node"))
		.OnClicked_Lambda([&]() ->FReply {
			_graph->beginTransaction();

			_graph->removeNode(_nodeHandle.index);

			_graph->endTransaction();

			_graph->unrealEditorSelection.Remove(_nodeHandle);

			onRegenerateChildren.ExecuteIfBound();

			return FReply::Handled();
		})
	];
}

void FGraphNode_DetailsCustomNodeBuilder::GenerateChildContent(IDetailChildrenBuilder& structBuilder)
{
	FGraph * graph = _graph;

	if (!graph)
		return;

	if (!_nodeHandle)
		return;

	FGraphNode& node = graph->node(_nodeHandle);

	if (!node.isValid())
		return;

	TSharedRef<FStructOnScope> nodeData(new FStructOnScope(FGraphNode::StaticStruct(), (uint8*)&node));

	IDetailPropertyRow * structureRow = structBuilder.AddExternalStructure(nodeData);


	TSharedRef<class FGraphNode_components_DetailCustonNodeBuilder> componentsBuilder = MakeShareable(new FGraphNode_components_DetailCustonNodeBuilder());
	componentsBuilder->_graph = graph;
	componentsBuilder->_nodeHandle = node.handle();

	structBuilder.AddCustomBuilder(componentsBuilder);
}

void FGraphNode_DetailsCustomNodeBuilder::Tick(float DeltaTime)
{
	throw std::logic_error("The method or operation is not implemented.");
}

bool FGraphNode_DetailsCustomNodeBuilder::RequiresTick() const
{
	return false;
}

bool FGraphNode_DetailsCustomNodeBuilder::InitiallyCollapsed() const
{
	return false;
}

FName FGraphNode_DetailsCustomNodeBuilder::GetName() const
{
	return FName();
}

// ----------------------------------------------------------------------
// FGraphNode_components - details custom node builder - can access the graph to display the node's components
// ----------------------------------------------------------------------

void FGraphNode_components_DetailCustonNodeBuilder::SetOnRebuildChildren(FSimpleDelegate InOnRegenerateChildren)
{
	onRegenerateChildren = InOnRegenerateChildren;
}

void FGraphNode_components_DetailCustonNodeBuilder::GenerateHeaderRowContent(FDetailWidgetRow& nodeRow)
{
	const bool bDisplayResetToDefault = false;
	const FText DisplayNameOverride = FText::GetEmpty();
	const FText DisplayToolTipOverride = FText::GetEmpty();

	TSharedPtr<SButton> loadComponentsButton;

	nodeRow.NameContent()
	[
		SNew(STextBlock)
		.Text(LOCTEXT("Components", "Components"))
	]
	.ValueContent()
	[
		SAssignNew(loadComponentsButton, SButton)
		.Text(LOCTEXT("Load Components", "Load Components"))
		.OnClicked_Lambda([&, loadComponentsButton]() ->FReply
		{
			_shouldLoadComponents = true;
			onRegenerateChildren.ExecuteIfBound();

			// and hide it
			_loadComponentsButton->SetEnabled(false);
			_loadComponentsButton->SetVisibility(EVisibility::Hidden);

			return FReply::Handled();
		})
	];

	_loadComponentsButton = loadComponentsButton.Get();
}

void FGraphNode_components_DetailCustonNodeBuilder::GenerateChildContent(IDetailChildrenBuilder& structBuilder)
{
	if (!_graph || !_shouldLoadComponents)
		return;

	FGraph * graph = _graph;

	FGraphNode& node = graph->node(_nodeHandle);

	FGraphNodeHandle nodeHandle = _nodeHandle;
	
	// add our components
	{
		IDetailCategoryBuilder * categoryBuilder = &structBuilder.GetParentCategory();

		for (auto componentType : node.components)
		{
			FGraphObject * component = node.component(*graph, componentType);
			UScriptStruct * scriptStruct = FGraphObject::componentStruct(componentType);
			TSharedRef<FStructOnScope> structData(new FStructOnScope(scriptStruct, (uint8*)component));

			if (!component || !scriptStruct)
				continue;

			IDetailPropertyRow * structureRow = structBuilder.AddExternalStructure(structData);

			structureRow->ShouldAutoExpand(false);

			FText name = FText::FromName(scriptStruct->GetFName());

			structureRow->CustomWidget(true).NameContent()
			[
				SNew(STextBlock).Text(name)
			]
			.ValueContent()
			[
				SNew(SButton)
				.Text(LOCTEXT("Remove Component", "Remove Component"))
				.OnClicked_Lambda([&, nodeHandle, graph, componentType]() ->FReply {

					FGraphNode& node = graph->node(nodeHandle);
					node.removeComponent(*graph, componentType);

					_shouldLoadComponents = true;
					onRegenerateChildren.ExecuteIfBound();

					return FReply::Handled();
				})
			];
		}
	}

	// create a popup menu of component classes to create
	{
		// populate our menu delegate
		TSharedPtr< const FUICommandList > commandsEmpty;
		FMenuBuilder componentTypesBuilder(true, commandsEmpty);
		{
			TArray<UScriptStruct*> componentStructs;
			for (TObjectIterator<UScriptStruct> it; it; ++it)
			{
				if (it->IsChildOf(FGraphObject::StaticStruct()))
					componentStructs.Add(*it);
			}

			for (UScriptStruct * componentStruct : componentStructs)
			{
				ComponentType componentType = FGraphObject::componentType(componentStruct);

				// ----------------------
				// Create a component
				// ----------------------
				FUIAction action(FExecuteAction::CreateLambda([&, nodeHandle, graph, componentType]() {
					graph->beginTransaction();

					FGraphNode& localNode = graph->node(nodeHandle);

					localNode.addComponent(*graph, componentType);

					graph->endTransaction();

					_shouldLoadComponents = true;
					onRegenerateChildren.ExecuteIfBound();
				}));

				componentTypesBuilder.AddMenuEntry(
					FText::FromName(componentStruct->GetFName()),
					FText::FromName(componentStruct->GetFName()),
					FSlateIcon(),
					action,
					NAME_None,
					EUserInterfaceActionType::Button);
			}
		}

		// build the menu widget
		structBuilder.AddCustomRow(LOCTEXT("Add Component", "Add Component"))
		.ValueContent()[
			SNew(SComboButton)
			.ButtonContent()
			[
				SNew(STextBlock)
				.Text(LOCTEXT("Add Component", "Add Component"))
			]
			.MenuContent()
			[
				componentTypesBuilder.MakeWidget()
			]
		];
	}
}

void FGraphNode_components_DetailCustonNodeBuilder::Tick(float DeltaTime)
{

}

bool FGraphNode_components_DetailCustonNodeBuilder::RequiresTick() const
{
	return false;
}

bool FGraphNode_components_DetailCustonNodeBuilder::InitiallyCollapsed() const
{
	return true;
}

FName FGraphNode_components_DetailCustonNodeBuilder::GetName() const
{
	return FName();
}


#undef LOCTEXT_NAMESPACE

#endif // WITH_EDITOR