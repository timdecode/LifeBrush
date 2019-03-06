// Copyright 2019, Timothy Davison. All rights reserved.

#pragma once

#if WITH_EDITOR

#include "IPropertyTypeCustomization.h"
#include "IDetailCustomization.h"
#include "IDetailCustomNodeBuilder.h"

#include "ShipEditorSimulation/Graph.h"

class IPropertyHandle;

struct FGraph;

// ----------------------------------------------------------------------
// FGraph - customization
// ----------------------------------------------------------------------

class FGraphPropteryCustomization : public IPropertyTypeCustomization
{
public:

	static TSharedRef<IPropertyTypeCustomization> MakeInstance();

	virtual void CustomizeHeader(TSharedRef<class IPropertyHandle> StructPropertyHandle, class FDetailWidgetRow& HeaderRow, IPropertyTypeCustomizationUtils& StructCustomizationUtils) override;
	virtual void CustomizeChildren(TSharedRef<class IPropertyHandle> StructPropertyHandle, class IDetailChildrenBuilder& StructBuilder, IPropertyTypeCustomizationUtils& StructCustomizationUtils) override;

protected:
	void OnChildValueChanged();
};

// ----------------------------------------------------------------------
// FGraphNode - customization - without components
// ----------------------------------------------------------------------

// Property customization for FGraphNode, but it doesn't show the components property.
class FGraphNodePropteryCustomization : public IPropertyTypeCustomization
{
public:

	static TSharedRef<IPropertyTypeCustomization> MakeInstance();

	virtual void CustomizeHeader(TSharedRef<class IPropertyHandle> StructPropertyHandle, class FDetailWidgetRow& HeaderRow, IPropertyTypeCustomizationUtils& StructCustomizationUtils) override;
	virtual void CustomizeChildren(TSharedRef<class IPropertyHandle> StructPropertyHandle, class IDetailChildrenBuilder& StructBuilder, IPropertyTypeCustomizationUtils& StructCustomizationUtils) override;

	// Get's the graph, if this property customization as embedded in a FGraph::allNodes array.
	FGraph * getGraph(TSharedRef<class IPropertyHandle> structPropertyHandle);
};

// ----------------------------------------------------------------------
// FGraphNode - details custom node builder - (can remove nodes from a graph)
// ----------------------------------------------------------------------

// Custom builder for FGraphNode that displays a FGraphNodePropertyCustomization and FGraphNode_components_DetailCustonNodeBuilder.
// It's linked to the graph, so it has access to the information that it needs to display the components. The FGraphNodePropertyCustomization
// does not. We could derive that information, however, using IDetailsCustomNodeBuilder also allows us to efficiently reload contents.
class FGraphNode_DetailsCustomNodeBuilder : 
	public IDetailCustomNodeBuilder,
	public TSharedFromThis<FGraphNode_DetailsCustomNodeBuilder>
{
public:
	FGraph * _graph;
	FGraphNodeHandle _nodeHandle;

	FSimpleDelegate onRegenerateChildren;

public:
	virtual void SetOnRebuildChildren(FSimpleDelegate InOnRegenerateChildren) override;

	virtual void GenerateHeaderRowContent(FDetailWidgetRow& NodeRow) override;
	virtual void GenerateChildContent(IDetailChildrenBuilder& ChildrenBuilder) override;

	virtual void Tick(float DeltaTime) override;
	virtual bool RequiresTick() const override;
	virtual bool InitiallyCollapsed() const override;

	virtual FName GetName() const override;

protected:
};

// ----------------------------------------------------------------------
// FGraphNode_components - details custom node builder - can access the graph to display the node's components
// ----------------------------------------------------------------------

// Custom builder for a FGraphNode that displays the components. The components are not displayed until the Load Components
// button is pressed.
class FGraphNode_components_DetailCustonNodeBuilder : 
	public IDetailCustomNodeBuilder, 
	public TSharedFromThis<FGraphNode_components_DetailCustonNodeBuilder>
{
public:
	FGraph * _graph;
	FGraphNodeHandle _nodeHandle;

	FSimpleDelegate onRegenerateChildren;

public:
	virtual void SetOnRebuildChildren(FSimpleDelegate InOnRegenerateChildren) override;

	virtual void GenerateHeaderRowContent(FDetailWidgetRow& NodeRow) override;
	virtual void GenerateChildContent(IDetailChildrenBuilder& ChildrenBuilder) override;

	virtual void Tick(float DeltaTime) override;
	virtual bool RequiresTick() const override;
	virtual bool InitiallyCollapsed() const override;

	virtual FName GetName() const override;

protected:
	bool _shouldLoadComponents = false;
	SButton * _loadComponentsButton; // this is a ptr on purpose and not a shared ptr, so we can use it in the lambda
};



#endif // WITH_EDITOR