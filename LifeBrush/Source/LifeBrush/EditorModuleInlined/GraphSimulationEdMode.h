// Copyright 2019 Timothy Davison, Inc. All Rights Reserved.
// Based on my code for Ship Editor.

#pragma once

#if WITH_EDITOR

class UGraphComponent;
class UMappedInstancedStaticMesh;
class IStructureDetailsView;

#include "UnrealEd.h" 
#include "Editor.h"
#include "ShipEditorSimulation/Graph.h"

class FGraphSimulationToolKit;

class FGraphSimulationEdMode : public FEdMode
{
public:
	const static FEditorModeID EM_GraphSimulationEdModeId;

public:
	FGraphSimulationEdMode();
	virtual ~FGraphSimulationEdMode();

public:
	// FEdMode interface
	virtual void Enter() override;
	virtual void Exit() override;

	virtual bool InputKey(FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey, EInputEvent InEvent) override;
	virtual bool HandleClick(FEditorViewportClient* InViewportClient, HHitProxy *HitProxy, const FViewportClick &Click) override;
	
	virtual bool InputDelta(FEditorViewportClient* InViewportClient, FViewport* InViewport, FVector& InDrag, FRotator& InRot, FVector& InScale) override;


	
	virtual bool IsSelectionAllowed(AActor* InActor, bool bInSelection) const override;
	virtual void SelectionChanged() override;

	virtual bool AllowWidgetMove()  override { return true; }

	virtual FVector GetWidgetLocation() const override;
	virtual bool UsesTransformWidget() const override;
	virtual bool UsesTransformWidget( FWidget::EWidgetMode CheckMode ) const override;
	virtual bool UsesToolkits() const override;
	virtual bool ShouldDrawWidget() const override { return true; }
	
public:	
	void clearSelection();
	void addOrRemoveToSelection(FGraphNodeHandle handle);

	void getRay(FEditorViewportClient* InViewportClient, FVector& start_out, FVector& dir_out, FVector& end_out);
	bool focusSelection(FEditorViewportClient* InViewportClient);

public:
	UGraphComponent * selectedGraphComponent = nullptr;
	FGraphNodeHandle lastSelectedNode = FGraphNodeHandle::null;

	TSharedPtr<FGraphSimulationToolKit> graphToolKit;

protected:
	bool _shouldProcessInputDelta = false;
};

// We need a toolkit and a separate widget (from the property module)
class FGraphSimulationToolKit : public FModeToolkit
{
public:

	FGraphSimulationToolKit(class FGraphSimulationEdMode * InOwningMode);

	virtual void RegisterTabSpawners(const TSharedRef<class FTabManager>& TabManager) override;
	virtual void UnregisterTabSpawners(const TSharedRef<class FTabManager>& TabManager) override;

	/** Initializes the geometry mode toolkit */
	virtual void Init(const TSharedPtr< class IToolkitHost >& InitToolkitHost) override;

	/** IToolkit interface */
	virtual FName GetToolkitFName() const override;
	virtual FText GetBaseToolkitName() const override;
	virtual class FEdMode* GetEditorMode() const override;
	virtual TSharedPtr<class SWidget> GetInlineContent() const override;

	void selectionChanged();

	TSharedPtr<class SVerticalBox> constructContentView(TSharedPtr<IDetailsView> details);
private:
	class FGraphSimulationEdMode * editorMode;
	TSharedPtr<class SWidget> widget;

	//TSharedPtr<IStructureDetailsView> detailView;

	TSharedPtr<SWidget> contentView;
	TSharedPtr<IDetailsView> detailView;
};

#endif
