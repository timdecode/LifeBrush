//
//  TCODSEditorMode.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-08-22.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

/*
 Some nice resources: https://www.youtube.com/watch?t=2465&v=zg_VstBxDi8
 */

#pragma once

#if WITH_EDITOR

#include <utility>

#include "UnrealEd.h"
#include "Editor/UnrealEd/Public/Toolkits/BaseToolkit.h"
#include "Editor/UnrealEd/Public/EditorModes.h"
#include "Toolkits/ToolkitManager.h"
#include "Editor/LevelEditor/Public/LevelEditor.h"
#include "Editor/LevelEditor/Public/LevelEditorActions.h"
#include "Editor/UnrealEd/Public/BSPOps.h"

#include "Editor/BspMode/Public/IBspModeModule.h"

#include "Algorithm/Algorithm.h"

class URegionGrowingComponent;
class AElementActor;

class FTCODSEditorMode : public FEdMode
{
public:
    const static FEditorModeID EM_Geometry;
    
public:
    FTCODSEditorMode();
    virtual ~FTCODSEditorMode();
    
    virtual void Tick(FEditorViewportClient* ViewportClient,float DeltaTime) override;

    virtual bool UsesToolkits() const override;

    virtual void Enter() override;
    virtual void Exit() override;
    
    virtual bool MouseMove( FEditorViewportClient* ViewportClient, FViewport* Viewport, int32 x, int32 y ) override;
    virtual bool CapturedMouseMove( FEditorViewportClient* InViewportClient, FViewport* InViewport, int32 InMouseX, int32 InMouseY ) override;
    virtual bool StartTracking(FEditorViewportClient* InViewportClient, FViewport* InViewport) override;
    virtual bool EndTracking(FEditorViewportClient* InViewportClient, FViewport* InViewport) override;
    virtual bool InputKey( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey, EInputEvent InEvent ) override;
    
    bool LeftMouseUp( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey );
    bool LeftMouseDown( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey );
	bool sketch_mouseDown( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey );
	bool sketch_mouseUp( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey );
	bool LeftMouseClick( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey );

    void getRay( FEditorViewportClient* InViewportClient, FViewport* InViewport, FVector& start_out, FVector& dir_out, FVector& end_out );
   
    virtual bool Select( AActor* InActor, bool bInSelected ) override;
    virtual bool IsSelectionAllowed( AActor* InActor, bool bInSelection ) const override;

    enum class ClickMode { AddSingularity, AddSeed, Stretch, SketchExemplar };
	enum class SketchMode { None, Sketching, Erasing };

    
    void setClickMode(ClickMode mode)
    {
        _clickMode = mode;
        _isDragging = false;
    }
    
    ClickMode clickMode() { return _clickMode; }
    
	void convertToElements();
	void toExemplar();
	void clear();

private:
    void _startEditing();
    void _stopEditing();
    
    auto _vertexForHitPoint(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit) -> std::pair<bool, FVector>;
    
    void _addSeed(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit);
    void _addSingularity(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit);
    
    void _removeSeedOrSingularity(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit);
    
	void _selectRegionGrower( FEditorViewportClient* InViewportClient, FViewport* InViewport );

	URegionGrowingComponent* _hitRegionGrower( FEditorViewportClient* InViewportClient, FViewport* InViewport );
	void _strechClick( FEditorViewportClient* InViewportClient, FViewport* InViewport, UStaticMeshComponent* staticMesh, URegionGrowingComponent* regionGrower );
    bool _strechUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport );
    bool _strechEnd( FEditorViewportClient* InViewportClient, FViewport* InViewport );
    
    void _sketchClick( FEditorViewportClient* InViewportClient, FViewport* InViewport);
    bool _sketchUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport );
	bool _sketchSelectUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport );
	bool _sketchEnd( FEditorViewportClient* InViewportClient, FViewport* InViewport );

	bool _eraseUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport );

	void _updateSelection();

private:
    ClickMode _clickMode = ClickMode::AddSingularity;
    
    bool _isEditing;
    
    FIntPoint _leftMouseDownPoint;
    bool _canBeLeftClick;

    FPlane _clickPlane;
    bool _isDragging = false;

	SketchMode _sketchMode = SketchMode::None;

    FVector _lineTraceStart;
    URegionGrowingComponent * _regionGrower = nullptr;

	URegionGrowingComponent * _selectionLastGrower = nullptr;

	Algorithm::ExampleSelectionPtr exampleSelectionHandle = nullptr;

	std::unordered_set<AElementActor*> _exampleSelection_actors;

	bool _selectionDidChange = false;
};

class FTrivialConnectionsTool : public FModeTool
{
public:
    FTrivialConnectionsTool();
    
    virtual bool StartModify() override;
    virtual bool EndModify() override;
    
    virtual void StartTrans() override;
    virtual void EndTrans() override;
};

class FTrivialConnectionsToolkit : public FModeToolkit
{
public:

    FTrivialConnectionsToolkit(class FTCODSEditorMode * InOwningMode);
    
    virtual void RegisterTabSpawners(const TSharedRef<class FTabManager>& TabManager) override;
    virtual void UnregisterTabSpawners(const TSharedRef<class FTabManager>& TabManager) override;
    
    /** Initializes the geometry mode toolkit */
    virtual void Init(const TSharedPtr< class IToolkitHost >& InitToolkitHost) override;
    
    /** IToolkit interface */
    virtual FName GetToolkitFName() const override;
    virtual FText GetBaseToolkitName() const override;
    virtual class FEdMode* GetEditorMode() const override;
    virtual TSharedPtr<class SWidget> GetInlineContent() const override;
    
private:
    class FTCODSEditorMode * editorMode;
    TSharedPtr<class STrivialConnectionWidget> widget;
};

class STrivialConnectionWidget : public SCompoundWidget
{
public:
    SLATE_BEGIN_ARGS( STrivialConnectionWidget ) {}
    SLATE_END_ARGS();
    
    ~STrivialConnectionWidget();
    
    void Construct(const FArguments& InArgs, TSharedRef<FTrivialConnectionsToolkit> InParentToolkit, class FTCODSEditorMode* InOwningMode);
    
private:
    ECheckBoxState IsRadioChecked( FTCODSEditorMode::ClickMode buttonID ) const
    {
        if( editorMode == nullptr )
            return ECheckBoxState::Unchecked;
        
        return (editorMode->clickMode() == buttonID) ? ECheckBoxState::Checked : ECheckBoxState::Unchecked;
    }
    
private:
    FTCODSEditorMode * editorMode;
};



#endif // WITH_EDITOR
