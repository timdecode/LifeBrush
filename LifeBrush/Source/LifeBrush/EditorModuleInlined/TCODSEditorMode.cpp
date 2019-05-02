//
//  TCODSEditorMode.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-08-22.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "RuntimeMeshComponent.h"

#if WITH_EDITOR

#include "TCODSEditorMode.h"
#include <limits>
#include "StaticMeshResources.h"
#include "ElementActor.h"

#ifdef __APPLE__
#include <dlfcn.h>
#include <cstring>
#include <string>
#endif

#include "LifeBrush/RegionGrowingComponent.h"


// -----------------------------------------------------------------------------
/// Module
// -----------------------------------------------------------------------------


const FEditorModeID  FTCODSEditorMode::EM_Geometry(TEXT("EM_FTCODSEditorMode"));



// -----------------------------------------------------------------------------
/// Editor Mode
// -----------------------------------------------------------------------------

FTCODSEditorMode::FTCODSEditorMode()
{
    Tools.Add( new FTrivialConnectionsTool());
    SetCurrentTool(MT_GeometryModify);
}

FTCODSEditorMode::~FTCODSEditorMode()
{
    
}

bool FTCODSEditorMode::UsesToolkits() const
{
    return true;
}

void FTCODSEditorMode::Enter()
{
    FEdMode::Enter();
    
    if( !Toolkit.IsValid() )
    {
        Toolkit = MakeShareable(new FTrivialConnectionsToolkit(this));
        Toolkit->Init(Owner->GetToolkitHost());
    }
}

void FTCODSEditorMode::Exit()
{
    if( Toolkit.IsValid() )
    {
        FToolkitManager::Get().CloseToolkit(Toolkit.ToSharedRef());
        Toolkit.Reset();
    }
    
    FEdMode::Exit();
}

void FTCODSEditorMode::Tick(FEditorViewportClient* ViewportClient,float DeltaTime)
{
}


bool FTCODSEditorMode::MouseMove( FEditorViewportClient* ViewportClient, FViewport* Viewport, int32 x, int32 y )
{
    // We only care about perspective viewports
    if( ViewportClient->IsPerspective() )
    {

    }
    
    return false;
}

bool FTCODSEditorMode::CapturedMouseMove( FEditorViewportClient* InViewportClient, FViewport* InViewport, int32 InMouseX, int32 InMouseY )
{
    // We only care about perspective viewports
    if( InViewportClient->IsPerspective() && InViewportClient->EngineShowFlags.ModeWidgets )
    {
        if(_isDragging && _clickMode == ClickMode::Stretch)
            return _strechUpdate( InViewportClient, InViewport );
        else if(_clickMode == ClickMode::SketchExemplar && _sketchMode == SketchMode::Sketching)
            return _sketchUpdate( InViewportClient, InViewport );
		else if(_clickMode == ClickMode::SketchExemplar && _sketchMode == SketchMode::Erasing)
			return _eraseUpdate( InViewportClient, InViewport );


        FSceneViewFamilyContext viewFamily( FSceneViewFamily::ConstructionValues(
                                                                                 InViewportClient->Viewport,
                                                                                 InViewportClient->GetScene(),
                                                                                 InViewportClient->EngineShowFlags)
                                           .SetRealtimeUpdate( InViewportClient->IsRealtime() ));
        FSceneView* view = InViewportClient->CalcSceneView( &viewFamily );
        FViewportCursorLocation mouseViewportRay( view, (FEditorViewportClient*)InViewport->GetClient(), InMouseX, InMouseY );
        
        // capture mouse clicks
        FIntPoint mousePoint;
        InViewport->GetMousePos(mousePoint);
        
        if( (mousePoint - _leftMouseDownPoint).SizeSquared() > 4 )
            _canBeLeftClick = false;
        
        return true;
    }
    
    return false;
}

/** FEdMode: Called when a mouse button is pressed */
bool FTCODSEditorMode::StartTracking(FEditorViewportClient* InViewportClient, FViewport* InViewport)
{
    return true;
}



/** FEdMode: Called when the a mouse button is released */
bool FTCODSEditorMode::EndTracking(FEditorViewportClient* InViewportClient, FViewport* InViewport)
{
    return true;
}

/** FEdMode: Called when a key is pressed */
bool FTCODSEditorMode::InputKey( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey, EInputEvent InEvent )
{
    if( InEvent == IE_Pressed && InKey == EKeys::LeftMouseButton )
    {
        InViewport->GetMousePos(_leftMouseDownPoint);
        _canBeLeftClick = true;
    }
    
	if(_clickMode == ClickMode::SketchExemplar)
	{
		bool consume = false;

		if(InKey == EKeys::LeftMouseButton && InEvent == IE_Pressed)
		{
			_sketchMode = SketchMode::Sketching;

			sketch_mouseDown( InViewportClient, InViewport, InKey );

			consume = true;
		}
		else if(InKey == EKeys::LeftMouseButton && InEvent == IE_Released)
		{
			_sketchMode = SketchMode::None;

			sketch_mouseUp( InViewportClient, InViewport, InKey );

			consume = true;
		}

		return consume;
	}
	else if( _clickMode == ClickMode::Stretch )
    {

        if(InEvent == IE_Pressed && InKey == EKeys::LeftMouseButton)
            return LeftMouseDown( InViewportClient, InViewport, InKey );
        if(InEvent == IE_Released && InKey == EKeys::LeftMouseButton)
            return LeftMouseUp( InViewportClient, InViewport, InKey );
    }

    if(InEvent == IE_Released && InKey == EKeys::LeftMouseButton && _canBeLeftClick)
        return LeftMouseClick( InViewportClient, InViewport, InKey );


    const bool isLeftButtonDown = ( InKey == EKeys::LeftMouseButton && InEvent != IE_Released ) || InViewport->KeyState( EKeys::LeftMouseButton );
    const bool isCtrlDown = ( ( InKey == EKeys::LeftControl || InKey == EKeys::RightControl ) && InEvent != IE_Released ) || InViewport->KeyState( EKeys::LeftControl ) || InViewport->KeyState( EKeys::RightControl );
    const bool isShiftDown = ( ( InKey == EKeys::LeftShift || InKey == EKeys::RightShift ) && InEvent != IE_Released ) || InViewport->KeyState( EKeys::LeftShift ) || InViewport->KeyState( EKeys::RightShift );
    const bool isAltDown = ( ( InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt ) && InEvent != IE_Released ) || InViewport->KeyState( EKeys::LeftAlt ) || InViewport->KeyState( EKeys::RightAlt );
    
    if( !isAltDown && isLeftButtonDown && InViewportClient->IsPerspective() && InViewportClient->EngineShowFlags.ModeWidgets)
    {
        const int32 HitX = InViewport->GetMouseX();
        const int32 HitY = InViewport->GetMouseY();
        const HHitProxy* HitProxy = InViewport->GetHitProxy(HitX, HitY);
        
        if (HitProxy && HitProxy->IsA(HActor::StaticGetType()))
        {
            const AActor* ClickedActor = (static_cast<const HActor*>(HitProxy))->Actor;

			return true;
            //USelection& SelectedActors = *Owner->GetSelectedActors();
            //if (SelectedActors.IsSelected(ClickedActor))
            //{
            //    // Clicked actor is currently selected, start painting.
            //    return true;
            //}
            //else
            //    return false;
        }
    }
    
    return false;
}

void FTCODSEditorMode::getRay( FEditorViewportClient* InViewportClient, FViewport* InViewport, FVector& start_out, FVector& dir_out, FVector& end_out )
{
    const int32 HitX = InViewport->GetMouseX();
    const int32 HitY = InViewport->GetMouseY();

    FSceneViewFamilyContext viewFamily( FSceneViewFamily::ConstructionValues( InViewportClient->Viewport,
        InViewportClient->GetScene(),
        InViewportClient->EngineShowFlags )
        .SetRealtimeUpdate( InViewportClient->IsRealtime() ) );
    FSceneView* view = InViewportClient->CalcSceneView( &viewFamily );
    FViewportCursorLocation mouseViewportRay( view, (FEditorViewportClient*)InViewport->GetClient(), HitX, HitY );

    // we have to translate into the components coordinate frame

    start_out = mouseViewportRay.GetOrigin();
    dir_out = mouseViewportRay.GetDirection();
    end_out = start_out + dir_out * HALF_WORLD_MAX;
}


bool FTCODSEditorMode::LeftMouseDown( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey )
{
    const bool isAltDown = (InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt) || InViewport->KeyState( EKeys::LeftAlt ) || InViewport->KeyState( EKeys::RightAlt );

    const bool isLeftButtonDown = (InKey == EKeys::LeftMouseButton) || InViewport->KeyState( EKeys::LeftMouseButton );
    const bool isRightButtonDown = (InKey == EKeys::RightMouseButton) || InViewport->KeyState( EKeys::RightMouseButton );
    const bool isShiftDown = (InKey == EKeys::LeftShift || InKey == EKeys::RightShift) || InViewport->KeyState( EKeys::LeftShift ) || InViewport->KeyState( EKeys::RightShift );

    if(!(isLeftButtonDown || isRightButtonDown))
        return false;

	

    if(InViewportClient->IsPerspective() && InViewportClient->EngineShowFlags.ModeWidgets)
    {
        const int32 HitX = InViewport->GetMouseX();
        const int32 HitY = InViewport->GetMouseY();
        const HHitProxy* HitProxy = InViewport->GetHitProxy( HitX, HitY );

        if(!(HitProxy && HitProxy->IsA( HActor::StaticGetType() )))
            return false;

        const AActor* ClickedActor = (static_cast<const HActor*>(HitProxy))->Actor;
        USelection& SelectedActors = *Owner->GetSelectedActors();

        if(!SelectedActors.IsSelected( ClickedActor ))
            return false;

        URegionGrowingComponent * regionGrower = ClickedActor->FindComponentByClass<URegionGrowingComponent>();

        if(regionGrower == nullptr)
            return false;

        UStaticMeshComponent * staticMesh = regionGrower->staticMeshComponent();
        if(staticMesh == nullptr)
            return false;

        const bool isAltDown = (InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt) || InViewport->KeyState( EKeys::LeftAlt ) || InViewport->KeyState( EKeys::RightAlt );

        if(!isShiftDown && _clickMode == ClickMode::Stretch)
        {
            _strechClick(InViewportClient, InViewport, staticMesh, regionGrower );
            return true;
        }

    }

    return false;
}

void FTCODSEditorMode::_selectRegionGrower( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	// if we have a region grower, check if we still hit it
	if(_regionGrower && _regionGrower->staticMeshComponent() && _regionGrower->staticMeshComponent()->IsVisible() )
	{
		URegionGrowingComponent * nearest = _hitRegionGrower( InViewportClient, InViewport );

		if(nearest == _regionGrower)
			return;
	}

	URegionGrowingComponent * grower = _hitRegionGrower( InViewportClient, InViewport );

	if(grower == _regionGrower)
		return;

	// make a new one
	_regionGrower = grower;

	exampleSelectionHandle = nullptr;


	if(_regionGrower)
	{
		exampleSelectionHandle = _regionGrower->getExampleSelection();

		_updateSelection();
	}
}

URegionGrowingComponent* FTCODSEditorMode::_hitRegionGrower( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	FVector traceStart;
	FVector traceDirection;
	FVector traceEnd;


	getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

	float minDistance = std::numeric_limits<float>::max();

	URegionGrowingComponent * nearest = nullptr;

	for(TObjectIterator<URegionGrowingComponent> itr; itr; ++itr)
	{
		URegionGrowingComponent * grower = *itr;

		if(grower->staticMeshComponent() && !grower->staticMeshComponent()->IsVisible())
			continue;

		auto meshInterface = grower->meshInterface();

		auto hit = meshInterface->getIntersectionAndFace( traceStart, traceDirection );

		if(!hit.first)
			continue;

		FVector hitPosition = hit.second.point;

		float dSqrd = FVector::DistSquared( traceStart, hitPosition );

		if(dSqrd < minDistance)
		{
			nearest = grower;
			minDistance = dSqrd;
		}
	}


	return nearest;
}

bool FTCODSEditorMode::sketch_mouseDown( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey )
{
	const bool isAltDown = (InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt) || InViewport->KeyState( EKeys::LeftAlt ) || InViewport->KeyState( EKeys::RightAlt );

	const bool isLeftButtonDown = (InKey == EKeys::LeftMouseButton) || InViewport->KeyState( EKeys::LeftMouseButton );
	const bool isRightButtonDown = (InKey == EKeys::RightMouseButton) || InViewport->KeyState( EKeys::RightMouseButton );
	const bool isShiftDown = (InKey == EKeys::LeftShift || InKey == EKeys::RightShift) || InViewport->KeyState( EKeys::LeftShift ) || InViewport->KeyState( EKeys::RightShift );
	const bool isExampleSelectionKeyDown = InViewport->KeyState( EKeys::Z );
		

	if( isRightButtonDown )
		return false;

	if( InViewportClient->IsPerspective() && InViewportClient->EngineShowFlags.ModeWidgets )
	{
		if( !isShiftDown && !isExampleSelectionKeyDown )
			_selectRegionGrower(InViewportClient, InViewport);

		if(!_regionGrower)
			return false;

		if(!isShiftDown )
		{
			_sketchClick( InViewportClient, InViewport );

			return true;
		}

	}

	return false;
}

bool FTCODSEditorMode::sketch_mouseUp( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey )
{
	_sketchMode = SketchMode::None;

	_sketchEnd( InViewportClient, InViewport );

	return true;
}


bool FTCODSEditorMode::LeftMouseUp( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey )
{
    if(_isDragging && _clickMode == ClickMode::Stretch)
        return _strechEnd( InViewportClient, InViewport);


    return false;
}

bool FTCODSEditorMode::LeftMouseClick( FEditorViewportClient* InViewportClient, FViewport* InViewport, FKey InKey)
{
    const bool isAltDown = ( InKey == EKeys::LeftAlt || InKey == EKeys::RightAlt ) || InViewport->KeyState( EKeys::LeftAlt ) || InViewport->KeyState( EKeys::RightAlt );
    
    const bool isLeftButtonDown = ( InKey == EKeys::LeftMouseButton) || InViewport->KeyState( EKeys::LeftMouseButton );
    const bool isRightButtonDown = ( InKey == EKeys::RightMouseButton ) || InViewport->KeyState( EKeys::RightMouseButton );
    const bool isShiftDown = ( InKey == EKeys::LeftShift || InKey == EKeys::RightShift ) || InViewport->KeyState( EKeys::LeftShift) || InViewport->KeyState( EKeys::RightShift);
    
    if( !(isLeftButtonDown || isRightButtonDown) )
        return false;
    
    if( !isAltDown && InViewportClient->IsPerspective() && InViewportClient->EngineShowFlags.ModeWidgets)
    {
        const int32 HitX = InViewport->GetMouseX();
        const int32 HitY = InViewport->GetMouseY();
        const HHitProxy* HitProxy = InViewport->GetHitProxy(HitX, HitY);
        
        if( !(HitProxy && HitProxy->IsA(HActor::StaticGetType())) )
            return false;
        
        const AActor* ClickedActor = (static_cast<const HActor*>(HitProxy))->Actor;
        //USelection& SelectedActors = *Owner->GetSelectedActors();
        //
        //if( !SelectedActors.IsSelected(ClickedActor) )
        //    return false;
        
        URegionGrowingComponent * regionGrower = ClickedActor->FindComponentByClass<URegionGrowingComponent>();
        
        if( regionGrower == nullptr )
            return false;
        
        UStaticMeshComponent * staticMesh = regionGrower->staticMeshComponent();
        if( staticMesh == nullptr )
            return false;
        
        FVector traceStart;
        FVector traceDirection;
        FVector traceEnd;

        getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );
        
        // our physics state could be invalid for the raycast, therefore we might need to do this:
        // https://github.com/EpicGames/UnrealEngine/blob/8a80b5541f69a79abf5855668f39e1d643717600/Engine/Source/Editor/MeshPaint/Private/MeshPaintStaticMeshAdapter.cpp
        // but in the meantime lets just force physics update
        staticMesh->RecreatePhysicsState();
        
        
        FHitResult traceResult(1.0f);
        
        FCollisionQueryParams params(FName(TEXT("KeyUp")), true);
        params.bReturnFaceIndex = true;
        
        if( !staticMesh->LineTraceComponent(traceResult, traceStart, traceEnd, params) )
            return false;

        if( !isShiftDown )
        {
            if( _clickMode == ClickMode::AddSingularity )
                _addSingularity(staticMesh, regionGrower, traceResult);
            else if( _clickMode == ClickMode::AddSeed )
                _addSeed(staticMesh, regionGrower, traceResult);
        }
        else
        {
			if(_clickMode == ClickMode::AddSingularity)
				_removeSeedOrSingularity( staticMesh, regionGrower, traceResult );

		}
    }
    
    return false;
}


/** FEdMode: Handling SelectActor */
bool FTCODSEditorMode::Select( AActor* InActor, bool bInSelected )
{
    TInlineComponentArray<URegionGrowingComponent*> regionGrowers;
    InActor->GetComponents<URegionGrowingComponent>(regionGrowers);
    for(const auto& grower : regionGrowers)
    {
        if (grower != nullptr)
        {
        }
    }
    
    return false;
}

/** FEdMode: Check to see if an actor can be selected in this mode - no side effects */
bool FTCODSEditorMode::IsSelectionAllowed( AActor* InActor, bool bInSelection ) const
{
    TInlineComponentArray<URegionGrowingComponent*> regionGrowers;
    URegionGrowingComponent * grower = InActor->FindComponentByClass<URegionGrowingComponent>();
    
    if( grower != nullptr )
        return true;
        
	if(InActor->IsA( AElementActor::StaticClass() ))
		return true;

    return false;
}

void FTCODSEditorMode::convertToElements()
{
	USelection& selection = *Owner->GetSelectedActors();

	TArray<AActor*> selectedActors;
	selection.GetSelectedObjects<AActor>( selectedActors );

	for(AActor * actor : selectedActors)
	{
		FTransform transform = actor->GetActorTransform();

		AElementActor * elementActor = GetWorld()->SpawnActorAbsolute<AElementActor>( AElementActor::StaticClass(), transform );

		if(actor->GetParentActor())
			elementActor->AttachToActor( actor->GetParentActor(), FAttachmentTransformRules::KeepRelativeTransform );

		// copy any static mesh properties
		UStaticMeshComponent * mesh = actor->FindComponentByClass<UStaticMeshComponent>();

		if(mesh)
		{
			UStaticMeshComponent * elementMesh = elementActor->GetStaticMeshComponent();

			elementMesh->SetStaticMesh( mesh->GetStaticMesh() );

			auto materials = mesh->GetMaterials();

			for(int32 i = 0; i < materials.Num(); ++i)
				elementMesh->SetMaterial( i, materials[i] );
		}
	}
}

void FTCODSEditorMode::toExemplar()
{
	if(_regionGrower)
		_regionGrower->generateExemplar();
}

void FTCODSEditorMode::clear()
{
	if(_regionGrower)
		_regionGrower->ClearOutput();
}

void FTCODSEditorMode::_startEditing()
{
    if( _isEditing )
        return;
    
    _isEditing = true;
}

void FTCODSEditorMode::_stopEditing()
{
    if( !_isEditing )
        return;
    
    _isEditing = false;
}

auto FTCODSEditorMode::_vertexForHitPoint(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit) -> std::pair<bool, FVector>
{
    if( staticMeshComponent == nullptr ) return {false, FVector::ZeroVector};
    
    UStaticMesh * staticMesh = staticMeshComponent->GetStaticMesh();
    
    if( staticMesh->RenderData == nullptr ) return {false, FVector::ZeroVector};
    if( staticMesh->RenderData->LODResources.Num() == 0 ) return {false, FVector::ZeroVector};
    
    FStaticMeshLODResources& resource = staticMesh->RenderData->LODResources[0];
    
    FPositionVertexBuffer& vertexBuffer = resource.PositionVertexBuffer;
    FRawStaticIndexBuffer& indexBuffer = resource.IndexBuffer;
    
    // we'll work in the coordinate frame of the static mesh
    const FVector hitPoint = staticMeshComponent->GetComponentTransform().InverseTransformPosition(hit.ImpactPoint);
    
    auto vertexCount = vertexBuffer.GetNumVertices();
    
    decltype(vertexCount) minI = 0;
    float minDistance = std::numeric_limits<float>::max();
    
    for( decltype(vertexCount) i = 0; i < vertexCount; ++i )
    {
        const FVector v = vertexBuffer.VertexPosition(i);
        
        const float d = FVector::DistSquared(hitPoint, v);
        if( d < minDistance )
        {
            minDistance = d;
            minI = i;
        }
    }
    
    FVector vOut = vertexBuffer.VertexPosition(minI);
    
    return {true, vOut};
}

void FTCODSEditorMode::_addSingularity(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit)
{
    auto vertex = _vertexForHitPoint(staticMeshComponent, regionGrower, hit);
    
    if( vertex.first )
        regionGrower->addSingularity(vertex.second, 1.0);
}

void FTCODSEditorMode::_addSeed(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit)
{
    auto vertex = _vertexForHitPoint(staticMeshComponent, regionGrower, hit);

    if( vertex.first )
        regionGrower->addSeed(vertex.second);
}

void FTCODSEditorMode::_removeSeedOrSingularity(UStaticMeshComponent *staticMeshComponent, URegionGrowingComponent * regionGrower, FHitResult hit)
{
    if( staticMeshComponent == nullptr ) return;
    
    UStaticMesh * staticMesh = staticMeshComponent->GetStaticMesh();
    
    if( staticMesh->RenderData == nullptr ) return;
    if( staticMesh->RenderData->LODResources.Num() == 0 ) return;
    
    // we'll work in the coordinate frame of the static mesh
    const FVector hitPoint = hit.ImpactPoint;
    float sphereRadius = 1.0f;
    
    if( regionGrower->seedDebugSpheres.Num() )
        sphereRadius = regionGrower->seedDebugSpheres[0]->Bounds.GetSphere().W;
    else if( regionGrower->singularityDebugSpheres.Num() )
        sphereRadius = regionGrower->singularityDebugSpheres[0]->Bounds.GetSphere().W;
    
    const FSphere hitSphere = FSphere(hit.ImpactPoint,sphereRadius);
    
    int i = 0;
    for( auto seed : regionGrower->seedDebugSpheres )
    {
        auto bounds = seed->Bounds;
        
        if( bounds.GetSphere().Intersects(hitSphere) )
        {
            regionGrower->removeSeedAt(i);
            return;
        }
        
        i++;
    }
    
    i = 0;
    for( auto singularity : regionGrower->singularityDebugSpheres )
    {
        auto bounds = singularity->Bounds;
        
        if( bounds.GetSphere().IsInside(hitPoint) )
        {
            regionGrower->removeSingularityAt(i);
            return;
        }
        
        i++;
    }
}

// -----------------------------------------------------------------------------
// Stretch
// -----------------------------------------------------------------------------

void FTCODSEditorMode::_strechClick( FEditorViewportClient* InViewportClient, FViewport* InViewport, UStaticMeshComponent * staticMesh, URegionGrowingComponent * regionGrower)
{
    FVector traceStart;
    FVector traceDirection;
    FVector traceEnd;

    getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

	_selectRegionGrower( InViewportClient, InViewport );

	if(!_regionGrower)
		return;
	auto meshInterface = _regionGrower->meshInterface();

	auto hit = meshInterface->getIntersectionAndFace( traceStart, traceDirection );

	if(!hit.first)
		return;

	FVector& hitPoint = hit.second.point;

    _clickPlane = FPlane( hitPoint, -traceDirection.GetSafeNormal() );

    _lineTraceStart = hitPoint;

    _isDragging = true;

    _regionGrower->startStretch( hitPoint );
}

bool FTCODSEditorMode::_strechUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
    FVector traceStart;
    FVector traceDirection;
    FVector traceEnd;

    getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

    // they will intersect, we're not parallel during the drag
    FVector point = FMath::LinePlaneIntersection( traceStart, traceEnd, _clickPlane );

    _regionGrower->updateStretch( point );

    return true;
}

bool FTCODSEditorMode::_strechEnd( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
    _isDragging = false;

    FVector traceStart;
    FVector traceDirection;
    FVector traceEnd;

    getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

    // they will intersect, we're not parallel during the drag
    FVector point = FMath::LinePlaneIntersection( traceStart, traceEnd, _clickPlane );

    _regionGrower->endStretch( point );


    return true; 
}


// -----------------------------------------------------------------------------
// Exemplar Painting
// -----------------------------------------------------------------------------

void FTCODSEditorMode::_sketchClick( FEditorViewportClient* InViewportClient, FViewport* InViewport)
{
	if(!_regionGrower)
		return;

	_sketchMode = SketchMode::Sketching;
    
    _regionGrower->startPaint();

	_selectionDidChange = false;
}

bool FTCODSEditorMode::_sketchUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	if(!_regionGrower || _sketchMode != SketchMode::Sketching )
		return false;

	if(InViewport->KeyState( EKeys::X ))
		return _eraseUpdate( InViewportClient, InViewport );
	else if(InViewport->KeyState( EKeys::Z ))
		return _sketchSelectUpdate( InViewportClient, InViewport );

    FVector traceStart;
    FVector traceDirection;
    FVector traceEnd;
    
    getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );
    
	auto meshInterface = _regionGrower->meshInterface();

	auto hit = meshInterface->getIntersectionAndFace( traceStart, traceDirection );

	if(!hit.first)
		return false;

	auto& surfacePoint = hit.second;

    _regionGrower->addBrushPoint(surfacePoint.point, surfacePoint.surfaceIndex);
    
    return true;
}

bool FTCODSEditorMode::_sketchSelectUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	if(!_regionGrower || _sketchMode != SketchMode::Sketching)
		return false;

	URegionGrowingComponent * grower = _hitRegionGrower( InViewportClient, InViewport );

	if(!grower)
		return true;

	FVector traceStart;
	FVector traceDirection;
	FVector traceEnd;

	getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

	auto meshInterface = grower->meshInterface();

	auto hit = meshInterface->getIntersectionAndFace( traceStart, traceDirection );

	if(!hit.first)
		return true;

	tcodsMeshInterface::SurfacePoint& surfacePoint = hit.second;

	UWorld * world = grower->GetWorld();

	TArray<FOverlapResult> overlaps;
	FCollisionObjectQueryParams queryParams(FCollisionObjectQueryParams::InitType::AllObjects);
	FCollisionShape collisionShape;
	collisionShape.SetSphere( _regionGrower->selectionRadius );

	world->OverlapMultiByObjectType( overlaps, surfacePoint.point, FQuat::Identity, queryParams, collisionShape );


	USelection * selection = Owner->GetSelectedActors();

	selection->BeginBatchSelectOperation();

	for(FOverlapResult& overlapResult : overlaps)
	{
		AActor * nearestActor = overlapResult.GetActor();

		if(nearestActor && nearestActor->IsA( AElementActor::StaticClass() ))
		{
			selection->Select( nearestActor, true );
		}

		_selectionDidChange = true;
	}

	selection->EndBatchSelectOperation( true );

	_selectionLastGrower = grower;


	return true;
}

void FTCODSEditorMode::_updateSelection()
{
	USelection * selection = Owner->GetSelectedActors();

	TArray<AElementActor*> elementActors;
	selection->GetSelectedObjects<AElementActor>( elementActors );

	std::vector<AElementActor*> actorsVector( elementActors.Num() );
	for(AElementActor * actor : elementActors)
		actorsVector.push_back( actor );

	_regionGrower->updateExampleSelection( actorsVector, 1.0f );
}

bool FTCODSEditorMode::_sketchEnd( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	if(!_regionGrower)
		return true;

    _regionGrower->endPaint();
    
	if(_selectionDidChange)
	{
		_updateSelection();
	}

    return true;
}



bool FTCODSEditorMode::_eraseUpdate( FEditorViewportClient* InViewportClient, FViewport* InViewport )
{
	if(!_regionGrower )
		return false;

	FVector traceStart;
	FVector traceDirection;
	FVector traceEnd;

	getRay( InViewportClient, InViewport, traceStart, traceDirection, traceEnd );

	auto meshInterface = _regionGrower->meshInterface();

	auto hit = meshInterface->getIntersectionAndFace( traceStart, traceDirection );

	if(!hit.first)
		return false;

	tcodsMeshInterface::SurfacePoint& surfacePoint = hit.second;

	_regionGrower->eraseAt( surfacePoint.point, surfacePoint.surfaceIndex );

	return true;
}



// -----------------------------------------------------------------------------
/// Tool
// -----------------------------------------------------------------------------

FTrivialConnectionsTool::FTrivialConnectionsTool()
{
    ID = MT_GeometryModify;
}

bool FTrivialConnectionsTool::StartModify()
{
    UE_LOG(LogTemp, Warning, TEXT("StartModify"));
    
    return false;
}

bool FTrivialConnectionsTool::EndModify()
{
    UE_LOG(LogTemp, Warning, TEXT("EndModify"));
    
    return true;
}

void FTrivialConnectionsTool::StartTrans()
{
    UE_LOG(LogTemp, Warning, TEXT("StartTrans"));
}

void FTrivialConnectionsTool::EndTrans()
{
    UE_LOG(LogTemp, Warning, TEXT("EndTrans"));
}

// -----------------------------------------------------------------------------
/// Toolkit
// -----------------------------------------------------------------------------

FTrivialConnectionsToolkit::FTrivialConnectionsToolkit(class FTCODSEditorMode* owningMode_in)
: editorMode(owningMode_in)
{
}

void FTrivialConnectionsToolkit::RegisterTabSpawners(const TSharedRef<class FTabManager>& TabManager)
{
    
}

void FTrivialConnectionsToolkit::UnregisterTabSpawners(const TSharedRef<class FTabManager>& TabManager)
{
    
}

/** Initializes the geometry mode toolkit */
void FTrivialConnectionsToolkit::Init(const TSharedPtr< class IToolkitHost >& InitToolkitHost)
{
    widget = SNew(STrivialConnectionWidget, SharedThis(this), editorMode);
    
    FModeToolkit::Init(InitToolkitHost);
}

/** IToolkit interface */
FName FTrivialConnectionsToolkit::GetToolkitFName() const
{
    return FName("TrivialConnectionsMode");
}

FText FTrivialConnectionsToolkit::GetBaseToolkitName() const
{
    return NSLOCTEXT("TCODS", "ToolkitName", "Trivial Connections");
}

FEdMode* FTrivialConnectionsToolkit::GetEditorMode() const
{
    return editorMode;
}

TSharedPtr<class SWidget> FTrivialConnectionsToolkit::GetInlineContent() const
{
    return widget;
}

// -----------------------------------------------------------------------------
/// Widget
// -----------------------------------------------------------------------------

STrivialConnectionWidget::~STrivialConnectionWidget()
{
    
}

void STrivialConnectionWidget::Construct(const FArguments& InArgs, TSharedRef<FTrivialConnectionsToolkit> InParentToolkit, class FTCODSEditorMode* owningMode_in)
{
    editorMode = owningMode_in;

    ChildSlot
    [
        SNew( SScrollBox )
        + SScrollBox::Slot()
        .Padding( 4.0f )
        [
            SNew(SVerticalBox)
            +SVerticalBox::Slot()
            .Padding(2.0f, 0.0f)
            .VAlign(VAlign_Top)
            [
                SNew(SHorizontalBox)
                +SHorizontalBox::Slot()
                .Padding(2.0f, 0.0f)
                .FillWidth(1)
                .HAlign(HAlign_Left)
                [
                    SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "Add", "Add"))
                ]
                +SHorizontalBox::Slot()
                .AutoWidth()
                .Padding(0.0f, 0.0f, 2.0f, 0.0f)
                .HAlign(HAlign_Right)
                [
                    SNew(SCheckBox)
                    .Style(FEditorStyle::Get(), "RadioButton")
                    .IsChecked( this, &STrivialConnectionWidget::IsRadioChecked, FTCODSEditorMode::ClickMode::AddSeed )
                    .OnCheckStateChanged_Lambda([this](ECheckBoxState newRadioState)->void
                    {
                        this->editorMode->setClickMode(FTCODSEditorMode::ClickMode::AddSeed);
                    })
                    .Content()
                    [
                        SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "Seeds", "Seeds"))
                    ]
                ]
                +SHorizontalBox::Slot()
                .AutoWidth()
                .Padding(0.0f, 0.0f, 10.0f, 0.0f)
                .HAlign(HAlign_Right)
                [
                    SNew(SCheckBox)
                    .Style(FEditorStyle::Get(), "RadioButton")
                    .IsChecked( this, &STrivialConnectionWidget::IsRadioChecked, FTCODSEditorMode::ClickMode::AddSingularity )
                    .OnCheckStateChanged_Lambda([this](ECheckBoxState newRadioState)->void
                    {
                        this->editorMode->setClickMode(FTCODSEditorMode::ClickMode::AddSingularity);
                    })
                    .Content()
                    [
                        SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "Singularities", "Singularities"))
                    ]
                ]
                +SHorizontalBox::Slot() 
                .AutoWidth()
                .Padding(0.0f, 0.0f, 10.0f, 0.0f)
                .HAlign(HAlign_Right)
                [
                    SNew(SCheckBox)
                    .Style(FEditorStyle::Get(), "RadioButton")
                    .IsChecked( this, &STrivialConnectionWidget::IsRadioChecked, FTCODSEditorMode::ClickMode::Stretch )
                    .OnCheckStateChanged_Lambda([this](ECheckBoxState newRadioState)->void
                    {
                        this->editorMode->setClickMode(FTCODSEditorMode::ClickMode::Stretch);
                    })
                    .Content()
                    [
                        SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "Strech", "Strech"))
                    ]
                ]
                +SHorizontalBox::Slot()
                .AutoWidth()
                .Padding(0.0f, 0.0f, 10.0f, 0.0f)
                .HAlign(HAlign_Right)
                [
                    SNew(SCheckBox)
                    .Style(FEditorStyle::Get(), "RadioButton")
                    .IsChecked( this, &STrivialConnectionWidget::IsRadioChecked, FTCODSEditorMode::ClickMode::SketchExemplar )
                    .OnCheckStateChanged_Lambda([this](ECheckBoxState newRadioState)->void
                    {
                        this->editorMode->setClickMode(FTCODSEditorMode::ClickMode::SketchExemplar);
                    })
                    .Content()
                    [
                        SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "Sketch", "Sketch"))
                    ]
                ]
            ]
            +SVerticalBox::Slot()
            .Padding(2.0f, 0.0f)
            .VAlign(VAlign_Top)
            [
                SNew(STextBlock).Text(NSLOCTEXT("STrivialConnectionWidget", "To Remove, Shift-Click", "To Remove, Shift-Click"))
            ]
			+ SVerticalBox::Slot()
			.Padding( 0.0f, 10.0f )
			.VAlign( VAlign_Top )
			[
				SNew( SButton ).Text( NSLOCTEXT( "STrivialConnectionWidget", "To Exemplar", "To Exemplar" ) )
				.OnClicked_Lambda( [this]()->FReply {
					this->editorMode->toExemplar();
					return FReply::Handled();
				})
			]
            + SVerticalBox::Slot()
            .Padding( 0.0f, 10.0f )
            .VAlign( VAlign_Top )
            [
                SNew( SButton ).Text( NSLOCTEXT( "STrivialConnectionWidget", "Convert to Elements", "Convert to Elements" ) )
                .OnClicked_Lambda( [this]()->FReply {
                    this->editorMode->convertToElements();
                    return FReply::Handled();
                })
            ]
            + SVerticalBox::Slot()
            .Padding( 0.0f, 10.0f )
            .VAlign( VAlign_Top )
            [
                SNew( SButton ).Text( NSLOCTEXT( "STrivialConnectionWidget", "Clear", "Clear" ) )
                .OnClicked_Lambda( [this]()->FReply {
                    this->editorMode->clear();
                    return FReply::Handled();
                })
            ]
         ]
     ];


//    SNew(STextBlock).Text(NSLOCTEXT("TCODS", "Whatever", "Whatever"))
//    ]
//    ];
    
}


#endif // WITH_EDITOR




