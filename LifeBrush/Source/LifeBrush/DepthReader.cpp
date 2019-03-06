// Copyright 2015, Timothy Davison. All rights reserved.

#include "LifeBrush.h"
#include "DepthReader.h"


// Sets default values for this component's properties
UDepthReader::UDepthReader()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

//    bWantsInitializeComponent = true;
    bIsActive = true;
    bTickInEditor = true;
    
	// ...
}


// Called when the game starts
void UDepthReader::BeginPlay()
{
	Super::BeginPlay();
}


// Called every frame
void UDepthReader::TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction )
{
	Super::TickComponent( DeltaTime, TickType, ThisTickFunction );
    
	_didUpdateThisFrame = false;
}

USceneCaptureComponent2D* UDepthReader::sceneCaptureComponent()
{
    AActor * actor = this->GetOwner();

    if( !actor )
        return nullptr;
    
    USceneCaptureComponent2D * capture = (USceneCaptureComponent2D*)actor->GetComponentByClass(USceneCaptureComponent2D::StaticClass());
    
    if( capture == nullptr && actor )
    {
        capture = NewObject<USceneCaptureComponent2D>(actor);
        capture->SetMobility(EComponentMobility::Movable);

        USceneComponent * scene = actor->GetRootComponent();
        
        if( scene )
            capture->AttachToComponent(scene, FAttachmentTransformRules::KeepRelativeTransform);
        else
            actor->SetRootComponent(capture);
        
        actor->AddInstanceComponent(capture);
    }
    
    return capture;
}

void UDepthReader::updateBuffer()
{
	USceneCaptureComponent2D * capture = sceneCaptureComponent();

	if(_didUpdateThisFrame || !renderTarget || !capture )
		return;

    buffer.Reset();

	_didUpdateThisFrame = true;

	_sceneViewAndFamily.init( GetWorld() );
    
    int32 sizeX = _sceneViewAndFamily.viewRect.Width();
    int32 sizeY = _sceneViewAndFamily.viewRect.Height();
    
    
//    if( sizeX > sizeY && sizeX > maxSize )
//    {
//        float aspectRatio = float(sizeY) / float(sizeX);
//    
//        sizeX = maxSize;
//        sizeY = maxSize * aspectRatio;
//    }
//    else if( sizeY > sizeX && sizeY > maxSize )
//    {
//        float aspectRatio = float(sizeX) / float(sizeY);
//        
//        sizeY = maxSize;
//        sizeX = maxSize * aspectRatio;
//    }
//    
    
    
	if(sizeX != renderTarget->SizeX)
		renderTarget->SizeX = sizeX;

	if(sizeY != renderTarget->SizeY)
		renderTarget->SizeY = sizeY;

	if(trackSceneView)
	{
		AActor * actor = GetOwner();

		FVector location;
		FRotator rotation;

		SceneViewAndFamily::viewLocationRotation( GetWorld(), location, rotation );

		actor->SetActorLocation( location );
		actor->SetActorRotation( rotation );
	}

	if(capture->TextureTarget == nullptr)
		capture->TextureTarget = renderTarget; 

	_surfaceWidth = renderTarget->GetSurfaceWidth();
	_surfaceHeight = renderTarget->GetSurfaceHeight();

    
    FTextureRenderTarget2DResource* textureResource = (FTextureRenderTarget2DResource*)renderTarget->Resource;
    

    textureResource->ReadPixels(buffer);
    
    // hmm, check this out
    //https://github.com/tk-master/VaOcean/blob/4.10/Source/VaOceanPlugin/Private/VaOceanSimulatorComponent.cpp
    // https://answers.unrealengine.com/questions/228348/how-to-normalise-the-scenedepth.html
    

    
    

}

float UDepthReader::depthAt(float x, float y)
{
    return 0;
}

float UDepthReader::depthOf(const FVector& worldPosition)
{
    FVector4 depth = _sceneViewAndFamily.viewProjectionMatrix.TransformFVector4(worldPosition);
    
    
    return depth.W;
}

float UDepthReader::depthAt(const FVector& worldPosition)
{
    FVector2D screenPosition;
    
    if( !_sceneViewAndFamily.didInit )
        return 0.0f;

    // will be $\in [0,1]$
    screenPosition = _sceneViewAndFamily.worldToNormalizedPoint(worldPosition);
    
    int xi = screenPosition.X * _surfaceWidth;
    int yi = (screenPosition.Y) * _surfaceHeight;
    
    int i = xi + yi * _surfaceWidth;
    
    if( i < 0 || i >= buffer.Num() )
        return 0.0f;
     
    FColor color = buffer[i];
    
	// http://aras-p.info/blog/2009/07/30/encoding-floats-to-rgba-the-final/
	FVector rgb( color.R, color.G, color.B );
	FVector divisor( 1.0/255.0, 1.0 / (255.0 * 255.0), 1.0 / (255.0 * 255.0 * 255.0) );

	float depth = FVector::DotProduct( rgb, divisor );
    
	// 120000 is from the shader to narrow precision on stuff close to the camera
    return depth * 12000.0f;
}

bool UDepthReader::inFrustum(const FVector& worldPosition)
{
    bool visible =  _sceneViewAndFamily.viewFrustum.IntersectSphere(worldPosition, 0.0f);
    
    return visible;
}



