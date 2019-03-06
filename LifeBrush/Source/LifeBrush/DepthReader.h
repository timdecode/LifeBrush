// Copyright 2015, Timothy Davison. All rights reserved.

#pragma once

#include "Components/ActorComponent.h"

#include "Utility.h"

#include "DepthReader.generated.h"


UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UDepthReader : public UActorComponent
{
	GENERATED_BODY()

public:	
	// Sets default values for this component's properties
	UDepthReader();

	// Called when the game starts
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction ) override;

    void updateBuffer();
    
    // $x \in [0,1], y \in [0,1]$
    float depthAt(float x, float y);
    float depthAt(const FVector& worldPosition);
    
    float depthOf(const FVector& worldPosition);
    
    bool inFrustum(const FVector& worldPosition);
public:
    UPROPERTY(Category = RegionGrowing, EditAnywhere)
    UTextureRenderTarget2D * renderTarget;
    
    UPROPERTY(Category = RegionGrowing, EditAnywhere)
    bool trackSceneView = true;
    
    UPROPERTY(Category = RegionGrowing, EditAnywhere)
    int32 maxSize = 200; // maximum width or height of the texture
    
    USceneCaptureComponent2D * sceneCaptureComponent();
private:
	bool _didUpdateThisFrame = false;
    
    TArray<FColor> buffer;
    SceneViewAndFamily _sceneViewAndFamily;
    
    float _surfaceWidth;
    float _surfaceHeight;
};
