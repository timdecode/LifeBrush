// Copyright 2015, Timothy Davison. All rights reserved.

#pragma once

#include "GameFramework/Pawn.h"
#include "PlaneNavigationPawn.generated.h"

UCLASS()
class LIFEBRUSH_API APlaneNavigationPawn : public APawn
{
	GENERATED_BODY()

public:
	// Sets default values for this pawn's properties
	APlaneNavigationPawn(const FObjectInitializer& ObjectInitializer);

	// Called when the game starts or when spawned
	virtual void BeginPlay() override;
	
	// Called every frame
	virtual void Tick( float DeltaSeconds ) override;

	// Called to bind functionality to input
	virtual void SetupPlayerInputComponent(class UInputComponent* InputComponent) override;

public:
    /** DefaultPawn movement component */
    UPROPERTY(Category = Pawn, VisibleAnywhere, BlueprintReadOnly, meta = (AllowPrivateAccess = "true"))
    UPawnMovementComponent* MovementComponent;
    
private:
    void setAnchor();
    void unsetAnchor();
    
    void lookX(float dx);
    void lookY(float dy);
    void zoom(float dz);
    
    void offsetByMousePosition();
    
    
    static FName movementComponentName;
    static FName collisionComponentName;
    
    FPlane anchorPlane;
    FVector lastPoint;
    
    FVector2D lastMouse;
    
    float mouseDelta = 0.0f;
    
    bool dragActive = false;
};
