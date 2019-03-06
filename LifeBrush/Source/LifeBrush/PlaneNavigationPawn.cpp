// Copyright 2015, Timothy Davison. All rights reserved.

#include "LifeBrush.h"
#include "PlaneNavigationPawn.h"

FName APlaneNavigationPawn::movementComponentName(TEXT("MovementComponent0"));
FName APlaneNavigationPawn::collisionComponentName(TEXT("CollisionComponent0"));

// Sets default values
APlaneNavigationPawn::APlaneNavigationPawn(const FObjectInitializer& ObjectInitializer)
: Super(ObjectInitializer)
{
    PrimaryActorTick.bCanEverTick = true;
    
    bCanBeDamaged = true;
    
    bReplicates = true;
    NetPriority = 3.0f;
    
    BaseEyeHeight = 0.0f;
    bCollideWhenPlacing = false;
    SpawnCollisionHandlingMethod = ESpawnActorCollisionHandlingMethod::AlwaysSpawn;
    
//    MovementComponent = CreateDefaultSubobject<UPaw/*nMovementComponent, UFloatingPawnMovement>(APlaneNavigationPawn::movementComponentName);
    
    RootComponent = CreateDefaultSubobject<USceneComponent>(TEXT("rootComponent0"));
}

// Called when the game starts or when spawned
void APlaneNavigationPawn::BeginPlay()
{
	Super::BeginPlay();
    
    auto controller = GetWorld()->GetFirstPlayerController();
    
    controller->bShowMouseCursor = true;
    controller->bEnableClickEvents = true;
    controller->bEnableMouseOverEvents = true;
    controller->SetIgnoreLookInput(true);
}

// Called every frame
void APlaneNavigationPawn::Tick( float DeltaTime )
{
	Super::Tick( DeltaTime );
    
    if( dragActive )
        offsetByMousePosition();
}

// Called to bind functionality to input
void APlaneNavigationPawn::SetupPlayerInputComponent(class UInputComponent* InputComponent)
{
	Super::SetupPlayerInputComponent(InputComponent);
    
    InputComponent->BindAxis("LookX", this, &APlaneNavigationPawn::lookX);
    InputComponent->BindAxis("LookY", this, &APlaneNavigationPawn::lookY);
    InputComponent->BindAxis("Zoom", this, &APlaneNavigationPawn::zoom);

    InputComponent->BindAction("LeftClick", IE_Pressed, this, &APlaneNavigationPawn::setAnchor);
    InputComponent->BindAction("LeftClick", IE_Released, this, &APlaneNavigationPawn::unsetAnchor);
}

void APlaneNavigationPawn::setAnchor()
{
    APlayerController * playerController = GetWorld()->GetFirstPlayerController();

    FVector2D mousePosition;
    playerController->GetMousePosition(mousePosition.X, mousePosition.Y);
    
    FHitResult hitResult;
    const bool traceComplex = false;
    if( playerController->GetHitResultAtScreenPosition(mousePosition, ECC_Visibility, traceComplex, hitResult) )
    {
        FVector mOriginA;
        FVector mDirectionA;
        UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition, mOriginA, mDirectionA);
        
        FVector hitA = hitResult.ImpactPoint;
        FPlane plane = FPlane(hitA, -mDirectionA);

        FVector mOriginB;
        FVector mDirectionB;
        UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition + FVector2D(1.0f, 0.0f), mOriginB, mDirectionB);
        
        FVector hitB = FMath::LinePlaneIntersection(mOriginB, mOriginB + mDirectionB, plane);

        mouseDelta = FVector::Dist(hitA, hitB);
        
        dragActive = true;
        lastMouse = mousePosition;
    }
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    
//    APlayerController * playerController = GetWorld()->GetFirstPlayerController();
//
//    FVector2D mousePosition;
//    playerController->GetMousePosition(mousePosition.X, mousePosition.Y);
//    
//    FHitResult hitResult;
//    const bool bTraceComplex = false;
//    if( playerController->GetHitResultAtScreenPosition(mousePosition, ECC_Visibility, bTraceComplex, hitResult) )
//    {
//        FVector hitPoint = hitResult.ImpactPoint;
//        FVector normal = {0.0f, 1.0f, 0.0f};
//        
//        anchorPlane = FPlane(hitPoint, normal);
//        lastPoint = hitPoint;
//        
//        
//        UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition, origin, direction)
//        
//        FVector2D m2 = mousePosition + FVector2D(1.0f, 0.0f);
//        FVector hit2 = FMath::LinePlaneIntersection(hitPoint, hitPoint + normal, anchorPlane);================789///////////////988=
//        
//        mouseDelta = FVector::Dist(hit2 - hitResult.hitPoint);
//    }
//    
//    
//    lastMouse = mousePosition;
//    dragActive = true;
}

void APlaneNavigationPawn::unsetAnchor()
{
    dragActive = false;
}

void APlaneNavigationPawn::offsetByMousePosition()
{
    APlayerController * playerController = GetWorld()->GetFirstPlayerController();
    
    FVector2D mousePosition;
    playerController->GetMousePosition(mousePosition.X, mousePosition.Y);

    FVector2D delta = mousePosition - lastMouse;
    
    
    FVector offset = {delta.X, 0.0f, delta.Y};
    offset *= mouseDelta;
    
    RootComponent->AddWorldOffset(offset);
    
    lastMouse = mousePosition;
    
    
    
//    FVector mOriginA;
//    FVector mDirectionA;
//    UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition, mOriginA, mDirectionA);
//    
//    FVector hitA = hitPoint.ImpactPoint;
//    FPlane plane = FPlane(hitA, -mDirection);
//    
//    FVector mOriginB;
//    FVector mDirectionB;
//    UGameplayStatics::DeprojectScreenToWorld(playerController, lastMouse, mOriginB, mDirectionB);
//    
//    
//    mouseDeleta = FVector::Dist(hitA, hitB);
    
    
//    APlayerController * playerController = GetWorld()->GetFirstPlayerController();
//    
//    FVector2D mousePosition;
//    playerController->GetMousePosition(mousePosition.X, mousePosition.Y);
//
//    if( mousePosition == lastMouse )
//        return;
//    
//    FVector m0 = {lastMouse.X, 0.0f, lastMouse.Y};
//    FVector m1 = {mousePosition.X, 0.0f, mousePosition.Y};
//    
//    FVector offset = m0 - m1;
//    
//    RootComponent->AddWorldOffset(-offset * mouseDelta);
//    
//    lastMouse = mousePosition;
    
//    FVector origin;
//    FVector direction;
//    if( UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition, origin, direction) )
//    {
//        FVector hit = FMath::LinePlaneIntersection(origin, origin + direction, anchorPlane);
//        
//        FVector offset = hit - lastPoint;
//        
//        RootComponent->AddWorldOffset(-offset);
//        
//        // deproject one more time with our updated position
//        UGameplayStatics::DeprojectScreenToWorld(playerController, mousePosition, origin, direction);
//        hit = FMath::LinePlaneIntersection(origin, origin + direction, anchorPlane);
//        
//        lastPoint = hit;
//        lastMouse = mousePosition;
//    }
}

void APlaneNavigationPawn::lookX(float dx)
{
}

void APlaneNavigationPawn::lookY(float dy)
{


}

void APlaneNavigationPawn::zoom(float dz)
{
    if( dz == 0.0f )
        return;
    
    FVector offset = {0.0f, dz, 0.0f};
    
    AddMovementInput( offset );
}