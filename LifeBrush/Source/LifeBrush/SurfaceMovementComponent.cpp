// Copyright 2015, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "tcodsMeshInterface.h"

#include "SurfaceMovementComponent.h"


// Sets default values for this component's properties
USurfaceMovementComponent::USurfaceMovementComponent()
{
	// Set this component to be initialized when the game starts, and to be ticked every frame.  You can turn these features
	// off to improve performance if you don't need them.
	PrimaryComponentTick.bCanEverTick = true;

	// ...
}


// Called when the game starts
void USurfaceMovementComponent::BeginPlay()
{
	Super::BeginPlay();

	// ...
    auto root = Cast<UPrimitiveComponent>(GetOwner()->GetRootComponent());

    FRotator rot = GetOwner()->GetActorRotation();
    hackStartRotation = rot.Quaternion();
}


// Called every frame
void USurfaceMovementComponent::TickComponent( float DeltaTime, ELevelTick TickType, FActorComponentTickFunction* ThisTickFunction )
{
	Super::TickComponent( DeltaTime, TickType, ThisTickFunction );

	// ...
    URegionGrowingComponent * regionGrower = _regionGrowingComponent();
    
    if( !regionGrower )
        return;
    
    if( !regionGrower->isMeshInterfaceReady() )
        return;
    
    auto meshInterface = regionGrower->meshInterface();
    
    FVector position = GetOwner()->GetActorLocation();
    
    tcodsMeshInterface::SurfacePoint nearest = meshInterface->nearestPointOnMesh(position);
    auto rotationAndNormal = meshInterface->rotationAndNormalAtIndex(nearest.surfaceIndex);
    auto frame = meshInterface->frameAtNearest(nearest);
    
    auto root = Cast<UPrimitiveComponent>(GetOwner()->GetRootComponent());
    
    if( !root )
        return;
    
    // Linear acceleration
    FVector normal = rotationAndNormal.second;
    
    FVector offset = nearest.point - position;
    FVector velocity = root->GetPhysicsLinearVelocity(TEXT("None"));
    
    FVector offset_ = offset.ProjectOnTo(normal);
    FVector velocity_ = velocity.ProjectOnTo(normal);
    
    FVector force = offset_ * 4.0f - velocity_ * 2.0f;
    
    root->AddForce(force, TEXT("None"), true);
    
    // Angular acceleration
//    FVector up = root->GetUpVector();
//    FVector angularOffset = -normal - up;
//    
//    root->AddTorque(-angularOffset, TEXT("None"), true);
    
    FQuat rotation = rotationAndNormal.first * hackStartRotation;
    
    rotation = rotation.Inverse();
    root->SetWorldRotation(rotation);
    
    // random walk
    t -= DeltaTime;
    
    if( t < 0.0f )
    {
        r1 = FMath::FRandRange(-maxAcceleration, maxAcceleration);
        r2 = FMath::FRandRange(-maxAcceleration, maxAcceleration);
        
        t = FMath::FRandRange(0.0f, maxWalkTime);
    }
    
    FVector e1 = frame.first;
    FVector e2 = frame.second;
    
    force = e1 * r1 + e2 * r2;
    
    root->AddForce(force, TEXT("None"), true);
    
    // align to surface
}

