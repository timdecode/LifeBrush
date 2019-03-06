// Copyright 2015, Timothy Davison. All rights reserved.

#include "LifeBrush.h"
#include "ContextActor.h"


// Sets default values
AContextActor::AContextActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = false;

}

// Called when the game starts or when spawned
void AContextActor::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AContextActor::Tick( float DeltaTime )
{
	Super::Tick( DeltaTime );

}

