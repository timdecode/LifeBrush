// Copyright (c) 2017 Timothy Davison. All rights reserved.

#pragma once

#include "Components/ActorComponent.h"

#include "Algorithm/ElementKDTree.h"

#include "Algorithm/ctpl_stl.h"

#include <map>
#include <memory>
#include <ctime>

#include "tcodsMeshInterface.h"

#include "CurveFittingComponent.generated.h"

class URuntimeMeshComponent;
class AElementActor;



UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UCurveFittingComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY() AActor * exemplar = nullptr;

	UFUNCTION( BlueprintCallable, Category = Generation ) void LoadExemplar();
	UFUNCTION( BlueprintCallable, Category = Generation ) void DoStep();

};
