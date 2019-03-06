// Fill out your copyright notice in the Description page of Project Settings.

#include "LifeBrush.h"
#include "CurveFittingComponent.h"
#include "ElementActor.h"
#include "ContextActor.h"
#include "Utility.h"
#include "tcodsMeshInterface.h"

#include <vector>
#include <map>
#include <cmath>
#include <unordered_map>
#include <ctime>

#include "Algorithm/Domain.h"
#include "Algorithm/Entity.hpp"

// unreal includes
#include "StaticMeshResources.h"
#include "RuntimeMeshComponent.h"
#include "RuntimeMeshLibrary.h"
//#include "InstancedStaticMeshComponent.h"


#include "Algorithm/tcods/MeshIO.h"
#include "Algorithm/tcods/Problem.h"

#if WITH_EDITOR
#include "LevelEditor.h"
#endif

// fuck windows http://stackoverflow.com/questions/118774/is-there-a-clean-way-to-prevent-windows-h-from-creating-a-near-far-macro
#undef near

void UCurveFittingComponent::LoadExemplar()
{

}

void UCurveFittingComponent::DoStep()
{

}
