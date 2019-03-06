// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "MeshInterfaceMode.generated.h"

UENUM(BlueprintType)
enum class EMeshInterfaceMode : uint8
{
	StaticMesh UMETA(DisplayName = "UStaticMeshComponent"),
	ChunkedMesh UMETA(DisplayName = "Chunked Mesh"),
};