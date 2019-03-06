// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"

#include "GraphSnapshot.generated.h"

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FGraphSnapshot
{
	GENERATED_BODY()

public:
	void snapshot(FGraph& graph);
	void restore(FGraph& graph);
public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	FGraph graph;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	int32 tickCount = 0;
};
