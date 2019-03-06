//
//  LifeBrush
//
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "Algorithm/SynthesisContext.h"
#include "ShipEditorSimulation/ObjectSimulation.h"

#include "ElementGenerator.generated.h"

struct FFlexSimulation;

UCLASS(ClassGroup = (Custom), EditInlineNew, DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UElementGenerator : public UObject
{
	GENERATED_BODY()

public:
	virtual ~UElementGenerator() {}

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
	{

	}

	virtual void detach()
	{

	}

	virtual void tick(float deltaT) {}

	virtual bool wantsFlex() { return false; }

	virtual void start() {}
	virtual void stop() {}

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void endBrushPath() {}

	virtual void beginEraserPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void eraseInRadiusAt(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void endEraserPath() {}

	// Only set if wantsFlex is true.
	FFlexSimulation * flexSimulation = nullptr;
};