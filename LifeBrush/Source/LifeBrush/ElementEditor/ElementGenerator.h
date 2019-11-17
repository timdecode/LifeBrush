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
class UFlexSimulationComponent;
class UDiscreteElementEditorComponent;

UCLASS(ClassGroup = (Custom), EditInlineNew, DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UElementGenerator : public UObject
{
	GENERATED_BODY()

public:
	virtual ~UElementGenerator() {}

	virtual void init(SynthesisContext * context, UGraphSimulationManager * simulationManager) {}

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) {}
	virtual void detach() {}

	virtual void tick(float deltaT) {}
	virtual void tickPaused(float deltaT) {}

	virtual bool wantsFlex() { return false; }

	virtual void start() {}
	virtual void stop() {}

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void endBrushPath() {}

	virtual void beginEraserPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void eraseInRadiusAt(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) {}
	virtual void endEraserPath() {}

	// these generators must receive an init before we do
	virtual std::vector<UClass*> dependencies() { return std::vector<UClass*>(); }

	FFlexSimulation * flexSimulation = nullptr;
	UFlexSimulationComponent * flexSimulationComponent = nullptr;
	UDiscreteElementEditorComponent * elementEditorComponent = nullptr;
	AActor * exemplarActor = nullptr;
};