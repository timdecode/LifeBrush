//
//  Created by Timothy Davison on 2018-12-28.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "ElementGenerator.h"

#include "Algorithm/Algorithm.h"

#include "StringGenerator.generated.h"

UCLASS(ClassGroup = (Custom), DefaultToInstanced, meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UStringGenerator : public UElementGenerator
{
	GENERATED_BODY()

public:
	SynthesisContext * _context;

	UPROPERTY()
	UGraphSimulationManager * _simulationManager;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float scaleFactor = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * material = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * staticMesh = nullptr;

protected:
	std::vector< PositionRadiusFace > _brushPoints;

	// local distance along the vector between the current and next point
	float _t = 0.0f;

	FGraphNodeHandle _lastStringElement;

public:
	virtual ~UStringGenerator() {}

	virtual void attach(SynthesisContext * context, UGraphSimulationManager * simulationManager) override;

	virtual void detach() override;

	virtual void tick(float deltaT) override;

	virtual bool wantsFlex()  override { return true; }

	virtual void beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface) override;
	virtual void endBrushPath() override;

protected:
	void _initPath();
	void _tick(float deltaT);
};

