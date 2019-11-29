// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/ObjectSimulation.h"

#include "Brownian.generated.h"

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSingleParticleBrownian : public FGraphObject
{
	GENERATED_BODY()

public:
	float time = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float dampening = 0.0f;
};

UCLASS(BlueprintType)
class LIFEBRUSH_API USingleParticleBrownianSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float minSpeed = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float maxSpeed = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float minTime = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float maxTime = 0.0f;

protected:
	FRandomStream rand;

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void tick(float deltaT) override;
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UGlobalParticleBrownianSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float minSpeed = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float maxSpeed = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float minTime = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	float maxTime = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite)
	bool enabled = false;


protected:
	FRandomStream rand;

	TArray<float> _times;

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void tick(float deltaT) override;
};