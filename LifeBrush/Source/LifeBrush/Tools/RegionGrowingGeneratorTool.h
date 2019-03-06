// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#include "RegionGrowingTool.h"

#include "ElementEditor/DiscreteElementEditorComponent.h"
#include "ElementEditor/RegionGrowingGenerator.h"

#include "RegionGrowingGeneratorTool.generated.h"



UCLASS(Blueprintable)
class UGeneratorTool : public UGenerativeBrushTool
{
	GENERATED_BODY()

public:
	UDiscreteElementEditorComponent * elementEditor = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	AActor * exemplarActor = nullptr;

	UPROPERTY()
	UElementGenerator * _generator = nullptr;

public:
	virtual ~UGeneratorTool();

	void init(FRGC_UToolInitProperties& initProperties);

	virtual void focused() override;
	virtual void loseFocus() override;

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void twoHandStart(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override {}
	virtual void twoHandEnd(UPrimitiveComponent * handA, UPrimitiveComponent * handB) override {}

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void tickTwoHand
	(
		float dt,
		UPrimitiveComponent * handA,
		UPrimitiveComponent * handB,
		FTransform transformA,
		FTransform transformB
	) override {}

protected:
	virtual void _tickOneHand_generate(float dt, UPrimitiveComponent * hand, FTransform lastTransform);
	virtual void _tickOneHand_erase(float dt, UPrimitiveComponent * hand, FTransform lastTransform);
};

UCLASS(Blueprintable)
class URegionGrowingGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~URegionGrowingGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);

	URegionGrowingGenerator * generator() { return Cast<URegionGrowingGenerator>(_generator); }

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

protected:
	virtual void _tickOneHand_generate(float dt, UPrimitiveComponent * hand, FTransform lastTransform);
	virtual void _tickOneHand_erase(float dt, UPrimitiveComponent * hand, FTransform lastTransform);
};