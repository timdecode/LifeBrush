//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#pragma once

#include "RegionGrowingGeneratorTool.h"

#include "ElementEditor/StringGenerator.h"

#include "FilamentGeneratorTool.generated.h"

class AFilamentPrototypeActor;
class UFilamentGenerator;

UCLASS(Blueprintable)
class UFilamentGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~UFilamentGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void faceDown_released() override;


protected:
	TArray<AFilamentPrototypeActor*> overlappingPrototypes(UPrimitiveComponent * hand, float radius);

	UFilamentGenerator * filamentGenerator();

	void _tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform);

protected:
	enum class Mode {
		Selecting,
		Painting
	};

	Mode _mode = Mode::Selecting;

	AFilamentPrototypeActor * _selection = nullptr;
};
