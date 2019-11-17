//
//  Created by Timothy Davison on 2019-10-25.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#pragma once

#include "RegionGrowingGeneratorTool.h"

#include "SwarmGeneratorTool.generated.h"

class USwarmGenerator;
class ASGPrototypeActor;

UCLASS(Blueprintable)
class USwarmGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~USwarmGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);

	USwarmGenerator * swarmGenerator();

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

protected:
	TArray<ASGPrototypeActor*> _overlappingPrototypes(UPrimitiveComponent * hand, float radius);

	void _tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform);

protected:
	enum class Mode {
		Selecting,
		Painting
	};

	Mode _mode = Mode::Selecting;

	ASGPrototypeActor * _selection = nullptr;
};