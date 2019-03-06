//
//  Created by Timothy Davison on 2018-12-28.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "RegionGrowingGeneratorTool.h"

#include "ElementEditor/StringGenerator.h"
#include "ElementEditor/CollagenGenerator.h"

#include "StringGeneratorTool.generated.h"

class USwarmGenerator;

UCLASS(Blueprintable)
class UStringGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~UStringGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);
};

UCLASS(Blueprintable)
class UCollagenGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~UCollagenGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);
};




UCLASS(Blueprintable)
class USwarmGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

public:
	virtual ~USwarmGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);

	USwarmGenerator * generator();

	virtual void faceUp_released(USceneComponent * interactionPoint /* = nullptr */) override;
};