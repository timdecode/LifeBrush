//
//  Created by Timothy Davison on 2018-12-28.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once

#include "RegionGrowingGeneratorTool.h"

#include "ElementEditor/StringGenerator.h"
#include "ElementEditor/SwarmGenerator.h"

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




