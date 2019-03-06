//  RegionGrowing
//
//  Created by Timothy Davison on 2018-11-28.
//  Copyright (c) Timothy Davison. All rights reserved.
//

#pragma once


#include "VRTool.h"
#include "RegionGrowingComponent.h"

struct FRGC_UToolInitProperties : public FUToolInitProperties
{
public:
	URegionGrowingComponent * regionGrowingComponent;
};
