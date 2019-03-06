//
//  Entity.hpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2016-01-22.
//  Copyright © 2016 Timothy Davison, Inc. All rights reserved.
//

#pragma once

#include "Eigen/Dense"

#include <stdint.h>
#include <vector>

struct Entity
{
public:
    std::vector<int32_t> elementIndices;
    int8_t entityID;
};