//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#pragma once

#include "Algorithm/PointCloud.hpp"
#include <vector>
#include <memory>
#include "Algorithm/SurfaceIndex.h"

struct LinearPath
{
public:
	// the path
	SurfacePositionCloud _brushPoints;
	std::vector<FQuat> _segmentRotations;
	std::unique_ptr<SurfacePositionDynamicCloud> _brushIndex;

	FVector _binormal;
	FVector _normal;

	float _runningLength = 0.0f;

	// The length between path segments
	float segmentLength = 0.5f; 

public:
	struct PositionRotation
	{
		FVector position = FVector::ZeroVector;
		FQuat rotation = FQuat::Identity;
		uint32_t segment = 0;
	};

public:
	LinearPath() { clear(); }

	// Adds a point to the path if it is segmentLength from the last point in the path.
	// \return Whether the point was added.
	bool conditionallyAddPoint(FVector point, float radius, FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface);

	void clear();

	float length();

	auto frameAlongPath(const FVector startPointOnPath, const uint32_t startSegment, const float distance, PositionRotation& resultFrame) -> bool;

	auto radiusSearch(FVector point, float radius) -> std::vector<std::pair<size_t, float>>;

	auto brushPoints()->std::vector<SurfacePosition>& { return _brushPoints.pts; }

	auto segmentRotation(uint32_t segment)->FQuat;

	// This is O(n), but we could make it O(log n) with some work.
	auto expensiveNearestPointOnPath(const FVector point, FVector& nearest, uint32_t& segment)->bool;
};