//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "LinearPath.h"

bool LinearPath::conditionallyAddPoint(FVector point, float radius, FSurfaceIndex surfaceIndex)
{
	if (!_brushIndex)
		clear();

	auto& points = _brushPoints.pts;

	// check if we can create the segment yet (when we are segmentLength away from the last one)
	if (!points.empty())
	{
		FVector last = points.back().position;

		if ((last - point).Size() < segmentLength)
			return false;
	}

	// update the running length
	if (points.size() >= 1)
	{
		_runningLength += (points.back().position - point).Size();
	}

	// remember the point
	points.emplace_back(point, radius, surfaceIndex);
	auto index = _brushPoints.kdtree_get_point_count() - 1;
	_brushIndex->addPoints(index, index);



	// create our starting element
	if (points.size() == 2)
	{
		// we only create the starting element once we have a segment, because then we can compute an orientation
		// Forward is (1,0,0)
		FVector tangent = (points[1].position - points[0].position).GetSafeNormal();

		// add the first segment rotation
		_normal = FVector::CrossProduct(FVector::ForwardVector, tangent).GetSafeNormal();
		_binormal = FVector::CrossProduct(tangent, _normal).GetSafeNormal();

		FQuat segmentRotation = FQuat::FindBetween(FVector::ForwardVector, tangent);

		FMatrix matrix = FRotationMatrix::MakeFromXY(tangent, FVector::ForwardVector);

		segmentRotation = matrix.ToQuat();

		_segmentRotations.push_back(segmentRotation);
	}
	// add a new segment
	else if (points.size() > 2)
	{
		// update our segment rotations
		auto last = points.end();

		auto n = points.size();

		FVector lastTangent = (points[n - 2].position - points[n - 3].position).GetSafeNormal();

		FVector curTangent = (points[n - 1].position - points[n - 2].position).GetSafeNormal();

		if (curTangent.IsZero())
			curTangent = lastTangent;

		// Frenet-Serret frames: https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas
		// https://github.com/valentingalea/android-3d-engine/blob/master/Source/VTrajectory.cpp
		_normal = FVector::CrossProduct(_binormal, curTangent).GetSafeNormal();
		_binormal = FVector::CrossProduct(curTangent, _normal).GetSafeNormal();

		FMatrix matrix = FRotationMatrix::MakeFromXY(curTangent, _binormal);

		FQuat segmentRotation = matrix.ToQuat();

		_segmentRotations.push_back(segmentRotation);
	}

	return true;
}

void LinearPath::clear()
{
	_runningLength = 0.0f;

	_brushPoints.pts.clear();
	_brushIndex.reset(new SurfacePositionDynamicCloud(3, _brushPoints, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */)));

	_segmentRotations.clear();
}

float LinearPath::length()
{
	return _runningLength;
}

auto LinearPath::frameAlongPath(const FVector startPoint, const uint32_t startSegment, const float distance, PositionRotation& resultFrame) -> bool
{
	if (_brushPoints.pts.size() < 2)
		return false;

	const int increment = distance >= 0.0f || distance == -0.0f ? 1 : -1;

	auto& points = _brushPoints.pts;

	FVector frame_p = startPoint;

	float d = FMath::Abs(distance);

	for (int segment = startSegment; increment == 1 ? segment + 2 < points.size() : segment - 1 >= 0; segment += increment)
	{
		const int i_a = increment == 1 ? segment : segment + 1;
		const int i_b = increment == 1 ? segment + 1 : segment;

		const FVector p_a = points[i_a].position;
		const FVector p_b = points[i_b].position;

		const float segmentLength = FVector::Dist(p_b, p_a);

		const float l = FVector::Dist(p_b, frame_p);
		const float t = (segmentLength - l) / segmentLength;

		// continue if we have distance to travel and we have a next segment
		if (l < d)
		{
			frame_p = p_b;

			d -= l;
		}
		// we finished walking
		else
		{
			FVector dir = (p_b - p_a).GetSafeNormal();

			resultFrame.position = dir * d + frame_p;
			resultFrame.rotation = _segmentRotations[segment];
			resultFrame.segment = segment;

			return true;
		}
	}

	return false;
}

auto LinearPath::radiusSearch(FVector point, float radius) -> std::vector<std::pair<size_t, float>>
{
	std::vector<std::pair<size_t, float>> points;

	if (_brushPoints.pts.size() == 0)
		return points;

	nanoflann::SearchParams searchParams;

	// this is squared radius!
	_brushIndex->radiusSearch(&point.X, radius * radius, points, searchParams);

	return points;
}

auto LinearPath::segmentRotation(uint32_t segment) ->FQuat
{
	return _segmentRotations[segment];
}

auto LinearPath::expensiveNearestPointOnPath(const FVector point, FVector& nearest, uint32_t& segment) -> bool
{
	FVector * last = nullptr;

	if (_brushPoints.pts.size() < 2)
		return false;

	int minSegment = -1;
	float minDistance = std::numeric_limits<float>::max();

	for (int i = 0; i + 1 < _brushPoints.pts.size(); ++i)
	{
		const FVector cur = _brushPoints.pts[i].position;
		const FVector next = _brushPoints.pts[i + 1].position;

		const FVector onSegment = FMath::ClosestPointOnSegment(point, cur, next);

		const float distance = FVector::Dist(onSegment, point);

		if (distance < minDistance)
		{
			minDistance = distance;
			minSegment = i;
		}
	}

	if (minSegment >= 0)
	{
		segment = minSegment;
		nearest = _brushPoints.pts[minSegment].position;
		return true;
	}
	else
		return false;
}
