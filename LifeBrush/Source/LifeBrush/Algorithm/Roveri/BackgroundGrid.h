#pragma once

#include "BackgroundPoint.h"
#include <vector>

class BackgroundGrid
{
public:
	BackgroundGrid(void);
	BackgroundGrid(std::vector<BackgroundPoint> _points);
	~BackgroundGrid(void);

	std::vector<BackgroundPoint> points;
	std::vector<BackgroundPoint> points_large;
	std::vector<int> points_large_taken;

	void generateBackgroundPoints(float gridSize, float cellSize, bool is_2d, float start_neighSize);

	float overlap_factor;

	std::vector<int> brush_indices;
};

