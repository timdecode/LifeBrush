#include "LifeBrush.h"

#include "stdafx.h"
#include "BackgroundGrid.h"

using namespace Eigen;
using namespace std;

BackgroundGrid::BackgroundGrid()
{
	overlap_factor = 1.65;// 1 means no overlap 
}

BackgroundGrid::BackgroundGrid(vector<BackgroundPoint> _points)
{
	points = _points;
}


BackgroundGrid::~BackgroundGrid(void)
{
}

void BackgroundGrid::generateBackgroundPoints(float halfGridSize, float cellSize, bool is_2d, float start_neighSize)
{
	float points_density = floor(2 * overlap_factor *halfGridSize/(cellSize));
	for (int i=0; i<=points_density; i++){
		for (int j=0; j<=points_density;j++){
			int depth=points_density;
			if (is_2d){
				depth=0;
			}
			for (int h=0; h<=depth;h++){
				float xx=-halfGridSize+2*halfGridSize/points_density*i;
				float yy=-halfGridSize+2*halfGridSize/points_density*j;
				float zz=-halfGridSize+2*halfGridSize/points_density*h;
				if (is_2d){
					zz=0;
				}
				Vector3f position(xx, yy, zz);
				Vector3f matching_position(0, 0, 0);
				float energy = 0;
				BackgroundPoint new_point(position, matching_position, energy, start_neighSize);
				points.push_back(new_point);
			}

		}

	}
	points_large = points;
	points_large_taken = vector<int>(points_large.size(), 0);
	points.clear();
	brush_indices.clear();
}