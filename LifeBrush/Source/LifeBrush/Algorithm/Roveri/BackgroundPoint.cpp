#include "LifeBrush.h"

#include "stdafx.h"
#include "BackgroundPoint.h"

using namespace Eigen;
using namespace std;


BackgroundPoint::BackgroundPoint(Vector3f _position, Vector3f _matching_position, float _energy, float start_neighSize)
{
	position = _position;
	matching_position = _matching_position;
	energy = _energy;
	orientation = Vector3f(0, 1, 0);
	scale = 1.0f;

	privatePointsCounter = 0;
	privateNeighSize = start_neighSize;
}


BackgroundPoint::~BackgroundPoint(void)
{
}
