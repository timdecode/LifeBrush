#include "LifeBrush.h"

#include "stdafx.h"
#include "Sample.h"

using namespace Eigen;
using namespace std;

Sample::Sample(Vector3f _position, vector<Vector3f> _attributes)
{
	position = _position;
	attributes = _attributes;
	last_k = 1;// -10; // -3 replicated in SimulationManager.cpp

	feature_id = -1;

	feature_offset = Vector3f(0, 0, 0);
}


Sample::~Sample(void)
{
}
