#include "LifeBrush.h"

#include "stdafx.h"
#include "OutputCloud.h"

using namespace Eigen;
using namespace std;

OutputCloud::OutputCloud(void)
{
}


OutputCloud::~OutputCloud(void)
{
}

void OutputCloud::addSample(Vector3f position, vector<Vector3f> attributes)
{
	Sample new_sample(position, attributes);
	cloud.samples.push_back(new_sample);
}

void OutputCloud::addSampleWithNormal(Vector3f position, Vector3f normal)
{
	vector<Vector3f> attributes;
	attributes.push_back(normal);
	addSample(position, attributes);
}

void OutputCloud::popSample()
{
	cloud.samples.pop_back();
}

void OutputCloud::deleteSample(int index)
{
	cloud.samples.erase(cloud.samples.begin() + index);
}