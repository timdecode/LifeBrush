#pragma once

#include "Cloud.h"

class OutputCloud
{
public:
	OutputCloud(void);
	~OutputCloud(void);

	void addSample( Eigen::Vector3f position, std::vector<Eigen::Vector3f> attributes);
	void addSampleWithNormal( Eigen::Vector3f position, Eigen::Vector3f normal);
	void popSample();
	void deleteSample(int index);

	Cloud cloud;

	std::vector<std::vector<int> > features_indices;
	
	std::vector<int> features_closest_bg_point;

};

