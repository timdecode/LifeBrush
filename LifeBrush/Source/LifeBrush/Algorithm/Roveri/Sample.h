#pragma once

#include <Eigen/Core>
#include <vector>

class Sample
{

public:
	Sample( Eigen::Vector3f _position, std::vector<Eigen::Vector3f> _attributes);
	~Sample(void);
	
	Eigen::Vector3f position;
	std::vector<Eigen::Vector3f> attributes;

	float hackScale = 1.0f;

	int last_k;

	int feature_id;

	Eigen::Vector3f feature_offset;
};

