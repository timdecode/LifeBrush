#pragma once

#include <Eigen/Core>
#include <vector>

class BackgroundPoint
{
public:
	BackgroundPoint( Eigen::Vector3f _position, Eigen::Vector3f _matching_position, float _energy, float start_neighSize);
	~BackgroundPoint(void);

	Eigen::Vector3f position;
	Eigen::Vector3f matching_position;
	float energy;
	Eigen::Vector3f orientation;
	float scale;

	std::vector<Eigen::Vector3f> inputCloud_positions;
	std::vector<Eigen::Vector3f> inputCloud_attributes;

	float privateNeighSize;
	int privatePointsCounter;

};

