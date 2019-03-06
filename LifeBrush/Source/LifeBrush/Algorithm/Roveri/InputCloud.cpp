#include "LifeBrush.h"

#include "stdafx.h"
#include "InputCloud.h"
#include <iostream>

using namespace Eigen;
using namespace std;

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

InputCloud::InputCloud(void)
{
	input_scale = 325;

	target_average_distance = 8.81f;
}


InputCloud::~InputCloud(void)
{
}

void InputCloud::addFeaturesSamples()
{

}



void InputCloud::loadTagAttribute()
{

}

void InputCloud::initCloud_manual()
{
	if(cloud.samples.size() == 0)
		return;

	input_scale = 1.0f;

	computeCenterOfMass();
	computeAverageSamplesDistance();
}

void InputCloud::initCloud_automatic()
{
	processModelWithScale();

	setInitialAverageDistance();

	if(abs( average_samples_distance - target_average_distance ) > 0.5) 
	{
		setInitialAverageDistance();
	}
}

void InputCloud::setInitialAverageDistance()
{
	if (average_samples_distance < target_average_distance){
		while(average_samples_distance < target_average_distance){
			input_scale+=1;
			processModelWithScale();
		}
	}else{
		while(average_samples_distance > target_average_distance){
			input_scale-=1;
			processModelWithScale();
		}
	}
}

void InputCloud::processModelWithScale()
{
	if(cloud.samples.size() == 0)
		return;

	AlignedBox3f box = unitizeSamples();

	Vector3f shift = box.center();
	float unitize_input_scale = 2.0f / box.sizes().maxCoeff();

	scaleSamples( input_scale );

	scales_for_features.push_back( unitize_input_scale * input_scale );
	shift_for_features.push_back( shift );

	computeCenterOfMass();
	computeAverageSamplesDistance();
}



Eigen::AlignedBox3f InputCloud::unitizeSamples()
{
	AlignedBox3f box = samplesBounds();

	Vector3f shift = box.center();
	float unitize_input_scale = 2.0f / box.sizes().maxCoeff();

	for(auto& sample : cloud.samples)
	{
		sample.position -= shift;
		sample.position *= unitize_input_scale;
	}

	return box;
}

void InputCloud::scaleSamples( float scaleFactor )
{
	for(auto& sample : cloud.samples)
	{
		sample.position *= scaleFactor;
	}
}


void InputCloud::assignInitialXOrientation(float angle)
{
	Matrix3f m;
	m = AngleAxisf(angle, Vector3f::UnitX())
		* AngleAxisf(0, Vector3f::UnitY())
		* AngleAxisf(0, Vector3f::UnitZ());
	for (int i=0; i<cloud.samples.size(); i++){
		Vector3f shiftedPosition = cloud.samples[i].position;
		shiftedPosition = m * shiftedPosition;
		cloud.samples[i].position = shiftedPosition;
		cloud.samples[i].attributes[0] = m * cloud.samples[i].attributes[0];
	}

	rotations_for_features.push_back(m);
}

void InputCloud::assignInitialYOrientation(float angle)
{
	Matrix3f m;
	m = AngleAxisf(0, Vector3f::UnitX())
		* AngleAxisf(angle, Vector3f::UnitY())
		* AngleAxisf(0, Vector3f::UnitZ());
	for (int i=0; i<cloud.samples.size(); i++){
		Vector3f shiftedPosition = cloud.samples[i].position;
		shiftedPosition = m * shiftedPosition;
		cloud.samples[i].position = shiftedPosition;
		cloud.samples[i].attributes[0] = m * cloud.samples[i].attributes[0];
	}

	rotations_for_features.push_back(m);
}

void InputCloud::assignInitialZOrientation(float angle)
{
	Matrix3f m;
	m = AngleAxisf(0, Vector3f::UnitX())
		* AngleAxisf(0, Vector3f::UnitY())
		* AngleAxisf(angle, Vector3f::UnitZ());
	for (int i=0; i<cloud.samples.size(); i++){
		Vector3f shiftedPosition = cloud.samples[i].position;
		shiftedPosition = m * shiftedPosition;
		cloud.samples[i].position = shiftedPosition;
		cloud.samples[i].attributes[0] = m * cloud.samples[i].attributes[0];
	}

	rotations_for_features.push_back(m);
}

Eigen::AlignedBox3f InputCloud::samplesBounds()
{
	assert( cloud.samples.size() > 0 );

	Eigen::AlignedBox3f box( cloud.samples[0].position );

	for(auto& sample : cloud.samples)
	{
		box.extend( sample.position );
	}

	return box;
}

void InputCloud::computeCenterOfMass()
{
	center_of_mass = Vector3f(0, 0, 0);
	for (int i=0;i<cloud.samples.size(); i++){
		center_of_mass = center_of_mass + cloud.samples[i].position;
	}
	center_of_mass = center_of_mass/cloud.samples.size();
}

void InputCloud::computeAverageSamplesDistance()
{
	vector<float> minDistVector;
	for (int i=0; i<cloud.samples.size(); i++){
		int closestIndex=0; 
		float minDist = std::numeric_limits<float>::max();
		for (int j=0; j<cloud.samples.size(); j++){
			float newDist = (cloud.samples[i].position-cloud.samples[j].position).norm();
			if (newDist<minDist && j!=i){
				minDist=newDist;
			}
		}
		minDistVector.push_back(minDist);
	}
	float totalDist=0;
	for (int u=0; u< minDistVector.size(); u++){
		totalDist+=minDistVector[u];
	}
	totalDist = totalDist/minDistVector.size();
	average_samples_distance = totalDist;
	//cout << "Average distance: " << average_samples_distance << endl;
}




void InputCloud::precomputeInputEnergy(float omega)
{
	precomputedInputEnergiesEigen = MatrixXf(cloud.samples.size(), cloud.samples.size());

	for(int i=0;i<precomputedInputEnergies.size();i++){
		precomputedInputEnergies[i].clear();
	}
	precomputedInputEnergies.clear();

	float dividend = 1.0f/(2.0*omega*omega);
	for (int i = 0; i<cloud.samples.size(); i++){
		vector<float> input_energies;
		for (int j = 0; j<cloud.samples.size(); j++){
			Vector3f dist_vec = cloud.samples[i].position - cloud.samples[j].position;
			float dist_value = dist_vec.dot(dist_vec);
			float pow_e = -dist_value * dividend;
			float coeff = exp1(pow_e);
			input_energies.push_back(coeff);

			precomputedInputEnergiesEigen(i, j) = coeff;
		}
		precomputedInputEnergies.push_back(input_energies);
	}
}