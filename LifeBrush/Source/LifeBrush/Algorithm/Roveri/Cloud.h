#pragma once

#include "Sample.h"

#include "Eigen/Dense"

#include <fstream>
#include <vector>
#include <iostream>

class Cloud
{
public:
	Cloud();
	Cloud(std::vector<Sample> _samples);
	~Cloud(void);

	void exportCloudWithNormals(char* file_name);
	void exportCloudWithNormalsTag(char* file_name);
	void exportCloudWithNormalsFeatures(char* file_name);
	void exportCloudWithNormalsTagFeatures(char* file_name);
	void exportCloudWithNormalsAndEdges(char* file_name);
	void exportCloudWithNormalsAndSplats(char* file_name);
	void exportCloudWithNormalsAndDTE(char* file_name, std::vector<Eigen::Vector3f> discreteElementsAttributes, int is_output);
	void exportCloudSGP(char* file_name, std::vector<Eigen::Vector3f> colorsAttributes);

	std::vector<Sample> samples;

	Eigen::Vector3f first_color_tag; // REPLICATED IN SIMULATIONMANAGER
	Eigen::Vector3f second_color_tag; // REPLICATED IN SIMULATIONMANAGER
};

