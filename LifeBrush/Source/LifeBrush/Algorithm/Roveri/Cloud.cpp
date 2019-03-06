#include "LifeBrush.h"

#include "stdafx.h"

#include <fstream>
#include <vector>
#include <iostream>
#include <string>

#include "Cloud.h"

using namespace Eigen;
using namespace std;

Cloud::Cloud()
{
	first_color_tag = Vector3f(1, 0, 0); // Vector3f(1, 1, 1) REPLICATED IN SIMULATIONMANAGER
	second_color_tag = Vector3f(0, 0, 1); // Vector3f(-1, -1, -1) REPLICATED IN SIMULATIONMANAGER
}

Cloud::Cloud(vector<Sample> _samples)
{
	samples = _samples;
	first_color_tag = Vector3f(1, 0, 0); // Vector3f(1, 1, 1) REPLICATED IN SIMULATIONMANAGER
	second_color_tag = Vector3f(0, 0, 1); // Vector3f(-1, -1, -1) REPLICATED IN SIMULATIONMANAGER
}


Cloud::~Cloud(void)
{
}

void Cloud::exportCloudSGP(char* file_name, vector<Vector3f> colorsAttributes){
	vector<int> colors_attributes;
	for (int i = 0; i < samples.size(); i ++){
		float minDist = 10000000000;
		int minIndex = 0;
		for (int j = 0; j < colorsAttributes.size(); j ++){
			float newDist = (colorsAttributes[j] - samples[i].attributes[1]).norm();
			if (newDist < minDist){
				minDist = newDist;
				minIndex = j;
			}
		}
		colors_attributes.push_back(minIndex);
	}

	for (int j = 0; j < colorsAttributes.size(); j ++){
		vector<int> current_samples_to_consider;
		for (int i = 0; i < samples.size(); i ++){
			if (colors_attributes[i] == j){
				current_samples_to_consider.push_back(i);
			}
		}

		if (current_samples_to_consider.size() == 0){
			continue;
		}

		std::string filename = "outputSGP_" + to_string(j) + ".ply";
		
		std::ofstream myfile (filename);
		if (myfile.is_open())
		{
			myfile << "ply\n";
			myfile << "format ascii 1.0\n";
			myfile << "comment VCGLIB generated\n";
			myfile << "element vertex " << current_samples_to_consider.size()<< "\n";
			myfile << "property float x\n";
			myfile << "property float y\n";
			myfile << "property float z\n";
			myfile << "property float nx\n";
			myfile << "property float ny\n";
			myfile << "property float nz\n";
			myfile << "element face 0\n";
			myfile << "property list uchar int vertex_indices\n";
			myfile << "end_header\n";
			for (int ii=0; ii<current_samples_to_consider.size(); ii++){
				int i = current_samples_to_consider[ii];
				myfile << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
			myfile.close();
			cout << "Generated file "<<  file_name  << endl;;
		}

	}
}


void Cloud::exportCloudWithNormalsAndDTE(char* file_name, vector<Vector3f> discreteElementsAttributes, int is_output)
{
	vector<int> samples_dte_attributes;
	for (int i = 0; i < samples.size(); i ++){
		float minDist = 10000000000;
		int minIndex = 0;
		for (int j = 0; j < discreteElementsAttributes.size(); j ++){
			float newDist = (discreteElementsAttributes[j] - samples[i].attributes[1]).norm();
			if (newDist < minDist){
				minDist = newDist;
				minIndex = j;
			}
		}
		samples_dte_attributes.push_back(minIndex);
	}

	for (int j = 0; j < discreteElementsAttributes.size(); j ++){
		vector<int> current_samples_to_consider;
		for (int i = 0; i < samples.size(); i ++){
			if (samples_dte_attributes[i] == j){
				current_samples_to_consider.push_back(i);
			}
		}

		if (current_samples_to_consider.size() == 0){
			continue;
		}

		string filename = "outputDTEs_" + to_string(j) + ".ply";
		if (is_output == 0){
			filename = "inputDTEs_" + to_string(j) + ".ply";
		}
		ofstream myfile (filename);
		if (myfile.is_open())
		{
			myfile << "ply\n";
			myfile << "format ascii 1.0\n";
			myfile << "comment VCGLIB generated\n";
			myfile << "element vertex " << current_samples_to_consider.size() * 3 << "\n";
			myfile << "property float x\n";
			myfile << "property float y\n";
			myfile << "property float z\n";
			myfile << "property float nx\n";
			myfile << "property float ny\n";
			myfile << "property float nz\n";
			myfile << "element face " << current_samples_to_consider.size() <<"\n";
			myfile << "property list uchar int vertex_indices\n";
			myfile << "end_header\n";
			for (int ii=0; ii<current_samples_to_consider.size(); ii++){
				int i = current_samples_to_consider[ii];
				

				Vector3f randv = Vector3f((rand()%10), (rand()%10), (rand()%10));
				randv.normalize();
				Vector3f tangv = samples[i].attributes[0].normalized().cross(randv);
				Vector3f tangv2 = samples[i].attributes[0].normalized().cross(tangv);

				Vector3f newPos1 = samples[i].position + tangv.normalized() * 0.1f;
				myfile << newPos1[0] << " " << newPos1[1] << " " << newPos1[2] 
					<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
				Vector3f newPos2 = samples[i].position + tangv2.normalized() * 0.1f;
				myfile << newPos2[0] << " " << newPos2[1] << " " << newPos2[2] 
					<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
				Vector3f newPos3 = samples[i].position; 
				myfile << newPos3[0] << " " << newPos3[1] << " " << newPos3[2] 
					<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
			for (int i=0; i<current_samples_to_consider.size(); i++){
				myfile << "3 " << 3 * i + 0 << " " << 3 * i + 1 << " " << 3 * i + 2   <<"\n" ;
			}
			myfile.close();
			cout << "Generated file "<<  file_name  << endl;;
		}

	}

}

void Cloud::exportCloudWithNormalsAndSplats(char* file_name)
{
	ofstream myfile ("outputSplats.ply");
	if (myfile.is_open())
	{
		myfile << "ply\n";
		myfile << "format ascii 1.0\n";
		myfile << "comment VCGLIB generated\n";
		myfile << "element vertex " << samples.size() * 3 << "\n";
		myfile << "property float x\n";
		myfile << "property float y\n";
		myfile << "property float z\n";
		myfile << "property float nx\n";
		myfile << "property float ny\n";
		myfile << "property float nz\n";
		myfile << "element face " << samples.size() <<"\n";
		myfile << "property list uchar int vertex_indices\n";
		myfile << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			
			Vector3f randv = Vector3f((rand()%10), (rand()%10), (rand()%10));
			randv.normalize();
			Vector3f tangv = samples[i].attributes[0].normalized().cross(randv);
			Vector3f tangv2 = samples[i].attributes[0].normalized().cross(tangv);

			Vector3f newPos1 = samples[i].position + tangv.normalized() * 0.1f;
			myfile << newPos1[0] << " " << newPos1[1] << " " << newPos1[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			Vector3f newPos2 = samples[i].position + tangv2.normalized() * 0.1f;
			myfile << newPos2[0] << " " << newPos2[1] << " " << newPos2[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			Vector3f newPos3 = samples[i].position; 
			myfile << newPos3[0] << " " << newPos3[1] << " " << newPos3[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;

		}
		for (int i=0; i<samples.size(); i++){
			myfile << "3 " << 3 * i + 0 << " " << 3 * i + 1 << " " << 3 * i + 2   <<"\n" ;
		}
		myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";
}


void Cloud::exportCloudWithNormalsAndEdges(char* file_name)
{
	ofstream myfile ("exportedOutputEdges.obj");
	if (myfile.is_open())
	{
		
		for (int i=0; i<samples.size(); i++){
			myfile << "v " << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] <<"\n" ;
		}
		

		for (int i=0; i<samples.size(); i++){
			float mind1=10000000;
			float mind2=10000000;
			int ind1=0;
			int ind2=0;
			for (int j=0; j<samples.size(); j++){
				float newD=sqrt(pow( samples[i].position[0]- samples[j].position[0],2)+pow( samples[i].position[1]- samples[j].position[1],2)+pow( samples[i].position[2]- samples[j].position[2],2));
				if (newD<mind1 && j!=i){
					mind1=newD;
					ind1=j;
				}
			}
			for (int j=0; j<samples.size(); j++){
				float newD=sqrt(pow( samples[i].position[0]- samples[j].position[0],2)+pow( samples[i].position[1]- samples[j].position[1],2)+pow( samples[i].position[2]- samples[j].position[2],2));
				float otherDist=sqrt(pow( samples[ind1].position[0]- samples[j].position[0],2)+pow( samples[ind1].position[1]- samples[j].position[1],2)+pow( samples[ind1].position[2]- samples[j].position[2],2));
				if (newD<mind2 && j!=ind1 && j!=i && otherDist>newD){
					mind2=newD;
					ind2=j;
				}
			}
			// cout << mind1 << " " << mind2<< endl;
			if (mind1<19 || true){
				myfile << "l " << ind1 + 1 << " " << i + 1<< "\n" ;
			}
			if (mind2<19 || true){
				myfile << "l " << ind2 + 1<< " " << i + 1<< "\n" ;
			}
		}

	myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";


}


void Cloud::exportCloudWithNormals(char* file_name)
{
	ofstream myfile (file_name);
	if (myfile.is_open())
	{
		myfile << "ply\n";
		myfile << "format ascii 1.0\n";
		myfile << "comment VCGLIB generated\n";
		myfile << "element vertex " << samples.size() << "\n";
		myfile << "property float x\n";
		myfile << "property float y\n";
		myfile << "property float z\n";
		myfile << "property float nx\n";
		myfile << "property float ny\n";
		myfile << "property float nz\n";
		myfile << "element face 0\n";
		myfile << "property list uchar int vertex_indices\n";
		myfile << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			myfile << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
			<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
		}
		myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	if (1){
		ofstream myfile ("text_for_matlab.txt");
		if (myfile.is_open())
		{
			for (int i=0; i<samples.size(); i++){
				myfile << samples[i].position[0] << " " << samples[i].position[1]  <<"\n" ;
			}
			myfile.close();
			cout << "Generated file "<<  "text_for_matlab.txt"  << endl;;
		}
		else cout << "Unable to open file";
	}
}


void Cloud::exportCloudWithNormalsTag(char* file_name)
{
	ofstream myfile ("first.ply");
	if (myfile.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_first < dist_from_second && dist_from_first < dist_from_third){	
				count ++;
			}
		}

		myfile << "ply\n";
		myfile << "format ascii 1.0\n";
		myfile << "comment VCGLIB generated\n";
		myfile << "element vertex " << count << "\n";
		myfile << "property float x\n";
		myfile << "property float y\n";
		myfile << "property float z\n";
		myfile << "property float nx\n";
		myfile << "property float ny\n";
		myfile << "property float nz\n";
		myfile << "element face 0\n";
		myfile << "property list uchar int vertex_indices\n";
		myfile << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_first < dist_from_second && dist_from_first < dist_from_third){	
				myfile << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	ofstream myfile2 ("second.ply");
	if (myfile2.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_second < dist_from_first && dist_from_second < dist_from_third){	
				count ++;
			}
		}

		myfile2 << "ply\n";
		myfile2 << "format ascii 1.0\n";
		myfile2 << "comment VCGLIB generated\n";
		myfile2 << "element vertex " << count << "\n";
		myfile2 << "property float x\n";
		myfile2 << "property float y\n";
		myfile2 << "property float z\n";
		myfile2 << "property float nx\n";
		myfile2 << "property float ny\n";
		myfile2 << "property float nz\n";
		myfile2 << "element face 0\n";
		myfile2 << "property list uchar int vertex_indices\n";
		myfile2 << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_second < dist_from_first && dist_from_second < dist_from_third){	
				myfile2 << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile2.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	ofstream myfile3 ("third.ply");
	if (myfile3.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_third < dist_from_first && dist_from_third < dist_from_second ){	
				count ++;
			}
		}

		myfile3 << "ply\n";
		myfile3 << "format ascii 1.0\n";
		myfile3 << "comment VCGLIB generated\n";
		myfile3 << "element vertex " << count << "\n";
		myfile3 << "property float x\n";
		myfile3 << "property float y\n";
		myfile3 << "property float z\n";
		myfile3 << "property float nx\n";
		myfile3 << "property float ny\n";
		myfile3 << "property float nz\n";
		myfile3 << "element face 0\n";
		myfile3 << "property list uchar int vertex_indices\n";
		myfile3 << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_third < dist_from_first && dist_from_third < dist_from_second){	
				myfile3 << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile3.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";
}

void Cloud::exportCloudWithNormalsTagFeatures(char* file_name)
{
	ofstream myfile ("first.ply");
	if (myfile.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_first < dist_from_second && dist_from_first < dist_from_third && samples[i].feature_id == -1 ){	
				count ++;
			}
		}

		myfile << "ply\n";
		myfile << "format ascii 1.0\n";
		myfile << "comment VCGLIB generated\n";
		myfile << "element vertex " << count << "\n";
		myfile << "property float x\n";
		myfile << "property float y\n";
		myfile << "property float z\n";
		myfile << "property float nx\n";
		myfile << "property float ny\n";
		myfile << "property float nz\n";
		myfile << "element face 0\n";
		myfile << "property list uchar int vertex_indices\n";
		myfile << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_first < dist_from_second && dist_from_first < dist_from_third && samples[i].feature_id == -1){	
				myfile << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	ofstream myfile2 ("second.ply");
	if (myfile2.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_second < dist_from_first && dist_from_second < dist_from_third){	
				count ++;
			}
		}

		myfile2 << "ply\n";
		myfile2 << "format ascii 1.0\n";
		myfile2 << "comment VCGLIB generated\n";
		myfile2 << "element vertex " << count << "\n";
		myfile2 << "property float x\n";
		myfile2 << "property float y\n";
		myfile2 << "property float z\n";
		myfile2 << "property float nx\n";
		myfile2 << "property float ny\n";
		myfile2 << "property float nz\n";
		myfile2 << "element face 0\n";
		myfile2 << "property list uchar int vertex_indices\n";
		myfile2 << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_second < dist_from_first && dist_from_second < dist_from_third){	
				myfile2 << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile2.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	ofstream myfile3 ("third.ply");
	if (myfile3.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_third < dist_from_first && dist_from_third < dist_from_second ){	
				count ++;
			}
		}

		myfile3 << "ply\n";
		myfile3 << "format ascii 1.0\n";
		myfile3 << "comment VCGLIB generated\n";
		myfile3 << "element vertex " << count << "\n";
		myfile3 << "property float x\n";
		myfile3 << "property float y\n";
		myfile3 << "property float z\n";
		myfile3 << "property float nx\n";
		myfile3 << "property float ny\n";
		myfile3 << "property float nz\n";
		myfile3 << "element face 0\n";
		myfile3 << "property list uchar int vertex_indices\n";
		myfile3 << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			float dist_from_first = (samples[i].attributes[1] - first_color_tag).norm();
			float dist_from_second = (samples[i].attributes[1] - second_color_tag).norm();
			float dist_from_third = (samples[i].attributes[1] - Vector3f(1, -1, 1)).norm();
			if (dist_from_third < dist_from_first && dist_from_third < dist_from_second){	
				myfile3 << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile3.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";




	ofstream myfile2_ ("featt.ply");
	if (myfile2_.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			if (samples[i].feature_id != -1 ){	
				count ++;
			}
		}

		myfile2_ << "ply\n";
		myfile2_ << "format ascii 1.0\n";
		myfile2_ << "comment VCGLIB generated\n";
		myfile2_ << "element vertex " << count << "\n";
		myfile2_ << "property float x\n";
		myfile2_ << "property float y\n";
		myfile2_ << "property float z\n";
		myfile2_ << "property float nx\n";
		myfile2_ << "property float ny\n";
		myfile2_ << "property float nz\n";
		myfile2_ << "element face 0\n";
		myfile2_ << "property list uchar int vertex_indices\n";
		myfile2_ << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			if (samples[i].feature_id != -1 ){		
				myfile2_ << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile2_.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";


}

void Cloud::exportCloudWithNormalsFeatures(char* file_name)
{
	ofstream myfile ("first.ply");
	if (myfile.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			
			if (samples[i].feature_id == -1 ){	
				count ++;
			}
		}

		myfile << "ply\n";
		myfile << "format ascii 1.0\n";
		myfile << "comment VCGLIB generated\n";
		myfile << "element vertex " << count << "\n";
		myfile << "property float x\n";
		myfile << "property float y\n";
		myfile << "property float z\n";
		myfile << "property float nx\n";
		myfile << "property float ny\n";
		myfile << "property float nz\n";
		myfile << "element face 0\n";
		myfile << "property list uchar int vertex_indices\n";
		myfile << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			if (samples[i].feature_id == -1 ){	
				myfile << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

	ofstream myfile2 ("second.ply");
	if (myfile2.is_open())
	{
		int count = 0;
		for (int i=0; i<samples.size(); i++){
			if (samples[i].feature_id != -1 ){	
				count ++;
			}
		}

		myfile2 << "ply\n";
		myfile2 << "format ascii 1.0\n";
		myfile2 << "comment VCGLIB generated\n";
		myfile2 << "element vertex " << count << "\n";
		myfile2 << "property float x\n";
		myfile2 << "property float y\n";
		myfile2 << "property float z\n";
		myfile2 << "property float nx\n";
		myfile2 << "property float ny\n";
		myfile2 << "property float nz\n";
		myfile2 << "element face 0\n";
		myfile2 << "property list uchar int vertex_indices\n";
		myfile2 << "end_header\n";
		for (int i=0; i<samples.size(); i++){
			if (samples[i].feature_id != -1 ){		
				myfile2 << samples[i].position[0] << " " << samples[i].position[1] << " " << samples[i].position[2] 
				<< " " << samples[i].attributes[0][0] << " " << samples[i].attributes[0][1]  << " " << samples[i].attributes[0][2]  <<"\n" ;
			}
		}
		myfile2.close();
		cout << "Generated file "<<  file_name  << endl;;
	}
	else cout << "Unable to open file";

}

