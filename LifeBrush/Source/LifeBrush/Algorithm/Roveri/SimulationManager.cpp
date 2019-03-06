#include "LifeBrush.h"

#include "stdafx.h"
#include "SimulationManager.h"
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <ctime>

using namespace Eigen;
using namespace std;

SimulationManager::SimulationManager()
{
	

	//------------------------------------------------INSTRUCTIONS TO RUN:---------------------------------------------------

	 // 1. Press 'p' on the input window
	 // 2. draw a single stroke on the output window, press 'F' to generate the spline, press 'f' to load the background grid
	 // 3. manually add first points with '0' on the stroke or with '9' along the stroke.
	 // 4. if a mesh is present in the folder BackgroundMesh, do not draw the stroke and just press 'f' in step 2.
	 // 5. press 'c' to export

	//-----------------------------------------------------------------------------------------------------------------------




	 is_2d = false;

	 allow_tag_optimization = true;

	 allow_rotation = true;
	 allow_scale = false;

	 disk_render_size = 0.65f;
	 output_canvas_size = 750;//750
	 neighSize = 110; //35
	 brushSize = 100;

	 exp_approx_thres = 10;

	 assignStep_precision = 0.005f;
	 assignStep_max_iterations = 10;
	 assignStep_eps = 0.031f;
	 matchingStep_precision = 0.005f;
	 matchingStep_max_iterations = 10;
	 matchingStep_eps = 0.031f;
	 assignStepAttribute_precision = 0.005f;
	 assignStepAttribute_max_iterations = 5;
	 assignStepAttribute_eps = 0.0015;

	 draw_background_grid = 1;
	 draw_background_energies = 0;
	 draw_output_points = 1;
	 draw_matching_points = 1;

//	 use_periodic_borders = 1;
	 periodic_border_type = PeriodicBorderType::Radius;


	 allow_add_points = 1;
	 allow_delete_points = 1;

	 use_input_precomputation = 0;
	 use_output_precomputation = 0;

	 if (allow_scale){
		 use_input_precomputation = 0;
		 use_output_precomputation = 0;
	 }

	 show_average_energy = 1;

	 matching_position_new_seeds = 3;

	 accurate_add_delete = 0;

	 use_only_active_output_points = 0;

	 automatic_multiscale_opt = 0;

	 // SGP Scene:
	 // 0 - 3d points with normals, no color tag
	 // 1 - 3d points with normals, color tag



	 automatic_add_features = 0;

	 consider_attributes_in_position_optimization = 0;

	 show_mouse = 1;

	 sim_step = 0;

	 autoExport = false;

	 is_paused = false;

	 factorTS = 1;



	 scene = 1; // nothing -1, liana 0, details on branch 1, 3d chair 2, features wall 3, ivy 10 , carnival 11, smoke 12, DTE 13, text synt 14, pallini 15, waves 16, bone 17
}


SimulationManager::~SimulationManager(void)
{
}

void SimulationManager::clear()
{
	{
		float target_average_distance = inputCloud.target_average_distance;
		inputCloud = InputCloud();
		inputCloud.target_average_distance = target_average_distance;
	}

	outputCloud = OutputCloud();

	{
		float overlap_factor = backgroundGrid.overlap_factor;
		backgroundGrid = BackgroundGrid();
		backgroundGrid.overlap_factor = overlap_factor;
	}

	brush_line_points.clear();
	averageEnergies.clear();
	discreteElementsAttributes.clear();
}

void SimulationManager::init()
{
	SGPScene = 0;
	if(inputCloud.cloud.samples[0].attributes.size() > 1) {
		SGPScene = 1;
	}

	//-------------------------------------------------NOT USED STUFF:---------------------------------------------------------------

	first_color_tag = Vector3f( 1, 0, 0 ); // Vector3f(1, 1, 1) REPLICATED IN CLOUD
	second_color_tag = Vector3f( 0, 0, 1 ); // Vector3f(-1, -1, -1) REPLICATED IN CLOUD

	if(scene == 0) {
		for(int i = 0; i < inputCloud.cloud.samples.size(); i++) {
			if(i < 692) {
				inputCloud.cloud.samples[i].attributes.push_back( first_color_tag );
			}
			else {
				inputCloud.cloud.samples[i].attributes.push_back( second_color_tag );
			}
		}
		for(int i = 0; i < 4; i++) {
			int rand_index = rand() % (inputCloud.cloud.samples.size() - 1);
			inputCloud.cloud.samples[rand_index].attributes[1] = Vector3f( 1, -1, 1 );
		}
	}

	if(scene == 10 && false) {
		for(int i = 0; i < inputCloud.cloud.samples.size(); i++) {
			if(i < 257) {
				inputCloud.cloud.samples[i].attributes.push_back( first_color_tag );
			}
			else if(i < (899999)) {
				inputCloud.cloud.samples[i].attributes.push_back( second_color_tag );
			}
			else {
				float minDist = 1000000000;
				int minIndex = 0;
				for(int k = 0; k < inputCloud.cloud.samples.size(); k++) {
					float newDist = (inputCloud.cloud.samples[k].position - inputCloud.cloud.samples[i].position).norm();
					if(newDist < minDist) {
						minDist = newDist;
						minIndex = k;
					}
				}
				inputCloud.cloud.samples[i].attributes.push_back( inputCloud.cloud.samples[minIndex].attributes[1] );
			}
		}
	}
	if(scene == 13) {

		int numberOfDiscreteElements = 5;

		for(int u = 0; u < numberOfDiscreteElements; u++) {
			discreteElementsAttributes.push_back( Vector3f( rand() % 10 - 5, rand() % 10 - 5, rand() % 10 - 5 ).normalized() );
		}
		for(int i = 0; i < inputCloud.cloud.samples.size(); i++) {
			int randi = rand() % numberOfDiscreteElements;
			inputCloud.cloud.samples[i].attributes.push_back( discreteElementsAttributes[randi] );
		}
	}
}

void SimulationManager::autoExportCloud()
{

}




bool fisttimeever = true;

struct sort_pred {
    bool operator()(const std::pair<int,int> &left, const std::pair<int,int> &right) {
        return left.first < right.first;
    }
};



void SimulationManager::preprocessInput()
{

	dteBigArray.clear();
	for (int o = 0; o < outputCloud.cloud.samples.size(); o++){
		dteBigArray.push_back(vector<int>(outputCloud.cloud.samples.size(), -1));
	}

	vector<vector<int> > neighborsInputDTE;
	for (int p = 0; p < inputCloud.cloud.samples.size(); p++){

		vector<int> current_neighbors;

		for (int pp = 0; pp < inputCloud.cloud.samples.size(); pp++){
			float dist = (inputCloud.cloud.samples[p].position - inputCloud.cloud.samples[pp].position).norm();
			if (dist < neighSize && pp != p){
				current_neighbors.push_back(pp);
			}
		}

		neighborsInputDTE.push_back(current_neighbors);

	}
	
	kco_neighbors = vector<vector<int> >(inputCloud.cloud.samples.size());

	int kco = 7;

	#pragma omp parallel for
	for (int o = 0; o < inputCloud.cloud.samples.size(); o++){

		vector<pair<int,int> > kco_energies;

		for (int i = 0; i < inputCloud.cloud.samples.size(); i++){

			float _energ = 100000000000000;
			vector<int> _pairs;
			int _matchIndex = 0;

			vector<int> usedIndices;
			float totEnergy = 0;
			vector<int> currentPairs(inputCloud.cloud.samples.size(), 0);

			for (int ooo = 0; ooo < neighborsInputDTE[o].size(); ooo++){

				int oo = neighborsInputDTE[o][ooo];
				float minDist = 10000000000000;
				int minIndex = 0;
				for (int ii = 0; ii < inputCloud.cloud.samples.size(); ii++){

					if (ii == i){
						continue;
					}

					float newDist = ((inputCloud.cloud.samples[ii].position - inputCloud.cloud.samples[i].position) - (inputCloud.cloud.samples[oo].position - inputCloud.cloud.samples[o].position)).norm();
					if (newDist < minDist){
						bool busy = false;
						for (int u = 0; u < usedIndices.size(); u++){
							if (usedIndices[u] == ii){
								busy = true;
								break;
							}
						}
						if (!busy){
							minDist = newDist;
							minIndex = ii; 
						}
					}
				}
				usedIndices.push_back(minIndex);
				totEnergy += minDist;
				currentPairs[oo] = minIndex;
			}

			kco_energies.push_back(make_pair(totEnergy, i));

		}

		std::sort(kco_energies.begin(), kco_energies.end(), sort_pred());

		vector<int> kco_neighbors_temp;
		for (int y = 0; y < kco; y ++){
			if (y < kco_energies.size()){
				kco_neighbors_temp.push_back(kco_energies[y].second);
				//cout <<kco_energies[y].second << " ";
			}
		}
		//cout << endl;
		kco_neighbors[o] = kco_neighbors_temp;
		
	}
	cout <<"input samples "<< inputCloud.cloud.samples.size() << endl;
} 

void SimulationManager::multiScaleAdapt()
{
	if (automatic_multiscale_opt == 0){
		return;
	}

	if (sim_step % 10 == 0){
		for (int ii=0; ii<backgroundGrid.brush_indices.size(); ii++){
			int i_ = backgroundGrid.brush_indices[ii];
			if (closeOutputsIndicesGlobal[ii].size() > 5){
				if (abs((int)closeOutputsIndicesGlobal[ii].size() - backgroundGrid.points[i_].privatePointsCounter) < 4 ){
					if (0.9 * backgroundGrid.points[i_].privateNeighSize > 0.9*((backgroundGrid.points[0].position-backgroundGrid.points[1].position).norm())){
						backgroundGrid.points[i_].privateNeighSize = 0.9 * backgroundGrid.points[i_].privateNeighSize;
						//omega = 0.95 * omega;
					}
				}
				backgroundGrid.points[i_].privatePointsCounter = closeOutputsIndicesGlobal[ii].size();
			}
		}
	}
}




void SimulationManager::run()
{
	// Pause simulation
	if (is_paused){
		return;
	}

	// Increase step counter
	sim_step ++;

	// Auto export
	if (autoExport){
		if (sim_step % 20 == 0){
			autoExportCloud();
		}
	}

	// Adapt parameters for multiscale optimization
	multiScaleAdapt();

	// Find neighbors
	computeNeighbors();

	// Set indices of active samples
	vector<int> output_points_indices;
	if (use_only_active_output_points == 0){
		output_points_indices = generateOutputIndicesBasedOnAllIndices();
	} else{
		//output_points_indices = generateOutputIndicesBasedOnNonZeroK();
		output_points_indices = generateOutputIndicesBasedOnBrushPosition();
	}
	vector<int> matching_points_indices;
	if (use_only_active_output_points == 0){
		matching_points_indices = generateMatchingIndicesBasedOnAllIndices();
	} else{
		//matching_points_indices = generateMatchingIndicesBasedOnNonZeroK();
		matching_points_indices = generateMatchingIndicesBasedOnBrushPosition();
	}
	
	// 2 Steps iteration
	assignStep(output_points_indices);
	//assignFeaturesStep();
	matchingStep(matching_points_indices);
	
	// Add and remove samples
	if (allow_add_points || allow_delete_points){
		computeEnergy();
		if (allow_add_points){
			addStep();
		}
		if (automatic_add_features == 1){
			addFeaturesStep();
		}
		if (allow_delete_points){
			if (outputCloud.cloud.samples.size() > 10){
				deleteBadSamples();
			}
		}
		if (automatic_delete_features == 1){
			deleteFeaturesStep();
		}
	}

}

vector<int> SimulationManager::findInputNeighbors(int ii)
{
	//Find indices of input samples close to the matching point position of the background point ii

	int i_ = backgroundGrid.brush_indices[ii];
	float scaled_neigh_size = neighSize * backgroundGrid.points[i_].scale;
	vector<int> closeInputsIndices;
	for (int pp = 0; pp < inputCloud.cloud.samples.size(); pp++){
		float dist = (backgroundGrid.points[i_].matching_position - backgroundGrid.points[i_].inputCloud_positions[pp]).norm();
		//if (dist < scaled_neigh_size){
		if (dist < backgroundGrid.points[i_].privateNeighSize){
			closeInputsIndices.push_back(pp);
		}
	}
	return closeInputsIndices;
}

void SimulationManager::computeNeighbors()
{
	//Find indices of output samples close to the position of each background point

	closeOutputsIndicesGlobal.clear();
	closeInputsIndicesGlobal.clear();
	closeBrushIndicesGlobal.clear();

	closeBrushIndicesGlobal.resize(outputCloud.cloud.samples.size());

	for (int ii=0; ii<backgroundGrid.brush_indices.size(); ii++){
		int i_ = backgroundGrid.brush_indices[ii];

		float scaled_neigh_size = neighSize * backgroundGrid.points[i_].scale;

		vector<int> closeOutputsIndices;
		for (int pp = 0; pp < outputCloud.cloud.samples.size(); pp++){
			float dist = (backgroundGrid.points[i_].position - outputCloud.cloud.samples[pp].position).norm();
			//if (dist < scaled_neigh_size){
			if (dist < backgroundGrid.points[i_].privateNeighSize){
				closeOutputsIndices.push_back(pp);
				closeBrushIndicesGlobal[pp].push_back(ii);
			}
		}
		closeOutputsIndicesGlobal.push_back(closeOutputsIndices);

		vector<int> closeInputsIndices = findInputNeighbors(ii);
		closeInputsIndicesGlobal.push_back(closeInputsIndices);
	}

}

void SimulationManager::computeEnergy()
{
	// Compute the total energy

	vector<int> matching_points_indices;
	if (use_only_active_output_points == 0){
		matching_points_indices = generateMatchingIndicesBasedOnAllIndices();
	} else{
		//matching_points_indices = generateMatchingIndicesBasedOnNonZeroK();
		matching_points_indices = generateMatchingIndicesBasedOnBrushPosition();
	}

	#pragma omp parallel for
	for (int ii = 0; ii < matching_points_indices.size(); ii++){
		int i = matching_points_indices[ii];
		setEnergyToPoint(i);
	}

	if (show_average_energy == 1){
		computeAverageEnergy();
	}
}

void SimulationManager::computeAverageEnergy()
{
	// Compute the average energy

	float average_energy = 0;
	int active_bg_points_counter = 0;
	for (int ii = 0; ii < backgroundGrid.brush_indices.size(); ii++){
		int i = backgroundGrid.brush_indices[ii];
		if (backgroundGrid.points[i].energy != 0){
			average_energy += backgroundGrid.points[i].energy;
			active_bg_points_counter ++;
		}
	}
	if (active_bg_points_counter != 0){
		average_energy = average_energy / (float)active_bg_points_counter;
		averageEnergies.push_back(average_energy);
	}
}

void SimulationManager::setEnergyToPoint(int index)
{
	// Assign energy to background point index

	float energy_val = computeEnergyForPoint(index);
	int i = backgroundGrid.brush_indices[index];
	backgroundGrid.points[i].energy = energy_val;
}

float SimulationManager::computeEnergyForPoint(int index_brush_indices)
{
	// Compute the enrgy for the background point index_brush_indices

	if (closeOutputsIndicesGlobal[index_brush_indices].size() == 0){
		return 0;
	}

	int i = backgroundGrid.brush_indices[index_brush_indices];

	//float prefix = (float)sqrt(PI*omega*omega/2.0);
	float scaled_omega = omega * backgroundGrid.points[i].scale;
	float scaled_normalized_prefix = (float)sqrt(PI/2.0) * omega;// * scaled_omega;
	float scaled_dividend = 1.0f/(2.0*scaled_omega*scaled_omega);
	float scaled_local_exp_approx_thres = exp_approx_thres / scaled_dividend;
	float en = 0;

	//omp_set_nested(1);

	#pragma omp parallel for reduction(+:en)
	for (int jj=0; jj<closeInputsIndicesGlobal[index_brush_indices].size(); jj++){
		int j=closeInputsIndicesGlobal[index_brush_indices][jj];
		//#pragma omp parallel for
		for (int kk=0; kk<closeInputsIndicesGlobal[index_brush_indices].size(); kk++){
			int k=closeInputsIndicesGlobal[index_brush_indices][kk];

			if (use_input_precomputation){
				float new_en = inputCloud.precomputedInputEnergies[j][k];
				en = en + new_en;
			}
			else{
				Vector3f dist_vec =  backgroundGrid.points[i].inputCloud_positions[j] -  backgroundGrid.points[i].inputCloud_positions[k];
				float dist_value = dist_vec.dot(dist_vec);
				if (dist_value > scaled_local_exp_approx_thres){
					continue;
				}
				float coeff = exp1(-dist_value * scaled_dividend);
				if (consider_attributes_in_position_optimization){
					//float attributeCoeff = backgroundGrid.points[i].inputCloud_attributes[j].dot(backgroundGrid.points[i].inputCloud_attributes[k]); coeff *= attributeCoeff;
					float attributeCoeff = inputCloud.cloud.samples[j].attributes[1].dot(inputCloud.cloud.samples[k].attributes[1]); coeff *= attributeCoeff;
				}
				en = en + coeff;
			}
		}
	}
	#pragma omp parallel for reduction(+:en)
	for (int jj=0; jj<closeInputsIndicesGlobal[index_brush_indices].size(); jj++){
		int j=closeInputsIndicesGlobal[index_brush_indices][jj];

		Vector3f dist_vec_first_part = backgroundGrid.points[i].inputCloud_positions[j] - backgroundGrid.points[i].matching_position + backgroundGrid.points[i].position;

		//#pragma omp parallel for
		for (int kk=0; kk<closeOutputsIndicesGlobal[index_brush_indices].size(); kk++){
			int k=closeOutputsIndicesGlobal[index_brush_indices][kk];

			Vector3f dist_vec = dist_vec_first_part - outputCloud.cloud.samples[k].position;
			float dist_value = dist_vec.dot(dist_vec);

			if (dist_value > scaled_local_exp_approx_thres){
				continue;
			}
			float coeff = exp1(-dist_value * scaled_dividend);
			if (consider_attributes_in_position_optimization){
				//float attributeCoeff = backgroundGrid.points[i].inputCloud_attributes[j].dot(outputCloud.cloud.samples[k].attributes[0]); coeff *= attributeCoeff;
				float attributeCoeff = inputCloud.cloud.samples[j].attributes[1].dot(outputCloud.cloud.samples[k].attributes[1]); coeff *= attributeCoeff;
			}
			en = en - 2 * coeff;
		}
	}
	#pragma omp parallel for reduction(+:en)
	for (int kk=0; kk<closeOutputsIndicesGlobal[index_brush_indices].size(); kk++){
		int k=closeOutputsIndicesGlobal[index_brush_indices][kk];
	//	#pragma omp parallel for
		for (int yy=0; yy<closeOutputsIndicesGlobal[index_brush_indices].size(); yy++){
			int y=closeOutputsIndicesGlobal[index_brush_indices][yy];
			Vector3f dist_vec = outputCloud.cloud.samples[k].position - outputCloud.cloud.samples[y].position;

			float dist_value = dist_vec.dot(dist_vec);

			if (dist_value > scaled_local_exp_approx_thres){
				continue;
			}
			float coeff = exp1(-dist_value * scaled_dividend);
			if (consider_attributes_in_position_optimization){
				//float attributeCoeff = outputCloud.cloud.samples[y].attributes[0].dot(outputCloud.cloud.samples[k].attributes[0]); coeff *= attributeCoeff;
				float attributeCoeff = outputCloud.cloud.samples[y].attributes[1].dot(outputCloud.cloud.samples[k].attributes[1]); coeff *= attributeCoeff;
			}
			en = en + coeff;
		}
	}
		
	en = scaled_normalized_prefix * en;
	return en;

}

float SimulationManager::computeEnergyGivenFromOnePoint(int bg_point_index, int output_point_index)
{
	// Compute the energy contribution of the sample output_point_index for the background point bg_point_index

	int index = backgroundGrid.brush_indices[bg_point_index];

	float scaled_omega = omega * backgroundGrid.points[index].scale;
	float scaled_normalized_prefix = (float)sqrt(PI/2.0) * omega;// * scaled_omega;
	float scaled_dividend = 1.0f/(2.0*scaled_omega*scaled_omega);
	float scaled_local_exp_approx_thres = exp_approx_thres / scaled_dividend;
	float val = 0;

	for (int ii_=0; ii_<closeInputsIndicesGlobal[bg_point_index].size(); ii_++){
		int i_=closeInputsIndicesGlobal[bg_point_index][ii_];
		for (int jj_=0; jj_<1; jj_++){
			int j_=output_point_index;
			
			Vector3f dist_vec =  backgroundGrid.points[index].inputCloud_positions[i_] - backgroundGrid.points[index].matching_position - outputCloud.cloud.samples[j_].position + backgroundGrid.points[index].position;
			float dist_value = dist_vec.norm();
			dist_value *= dist_value;

			if (dist_value > scaled_local_exp_approx_thres){
				continue;
			}
			float coeff = exp1(-dist_value * scaled_dividend);
			if (consider_attributes_in_position_optimization){
				//float attributeCoeff = backgroundGrid.points[index].inputCloud_attributes[i_].dot(outputCloud.cloud.samples[j_].attributes[0]); coeff *= attributeCoeff;
				float attributeCoeff = inputCloud.cloud.samples[i_].attributes[1].dot(outputCloud.cloud.samples[j_].attributes[1]); coeff *= attributeCoeff;
			}
			val = val - 2 * coeff;
		}
	}

	for (int ii_=0; ii_<closeOutputsIndicesGlobal[bg_point_index].size(); ii_++){
		int i_=closeOutputsIndicesGlobal[bg_point_index][ii_];
		for (int jj_=0; jj_<1; jj_++){
			int j_=output_point_index;
			float optimization_coeff = 2;
			if (j_ == i_){
				optimization_coeff=1;
			}
			Vector3f dist_vec = outputCloud.cloud.samples[i_].position - outputCloud.cloud.samples[j_].position;
			float dist_value = dist_vec.norm();
			dist_value *= dist_value;
			if (dist_value > scaled_local_exp_approx_thres){
				continue;
			}
			float coeff = exp1(-dist_value * scaled_dividend);
			if (consider_attributes_in_position_optimization){
				//float attributeCoeff = outputCloud.cloud.samples[j_].attributes[0].dot(outputCloud.cloud.samples[i_].attributes[0]); coeff *= attributeCoeff;
				float attributeCoeff = outputCloud.cloud.samples[j_].attributes[1].dot(outputCloud.cloud.samples[i_].attributes[1]); coeff *= attributeCoeff;
			}
			val = val + coeff * optimization_coeff;
		}
	}

	val = scaled_normalized_prefix*val;
	return val;

}

Vector3f SimulationManager::pickNewRandomMatchingPosition(int bg_point_index)
{
	// Choose a new random candidate position

	int rand_index = rand() % (backgroundGrid.points[bg_point_index].inputCloud_positions.size() - 1);
	return backgroundGrid.points[bg_point_index].inputCloud_positions[rand_index];
}

void SimulationManager::matchingStep(vector<int> indices_to_consider)
{
	// Find the matching position for each background point with gradient descent

	float prefix = (float)sqrt(PI/2.0);
	//float dividend = 1.0f/(2.0*omega*omega);

	#pragma omp parallel for
	for (int iii = 0; iii < indices_to_consider.size(); iii++){

		int ii = indices_to_consider[iii];
		int i = backgroundGrid.brush_indices[ii];

		if (closeInputsIndicesGlobal[ii].size() == 0){
			continue;
		}

		float scaled_omega = omega * backgroundGrid.points[i].scale;
		float scaled_dividend = 1.0f / (2.0 * scaled_omega * scaled_omega);
		float scaled_prefix = prefix * (-2.0 / scaled_omega);
		float scaled_max_distance_allowed = inputCloud.max_distance_allowed * backgroundGrid.points[i].scale;

		float local_exp_approx_thres = exp_approx_thres/scaled_dividend;

		float minEnergy = 1000000000;
		Vector3f minPosition = backgroundGrid.points[i].matching_position;
		vector<int> minNeighbors = closeInputsIndicesGlobal[ii];
		for (int tries = 0; tries < matching_position_new_seeds + 1; tries ++){
			Vector3f old_pos(1000, 1000, 1000);
			Vector3f new_pos;
			if (tries == 0){
				new_pos = backgroundGrid.points[i].matching_position;
			}else{
				new_pos = pickNewRandomMatchingPosition(i); new_pos = new_pos + Vector3f(neighSize/3.0 * ((double) rand() / (RAND_MAX) * 2 - 1), neighSize/3.0 * ((double) rand() / (RAND_MAX)* 2 - 1), neighSize/3.0 * ((double) rand() / (RAND_MAX)* 2 - 1));
				vector<int> closeInputsIndices = findInputNeighbors(ii);
				closeInputsIndicesGlobal[ii].clear();
				closeInputsIndicesGlobal[ii] = closeInputsIndices;
			}

			//Vector3f old_pos(1000, 1000, 1000);
			//Vector3f new_pos = backgroundGrid.points[i].matching_position;
			int k = 0;
			while ((new_pos - old_pos).norm() > matchingStep_precision && k < matchingStep_max_iterations){
				old_pos = new_pos;
				Vector3f new_val(0, 0, 0);

				for (int jj=0; jj<closeInputsIndicesGlobal[ii].size(); jj++){
					int j=closeInputsIndicesGlobal[ii][jj];

					Vector3f dist_vec_first_part = backgroundGrid.points[i].inputCloud_positions[j] - old_pos + backgroundGrid.points[i].position;

					for (int kk=0; kk<closeOutputsIndicesGlobal[ii].size(); kk++){
						int k_=closeOutputsIndicesGlobal[ii][kk];

						Vector3f dist_vec = dist_vec_first_part - outputCloud.cloud.samples[k_].position;
						float dist_value = dist_vec.dot(dist_vec);

						if (dist_value > local_exp_approx_thres){
							continue;
						}
						float coeff = exp1(-dist_value*scaled_dividend);
						if (consider_attributes_in_position_optimization){
							//float attributeCoeff = backgroundGrid.points[i].inputCloud_attributes[j].dot(outputCloud.cloud.samples[k_].attributes[0]); coeff *= attributeCoeff; 
							float attributeCoeff = inputCloud.cloud.samples[j].attributes[1].dot(outputCloud.cloud.samples[k_].attributes[1]); coeff *= attributeCoeff; 
						}
						new_val = new_val + dist_vec * coeff;
					}
				}
				new_pos = old_pos - matchingStep_eps * new_val * scaled_prefix;// / backgroundGrid.points[i].scale; DOF!!
				if (is_2d){
					new_pos[2] = 0;
				}
				k++;
			}

			if (periodic_border_type == PeriodicBorderType::Radius)
			{
				if (new_pos.norm() > scaled_max_distance_allowed) {
					new_pos = new_pos.normalized() * scaled_max_distance_allowed;
				}
			}
			else if (periodic_border_type == PeriodicBorderType::Bounds)
			{
				float exteriorDistance = periodic_border_bounds.exteriorDistance(new_pos);

				float newLength = new_pos.norm() - exteriorDistance;

				new_pos = new_pos.normalized() * newLength;
			}

			backgroundGrid.points[i].matching_position = new_pos;

			float finalEn = -1;
			if (matching_position_new_seeds != 0){
				finalEn = computeEnergyForPoint(ii);
			}
			if (finalEn < minEnergy && finalEn != 0){
				minEnergy = finalEn;
				minPosition = new_pos;
				minNeighbors = closeInputsIndicesGlobal[ii];
			}
		}
		backgroundGrid.points[i].matching_position = minPosition;
		closeInputsIndicesGlobal[ii].clear();
		closeInputsIndicesGlobal[ii] = minNeighbors;

	}
}

void SimulationManager::assignFeaturesStep()
{
	// Assign step for the features

	float prefix = (float)sqrt(PI/2.0);
	#pragma omp parallel for
	for (int i_f = 0; i_f < outputCloud.features_indices.size(); i_f ++){
		Vector3f old_pos(1000, 1000, 1000);
		Vector3f new_pos(1, 1, 1);
		int k = 0;

		while ((new_pos - old_pos).norm() > assignStep_precision && k < assignStep_max_iterations){
			old_pos = outputCloud.cloud.samples[outputCloud.features_indices[i_f][0]].position;
			Vector3f new_val(0, 0, 0);

			for (int ii = 0; ii < outputCloud.features_indices[i_f].size(); ii ++){
				int i = outputCloud.features_indices[i_f][ii];

				for (int jjj=0; jjj<closeBrushIndicesGlobal[i].size(); jjj++){
				int jj=closeBrushIndicesGlobal[i][jjj];
				int j=backgroundGrid.brush_indices[jj];

				float scaled_omega = omega * backgroundGrid.points[j].scale;
				float scaled_dividend = 1.0f / ( 2.0 * scaled_omega * scaled_omega);
				float scaled_prefix = prefix * ( -2.0 / scaled_omega);// / backgroundGrid.points[j].scale; DOF!!

				float local_exp_approx_thres = exp_approx_thres/scaled_dividend;

				Vector3f dist_vec_first_part = - backgroundGrid.points[j].matching_position - old_pos - outputCloud.cloud.samples[i].feature_offset +  backgroundGrid.points[j].position;

				Vector3f scaled_new_val(0, 0, 0);

				/*
				MatrixXf dist_vec_first_part_matrix(3, closeInputsIndicesGlobal[jj].size());
				dist_vec_first_part_matrix = dist_vec_second_part_global[jj];
				dist_vec_first_part_matrix.colwise() += dist_vec_first_part;
				dist_vec_first_part_matrix = dist_vec_first_part_matrix * (VectorXf)(((((dist_vec_first_part_matrix.array().pow(2)).colwise().sum())* (-scaled_dividend)).array().exp()).transpose());
				scaled_new_val = scaled_new_val + dist_vec_first_part_matrix.rowwise().sum();
				*/
				
				//#pragma omp parallel for
				for (int dd=0; dd<closeInputsIndicesGlobal[jj].size(); dd++){
					int d=closeInputsIndicesGlobal[jj][dd];

					Vector3f dist_vec_second_part = backgroundGrid.points[j].inputCloud_positions[d];
					Vector3f dist_vec = dist_vec_second_part + dist_vec_first_part;
					float dist_value = dist_vec.dot(dist_vec);

					if (dist_value > local_exp_approx_thres){
						continue;
					}
					float coeff = exp1(-dist_value * scaled_dividend);
					if (consider_attributes_in_position_optimization){
						//float attributeCoeff = backgroundGrid.points[j].inputCloud_attributes[d].dot(outputCloud.cloud.samples[i].attributes[0]); coeff *= attributeCoeff;
						float attributeCoeff = inputCloud.cloud.samples[d].attributes[1].dot(outputCloud.cloud.samples[i].attributes[1]); coeff *= attributeCoeff;
					}
					scaled_new_val = scaled_new_val + dist_vec * coeff;
				}
				
				if (!use_output_precomputation || true){
					//#pragma omp parallel for
					for (int ee=0; ee<closeOutputsIndicesGlobal[jj].size(); ee++){
						int e = closeOutputsIndicesGlobal[jj][ee];
						Vector3f dist_vec = old_pos + outputCloud.cloud.samples[i].feature_offset - outputCloud.cloud.samples[e].position;
						float dist_value = dist_vec.dot(dist_vec);

						if (dist_value > local_exp_approx_thres){
							continue;
						}
						float coeff = exp1(-dist_value * scaled_dividend);
						if (consider_attributes_in_position_optimization){
							//float attributeCoeff = outputCloud.cloud.samples[e].attributes[0].dot(outputCloud.cloud.samples[i].attributes[0]); coeff *= attributeCoeff;
							float attributeCoeff = outputCloud.cloud.samples[e].attributes[1].dot(outputCloud.cloud.samples[i].attributes[1]); coeff *= attributeCoeff;
						}
						scaled_new_val = scaled_new_val + dist_vec * coeff;
					}
				}

				new_val = new_val + scaled_new_val * scaled_prefix;
			}


			new_pos = old_pos - assignStep_eps * new_val;
			if (is_2d){
				new_pos[2] = 0;
			}
			
			outputCloud.cloud.samples[outputCloud.features_indices[i_f][0]].position = new_pos;

			for (int yu = 0; yu < outputCloud.features_indices[i_f].size(); yu ++){
				int iu = outputCloud.features_indices[i_f][yu];
				outputCloud.cloud.samples[iu].position = outputCloud.cloud.samples[outputCloud.features_indices[i_f][0]].position + outputCloud.cloud.samples[iu].feature_offset;

				//outputCloud.cloud.samples[iu].attributes[1] = second_color_tag; //ocio
			}
			k++;

			}

		}

	}
}

void SimulationManager::assignStep(vector<int> indices_to_consider)
{
	// Assign step: move samples

	assignStepPosition(indices_to_consider);

	if( allow_tag_optimization ||
		(!is_2d && SGPScene == 1) )
		assignStepAttributeTag( indices_to_consider );

	if(!is_2d) {
		assignStepAttribute( indices_to_consider );
	}

	// Tim: disabled this old code so that we can use attribute assignment in 2D

	//if (!is_2d){
	//	assignStepAttribute(indices_to_consider);
	//	if (SGPScene == 1){
	//		// if tag attribute
	//		assignStepAttributeTag( indices_to_consider );
	//	}
	//}
}

void SimulationManager::assignStepPosition(vector<int> indices_to_consider)
{
	// Optimize samples position with gradient descent

	float prefix = (float)sqrt(PI/2.0);

	//omp_set_nested(1);
	#pragma omp parallel for
	for (int ii = 0; ii < indices_to_consider.size(); ii ++){

		int i = indices_to_consider[ii];

		Vector3f old_pos(1000, 1000, 1000);
		Vector3f new_pos(1, 1, 1);
		int k = 0;

		while ((new_pos - old_pos).norm() > assignStep_precision && k < assignStep_max_iterations){
			old_pos = outputCloud.cloud.samples[i].position;
			Vector3f new_val(0, 0, 0);

			//#pragma omp parallel for
			for (int jjj=0; jjj<closeBrushIndicesGlobal[i].size(); jjj++){
				int jj=closeBrushIndicesGlobal[i][jjj];
				int j=backgroundGrid.brush_indices[jj];

				float scaled_omega = omega * backgroundGrid.points[j].scale;
				float scaled_dividend = 1.0f / ( 2.0 * scaled_omega * scaled_omega);
				float scaled_prefix = prefix * ( -2.0 / scaled_omega);// / backgroundGrid.points[j].scale; DOF!!

				float local_exp_approx_thres = exp_approx_thres/scaled_dividend;

				Vector3f dist_vec_first_part = - backgroundGrid.points[j].matching_position - old_pos +  backgroundGrid.points[j].position;

				Vector3f scaled_new_val(0, 0, 0);

				/*
				MatrixXf dist_vec_first_part_matrix(3, closeInputsIndicesGlobal[jj].size());
				dist_vec_first_part_matrix = dist_vec_second_part_global[jj];
				dist_vec_first_part_matrix.colwise() += dist_vec_first_part;
				dist_vec_first_part_matrix = dist_vec_first_part_matrix * (VectorXf)(((((dist_vec_first_part_matrix.array().pow(2)).colwise().sum())* (-scaled_dividend)).array().exp()).transpose());
				scaled_new_val = scaled_new_val + dist_vec_first_part_matrix.rowwise().sum();
				*/
				
				//#pragma omp parallel for
				for (int dd=0; dd<closeInputsIndicesGlobal[jj].size(); dd++){
					int d=closeInputsIndicesGlobal[jj][dd];

					Vector3f dist_vec_second_part = backgroundGrid.points[j].inputCloud_positions[d];
					Vector3f dist_vec = dist_vec_second_part + dist_vec_first_part;
					float dist_value = dist_vec.dot(dist_vec);

					if (dist_value > local_exp_approx_thres)
					{
						continue;
					}

					float coeff = exp1(-dist_value * scaled_dividend);

					if (consider_attributes_in_position_optimization)
					{
						//float attributeCoeff = backgroundGrid.points[j].inputCloud_attributes[d].dot(outputCloud.cloud.samples[i].attributes[0]); coeff *= attributeCoeff;
						float attributeCoeff = inputCloud.cloud.samples[d].attributes[1].dot(outputCloud.cloud.samples[i].attributes[1]); 

						coeff *= attributeCoeff;
					}
					scaled_new_val = scaled_new_val + dist_vec * coeff;
				}
				
				if (!use_output_precomputation){
					//#pragma omp parallel for
					for (int ee=0; ee<closeOutputsIndicesGlobal[jj].size(); ee++){
						int e = closeOutputsIndicesGlobal[jj][ee];
						Vector3f dist_vec = old_pos - outputCloud.cloud.samples[e].position;
						float dist_value = dist_vec.dot(dist_vec);

						if (dist_value > local_exp_approx_thres){
							continue;
						}
						float coeff = exp1(-dist_value * scaled_dividend);
						if (consider_attributes_in_position_optimization){
							//float attributeCoeff = outputCloud.cloud.samples[e].attributes[0].dot(outputCloud.cloud.samples[i].attributes[0]); coeff *= attributeCoeff;
							float attributeCoeff = outputCloud.cloud.samples[e].attributes[1].dot(outputCloud.cloud.samples[i].attributes[1]); coeff *= attributeCoeff;
						}
						scaled_new_val = scaled_new_val + dist_vec * coeff;
					}
				}

				new_val = new_val + scaled_new_val * scaled_prefix;
			}

			new_pos = old_pos - assignStep_eps * new_val;
			if (is_2d ){ 
				new_pos[2] = 0;
			}
			
			outputCloud.cloud.samples[i].position = new_pos;
			k++;
		}

		
		if (k == 1){
			if (outputCloud.cloud.samples[i].last_k == 1){
				outputCloud.cloud.samples[i].last_k = -3; // -3 replicated in Sample.cpp
			}else if (outputCloud.cloud.samples[i].last_k == 0){
				outputCloud.cloud.samples[i].last_k = 0; 
			}else{
				outputCloud.cloud.samples[i].last_k +=1;
			}
		}else{
			outputCloud.cloud.samples[i].last_k = 1;
		}

	}
}

void SimulationManager::assignStepAttribute(vector<int> indices_to_consider)
{
	// Optimize samples first attribute with gradient descent

	float prefix = (float)sqrt(PI/2.0);

	#pragma omp parallel for
	for (int ii = 0; ii < indices_to_consider.size(); ii ++){

		int i = indices_to_consider[ii];

		Vector3f old_att(1000, 1000, 1000);
		Vector3f new_att = outputCloud.cloud.samples[i].attributes[0];
		int k = 0;

		while ((new_att - old_att).norm() > assignStepAttribute_precision && k < assignStepAttribute_max_iterations){
			old_att = new_att;
			Vector3f new_val(0, 0, 0);
			for (int jjj=0; jjj<closeBrushIndicesGlobal[i].size(); jjj++){
				int jj=closeBrushIndicesGlobal[i][jjj];
				int j=backgroundGrid.brush_indices[jj];

				float scaled_omega = omega * backgroundGrid.points[j].scale;
				float scaled_dividend = 1.0f / ( 2.0 * scaled_omega * scaled_omega);
				float scaled_prefix = prefix * scaled_omega;// / backgroundGrid.points[j].scale; DOF!!

				float local_exp_approx_thres = exp_approx_thres/scaled_dividend;

				Vector3f scaled_new_val(0, 0, 0);

				Vector3f dist_vec_first_part = - backgroundGrid.points[j].matching_position - outputCloud.cloud.samples[i].position +  backgroundGrid.points[j].position;

				for (int dd=0; dd<closeInputsIndicesGlobal[jj].size(); dd++){
					int d=closeInputsIndicesGlobal[jj][dd];

					Vector3f dist_vec = backgroundGrid.points[j].inputCloud_positions[d] + dist_vec_first_part;
					float dist_value = dist_vec.dot(dist_vec);

					if (dist_value > local_exp_approx_thres){
						continue;
					}
					float coeff = exp1(-dist_value * scaled_dividend);
					scaled_new_val = scaled_new_val - 2 * backgroundGrid.points[j].inputCloud_attributes[d] * coeff;
				}

				if (use_output_precomputation == 0){
					for (int ee=0; ee<closeOutputsIndicesGlobal[jj].size(); ee++){
						int e = closeOutputsIndicesGlobal[jj][ee];
						Vector3f dist_vec = outputCloud.cloud.samples[i].position - outputCloud.cloud.samples[e].position;

						float dist_value = dist_vec.dot(dist_vec);
					
						if (dist_value > local_exp_approx_thres){
							continue;
						}
						float coeff = exp1(-dist_value * scaled_dividend);
						scaled_new_val = scaled_new_val + outputCloud.cloud.samples[e].attributes[0] * coeff;
					}
				}

				new_val = new_val + scaled_new_val * scaled_prefix;
			}

			new_att = old_att - assignStepAttribute_eps * new_val;
			if (consider_attributes_in_position_optimization || true){ //ocioooo
				new_att.normalize();
			}
			outputCloud.cloud.samples[i].attributes[0] = new_att;


			k++;
		}
	}
}

void SimulationManager::assignStepAttributeTag(vector<int> indices_to_consider)
{
	// Optimize samples second attribute (tag) with gradient descent

	float prefix = (float)sqrt(PI/2.0);

	#pragma omp parallel for
	for (int ii = 0; ii < indices_to_consider.size(); ii ++){

		int i = indices_to_consider[ii];

		Vector3f old_att(1000, 1000, 1000);
		Vector3f new_att = outputCloud.cloud.samples[i].attributes[1];
		int k = 0;

		while ((new_att - old_att).norm() > assignStepAttribute_precision && k < assignStepAttribute_max_iterations){
			old_att = new_att;
			Vector3f new_val(0, 0, 0);
			for (int jjj=0; jjj<closeBrushIndicesGlobal[i].size(); jjj++){
				int jj=closeBrushIndicesGlobal[i][jjj];
				int j=backgroundGrid.brush_indices[jj];

				float scaled_omega = omega * backgroundGrid.points[j].scale;
				float scaled_dividend = 1.0f / ( 2.0 * scaled_omega * scaled_omega);
				float scaled_prefix = prefix * scaled_omega;// / backgroundGrid.points[j].scale; DOF!!

				float local_exp_approx_thres = exp_approx_thres/scaled_dividend;

				Vector3f scaled_new_val(0, 0, 0);

				Vector3f dist_vec_first_part = - backgroundGrid.points[j].matching_position - outputCloud.cloud.samples[i].position +  backgroundGrid.points[j].position;

				for (int dd=0; dd<closeInputsIndicesGlobal[jj].size(); dd++){
					int d=closeInputsIndicesGlobal[jj][dd];

					Vector3f dist_vec = backgroundGrid.points[j].inputCloud_positions[d] + dist_vec_first_part;
					float dist_value = dist_vec.dot(dist_vec);

					if (dist_value > local_exp_approx_thres){
						continue;
					}
					float coeff = exp1(-dist_value * scaled_dividend);
					scaled_new_val = scaled_new_val - 2 * inputCloud.cloud.samples[d].attributes[1] * coeff;
				}

				if (use_output_precomputation == 0){
					for (int ee=0; ee<closeOutputsIndicesGlobal[jj].size(); ee++){
						int e = closeOutputsIndicesGlobal[jj][ee];
						Vector3f dist_vec = outputCloud.cloud.samples[i].position - outputCloud.cloud.samples[e].position;

						float dist_value = dist_vec.dot(dist_vec);
					
						if (dist_value > local_exp_approx_thres){
							continue;
						}
						float coeff = exp1(-dist_value * scaled_dividend);
						scaled_new_val = scaled_new_val + outputCloud.cloud.samples[e].attributes[1] * coeff;
					}
				}

				new_val = new_val + scaled_new_val * scaled_prefix;
			}


			new_att = old_att - assignStepAttribute_eps * new_val;

			if (consider_attributes_in_position_optimization){
				new_att.normalize();
			}
			
			outputCloud.cloud.samples[i].attributes[1] = new_att;//.normalized();
			k++;
		}
	}
}


// --------------------------------- Generate indices in different ways -- START -------------------------------------------------

bool pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) 
{
  return firstElem.first < secondElem.first;
}

float SimulationManager::estimateOutputPointScale(int index)
{
	float total_scale = 0;
	for (int jjj=0; jjj<closeBrushIndicesGlobal[index].size(); jjj++){
		int jj=closeBrushIndicesGlobal[index][jjj];
		int j=backgroundGrid.brush_indices[jj];
		total_scale += backgroundGrid.points[j].scale;
	}
	total_scale = (float)total_scale/(float)closeBrushIndicesGlobal[index].size();
	return total_scale;
}

vector<int> SimulationManager::generateOutputIndicesBasedOnLowPointsDensity(vector<float> output_estimated_scales)
{
	vector<pair<int,int> > indPairs;
	for (int y = 0; y < outputCloud.cloud.samples.size(); y++){
		int myCount=0;
		for (int yy=0; yy<outputCloud.cloud.samples.size(); yy++){
			float dist = (outputCloud.cloud.samples[y].position - outputCloud.cloud.samples[yy].position).norm();
			if (dist < output_estimated_scales[y] * inputCloud.average_samples_distance * 2.5f){
				myCount++;
			}
		}
		indPairs.push_back(make_pair(myCount,y));
	}
	std::sort(indPairs.begin(), indPairs.end(), pairCompare);
	vector<int> myIndices;
	int pointsToCheck=300;
	for (int y=0; y<pointsToCheck; y++){
		if (y > outputCloud.cloud.samples.size()-1){
			break;
		}
		bool closeEnoughToStroke=false;
		for (int k=0;k<brush_line_points.size(); k++){
			float distanceFromStroke = (outputCloud.cloud.samples[indPairs[y].second].position - brush_line_points[k]).norm();
			if (distanceFromStroke < brushSize+neighSize){
				closeEnoughToStroke=true;
				break;
			}
		}
		if (!closeEnoughToStroke){
			pointsToCheck++;
		}
		else{
			myIndices.push_back(indPairs[y].second);
			myIndices.push_back(indPairs[y].second);
			myIndices.push_back(indPairs[y].second);
		}
	}

	return myIndices;
}
vector<int> SimulationManager::generateMatchingIndicesBasedOnAllIndices(){
	vector<int> matching_points_indices;
	for (int ii = 0; ii < backgroundGrid.brush_indices.size(); ii++){
		//int i = backgroundGrid.brush_indices[ii];
		matching_points_indices.push_back(ii);
	}

	return matching_points_indices;
}

vector<int> SimulationManager::generateOutputIndicesBasedOnAllIndices()
{
	vector<int> myIndices;
	for (int ii = 0; ii < outputCloud.cloud.samples.size(); ii ++){
		if (outputCloud.cloud.samples[ii].feature_id != -1){
			continue;
		}
		myIndices.push_back(ii);
	}
	
	return myIndices;
}

vector<int> SimulationManager::generateOutputIndicesBasedOnZeroK()
{
	vector<int> myIndices;
	for (int ii = 0; ii < outputCloud.cloud.samples.size(); ii ++){
		if (outputCloud.cloud.samples[ii].feature_id != -1){
			continue;
		}
		if (outputCloud.cloud.samples[ii].last_k == 0){
			myIndices.push_back(ii);
		}
	}
	
	return myIndices;
}

vector<int> SimulationManager::generateMatchingIndicesBasedOnNonZeroK()
{
	vector<int> matching_points_indices;
	for (int ii = 0; ii < backgroundGrid.brush_indices.size(); ii++){
		//int i = backgroundGrid.brush_indices[ii];

		for (int kk=0; kk<closeOutputsIndicesGlobal[ii].size(); kk++){
			int k=closeOutputsIndicesGlobal[ii][kk];
			if (outputCloud.cloud.samples[k].last_k != 0){
				matching_points_indices.push_back(ii);
				break;
			}
		}
	}

	return matching_points_indices;
}

vector<int> SimulationManager::generateMatchingIndicesBasedOnBrushPosition()
{
	float closestBrushLinePointIndex = 0;
	float minDist = 1000000000000;
	for (int i = 0; i < brush_line_points.size(); i ++){
		float newDist = (brush_line_points[i] - mousePos).norm();
		if (newDist < minDist){
			minDist = newDist;
			closestBrushLinePointIndex = i;
		}
	}

	vector<int> matching_points_indices;
	for (int ii = 0; ii < backgroundGrid.brush_indices.size(); ii++){
		int i = backgroundGrid.brush_indices[ii];

		float dist = (backgroundGrid.points[i].position - brush_line_points[closestBrushLinePointIndex]).norm();
		if (dist < 0*brushSize + 2*neighSize){
			matching_points_indices.push_back(ii);
		}
		
	}

	return matching_points_indices;
}

vector<int> SimulationManager::generateOutputIndicesBasedOnBrushPosition()
{
	float closestBrushLinePointIndex = 0;
	float minDist = 1000000000000;
	for (int i = 0; i < brush_line_points.size(); i ++){
		float newDist = (brush_line_points[i] - mousePos).norm();
		if (newDist < minDist){
			minDist = newDist;
			closestBrushLinePointIndex = i;
		}
	}

	vector<int> myIndices;
	for (int ii = 0; ii < outputCloud.cloud.samples.size(); ii ++){
		float dist = (outputCloud.cloud.samples[ii].position - brush_line_points[closestBrushLinePointIndex]).norm();
		if (dist < 0*brushSize + neighSize){
			myIndices.push_back(ii);
		}
	}
	
	return myIndices;
}

vector<int> SimulationManager::generateOutputIndicesBasedOnNonZeroK()
{
	vector<int> myIndices;
	for (int ii = 0; ii < outputCloud.cloud.samples.size(); ii ++){
		if (outputCloud.cloud.samples[ii].last_k != 0){
			myIndices.push_back(ii);
		}
	}
	
	return myIndices;
}

// --------------------------------- Generate indices in different ways -- END -------------------------------------------------


void SimulationManager::addFeatureAtPosition(Vector3f position, int feature_id)
{
	// Add all the samples belonging to the feature feature_id at position position

	float minDist = 10000000000;
	int closestBgIndex = 0;
	for (int i = 0; i < backgroundGrid.points.size(); i ++){
		float newDist = (backgroundGrid.points[i].position - position).norm();
		if (newDist < minDist){
			minDist = newDist;
			closestBgIndex = i;
		}
	}

	vector<int> new_indices;
	outputCloud.features_indices.push_back(new_indices);
	outputCloud.features_closest_bg_point.push_back(closestBgIndex);

	for (int ii = 0; ii < inputCloud.features_indices[feature_id].size(); ii ++){
		int i = inputCloud.features_indices[feature_id][ii];
		//Vector3f dist_vec = inputCloud.cloud.samples[i].position - inputCloud.cloud.samples[inputCloud.features_indices[feature_id][0]].position;
		Vector3f dist_vec = backgroundGrid.points[closestBgIndex].inputCloud_positions[i] - backgroundGrid.points[closestBgIndex].inputCloud_positions[inputCloud.features_indices[feature_id][0]];// inputCloud.cloud.samples[i].position - inputCloud.cloud.samples[inputCloud.features_indices[feature_id][0]].position;
		Vector3f new_position = dist_vec + position;
		//Vector3f new_normal = inputCloud.cloud.samples[i].attributes[0];
		Vector3f new_normal = backgroundGrid.points[closestBgIndex].inputCloud_attributes[i];// inputCloud.cloud.samples[i].attributes[0];
		outputCloud.addSampleWithNormal(new_position, new_normal);
		outputCloud.cloud.samples[outputCloud.cloud.samples.size() - 1].feature_id = feature_id;

		Vector3f input_offset = backgroundGrid.points[closestBgIndex].inputCloud_positions[i] - backgroundGrid.points[closestBgIndex].inputCloud_positions[inputCloud.features_indices[feature_id][0]];// inputCloud.cloud.samples[i].position - inputCloud.cloud.samples[inputCloud.features_indices[feature_id][0]].position;
		outputCloud.cloud.samples[outputCloud.cloud.samples.size() - 1].feature_offset = input_offset; 
	
		outputCloud.features_indices[outputCloud.features_indices.size() - 1].push_back(outputCloud.cloud.samples.size() - 1);
	}
}

void SimulationManager::addFeaturesStep()
{
	// Add features

	vector<Vector3f> new_features_candidates;
	vector<int> new_features_indices;
	vector<int> new_features_counter;

	for (int ii = 0; ii < backgroundGrid.brush_indices.size(); ii++){
		int i = backgroundGrid.brush_indices[ii];

		float distFromBorder = backgroundGrid.points[i].matching_position.norm(); 
		if (abs(distFromBorder - inputCloud.max_distance_allowed) < neighSize/2.0f){
			continue;
		}

		for (int frand = 0; frand < inputCloud.features_indices.size(); frand ++){
			int f = rand() %  inputCloud.features_indices.size();
			Vector3f dist_vec = backgroundGrid.points[i].matching_position - backgroundGrid.points[i].inputCloud_positions[inputCloud.features_indices[f][0]];// inputCloud.cloud.samples[inputCloud.features_indices[f][0]].position;
			float dist = dist_vec.norm();
			if (dist < neighSize){
				bool farEnough = true;
				Vector3f origin_output_pos = backgroundGrid.points[i].position - dist_vec;

				for (int ff = 0; ff < outputCloud.features_indices.size(); ff ++){
					Vector3f dist_vec_output = outputCloud.cloud.samples[outputCloud.features_indices[ff][0]].position - origin_output_pos;
					float norm_output = dist_vec_output.norm();
					if (norm_output < neighSize/1.0f){
						farEnough = false;
						break;
					}
				}

				if (farEnough){
					bool closeEnough = false;
					for (int h = 0; h < inputCloud.features_indices[f].size(); h ++){
						Vector3f output_pos = backgroundGrid.points[i].position - backgroundGrid.points[i].matching_position + backgroundGrid.points[i].inputCloud_positions[inputCloud.features_indices[f][h]];// inputCloud.cloud.samples[inputCloud.features_indices[f][h]].position;
						for (int u = 0; u < outputCloud.cloud.samples.size(); u ++){
							float out_dist = (output_pos - outputCloud.cloud.samples[u].position).norm();
							if (out_dist < inputCloud.average_samples_distance * 1.5f && outputCloud.cloud.samples[u].feature_id == -1){
								closeEnough = true;
								break;
							}
						}
					}

					if (closeEnough){

						bool merged = false;
						for (int y = 0; y < new_features_candidates.size(); y ++){
							float dist = (origin_output_pos - new_features_candidates[y]).norm();
							if (dist < neighSize/1.0f && f == new_features_indices[y]){
								new_features_counter[y] ++;
								merged = true;
								break;
							}
						}
						if (merged == false){
							new_features_candidates.push_back(origin_output_pos);
							new_features_indices.push_back(f);
							new_features_counter.push_back(1);
						}

						//addFeatureAtPosition(origin_output_pos, f);
						////automatic_add_features = 0;
						//return;
					}
				}

			}
		}
	}

	int max_count = 0;
	int max_index = 0;
	for (int u = 0; u < new_features_candidates.size(); u ++){
		if (new_features_counter[u] > max_count){
			max_count = new_features_counter[u];
			max_index = u;
		}
	}
	if (max_count > 2){
		addFeatureAtPosition(new_features_candidates[max_index], new_features_indices[max_index]);
	}
	
	//automatic_add_features = 0;		



}



void SimulationManager::addStep()
{
	// Add new samples

	if (outputCloud.cloud.samples.size() == 0){
		return;
	}

	int added_counter = 0;

	vector<float> output_estimated_scales;
	for (int y=0; y<outputCloud.cloud.samples.size(); y++){
		float estimated_scale = estimateOutputPointScale(y);
		output_estimated_scales.push_back(estimated_scale);
	}
	
	// Find active samples to use as seeds for new ones

	vector<int> myIndices;
	if (use_only_active_output_points == 0){
		myIndices = generateOutputIndicesBasedOnAllIndices();
	} else{
		//myIndices = generateOutputIndicesBasedOnNonZeroK();
		myIndices = generateOutputIndicesBasedOnBrushPosition();
	}


	for (int ii=0; ii<myIndices.size(); ii++){
		int i=myIndices[ii];

		// Ensure that the samples are not at the borders

		bool ok_to_add = false;
		if (outputDomain != 2){
			for (int b=0; b<brush_line_points.size(); b++){
				float dist = (outputCloud.cloud.samples[i].position - brush_line_points[b]).norm();
				float thres = brushSize + neighSize;
				if (0){
					thres = brushSize;
				}
				if (dist < thres){
					ok_to_add = true;
					break;
				}
			}
		} else if (outputDomain == 2){
			for (int b=0; b<backgroundGrid.points.size(); b++){
				float dist = (outputCloud.cloud.samples[i].position - backgroundGrid.points[b].position).norm();
				float thres = neighSize * 2; //no scale
				if (dist < thres){
					ok_to_add = true;
					break;
				}
			}
		}
		if (!ok_to_add){
			continue;
		}

		// Generate new random candidate position

		Vector3f new_pos = outputCloud.cloud.samples[i].position;
		Vector3f new_normal = outputCloud.cloud.samples[i].attributes[0];
		Vector3f rand_vec = Vector3f(rand()%100,rand()%100,rand()%100);
		Vector3f add_dir = rand_vec;// (new_normal).cross(rand_vec);
		add_dir = add_dir.normalized() * inputCloud.average_samples_distance * output_estimated_scales[i] * 2.0f;
		if (add_dir[0]!=add_dir[0] || add_dir[1]!=add_dir[1] || add_dir[2]!=add_dir[2]){
			continue;
		}
		if (rand() % 2){
			add_dir[0] *= -1;
		}
		if (rand() % 2){
			add_dir[1] *= -1;
		}
		if (rand() % 2){
			add_dir[2] *= -1;
		}
		new_pos = new_pos + add_dir;

		// Add the sample and update neighbors

		if (scene == 0 || scene == 10 || scene == 13 || SGPScene == 1){
			vector<Vector3f> attributes;

			Sample& sample_i = outputCloud.cloud.samples[i];
			attributes = sample_i.attributes;

			attributes[0] = new_normal;

			outputCloud.addSample(new_pos, attributes);

			Sample& sample = outputCloud.cloud.samples.back();
			sample.hackScale = sample_i.hackScale;
		}else{
			//cout << new_pos << " " << new_normal << " a"<< endl;
			outputCloud.addSampleWithNormal(new_pos, new_normal);
		}

		vector<int> output_points_neighbors_indices;
		vector<int> matching_points_neighbors_indices;
		vector<int> matching_points_neighbors_indices_toPass;
		vector<int> used_indices;
		vector<Vector3f> old_matching_point_positions;
		vector<Vector3f> old_output_point_positions;

		for (int rr=0; rr<backgroundGrid.brush_indices.size(); rr++){
			int r = backgroundGrid.brush_indices[rr];
			float dist = (new_pos - backgroundGrid.points[r].position).norm();
			
			if (dist < neighSize * backgroundGrid.points[r].scale){
				if (backgroundGrid.points[r].energy == 0){
					continue;
				}
				matching_points_neighbors_indices.push_back(rr);
				matching_points_neighbors_indices_toPass.push_back(r);
				closeOutputsIndicesGlobal[rr].push_back(outputCloud.cloud.samples.size()-1);
				used_indices.push_back(rr);
				old_matching_point_positions.push_back(backgroundGrid.points[r].matching_position);
			}
		}
		closeBrushIndicesGlobal.push_back(used_indices);

		// Check if the energy is smaller, otherwise remove the sample

		float energy_before = 0;
		float energy_after = 0;
		vector<float> new_energies;
		//vector<int> used_indices;

		for (int rrr=0; rrr<matching_points_neighbors_indices.size(); rrr++){
			int rr = matching_points_neighbors_indices[rrr];
			int r = backgroundGrid.brush_indices[rr];
			
			//closeOutputsIndicesGlobal[rr].push_back(outputCloud.cloud.samples.size()-1);
			float old_energy = backgroundGrid.points[r].energy;
			float additional_energy = 0;
			float new_energy = 0;

			if (accurate_add_delete == 0){
				additional_energy = computeEnergyGivenFromOnePoint(rr, outputCloud.cloud.samples.size() - 1);
				new_energy = old_energy + additional_energy;
			} else if (accurate_add_delete == 1){
				new_energy = computeEnergyForPoint(rr);
			}

			new_energies.push_back(new_energy);
			//used_indices.push_back(rr);
			energy_before += old_energy;
			energy_after += new_energy;
		}

		//closeBrushIndicesGlobal.push_back(used_indices);

		if (energy_after >= energy_before){
			for (int rr=0; rr<used_indices.size(); rr++){
				closeOutputsIndicesGlobal[used_indices[rr]].pop_back();
			}

			if (accurate_add_delete == 1){
				for (int rr = 0; rr < output_points_neighbors_indices.size(); rr ++){
					outputCloud.cloud.samples[output_points_neighbors_indices[rr]].position = old_output_point_positions[rr];
				}
				for (int rr = 0; rr < matching_points_neighbors_indices_toPass.size(); rr ++){
					backgroundGrid.points[matching_points_neighbors_indices_toPass[rr]].matching_position = old_matching_point_positions[rr];
				}
			}
			outputCloud.popSample();
			closeBrushIndicesGlobal.pop_back();
		}
		else{
			for (int rr=0; rr<used_indices.size(); rr++){
				int r = backgroundGrid.brush_indices[used_indices[rr]];
				backgroundGrid.points[r].energy = new_energies[rr];
			}
			added_counter++;
		}
	}
	
}

void SimulationManager::deleteFeaturesStep()
{
	// Delete all the samples belonging to one feature

	for (int ff = 0; ff < outputCloud.features_indices.size(); ff ++){
		//Vector3f dist_vec_output = outputCloud.cloud.samples[outputCloud.features_indices[ff][0]].position - origin_output_pos;
	}
}

void SimulationManager::deleteBadSamples()
{
	// Delete unnecessary samples which make the energy increase

	int deleted_counter = 0;

	vector<int> indices_to_delete;

	// Find active candidate samples

	vector<int> myIndices;
	if (use_only_active_output_points == 0){
		myIndices = generateOutputIndicesBasedOnAllIndices();
	} else{
		myIndices = generateOutputIndicesBasedOnNonZeroK();
	}
	//reverse(myIndices.begin(), myIndices.end());

	for (int ii = myIndices.size() - 1; ii > -1; ii --){
		int i = myIndices[ii];

		if (outputCloud.cloud.samples[i].feature_id != -1){
			continue;
		}

		float energy_before = 0;
		float energy_after = 0;
		vector<float> new_energies;

		vector<int> output_points_neighbors_indices;
		vector<Vector3f> old_matching_point_positions;
		vector<Vector3f> old_output_point_positions;
		
		// Compute the energy given by the current sample

		for (int jj = 0; jj < closeBrushIndicesGlobal[i].size(); jj ++){
			int rr = closeBrushIndicesGlobal[i][jj];
			int r = backgroundGrid.brush_indices[rr];

			float old_energy = backgroundGrid.points[r].energy;
			float additional_energy = 0;
			float new_energy = 0;

			if (accurate_add_delete == 0){
				additional_energy = computeEnergyGivenFromOnePoint(rr, i);
				new_energy = old_energy - additional_energy;
			}
			
			new_energies.push_back(new_energy);
			energy_before += old_energy;
			energy_after += new_energy;
		}
		if (energy_after > energy_before){
			
		}
		else{
			// Add the sample to the list of elements that need to be deleted

			indices_to_delete.push_back(i);
			for (int jj = 0; jj < closeBrushIndicesGlobal[i].size(); jj ++){
				int rr = closeBrushIndicesGlobal[i][jj];
				int r = backgroundGrid.brush_indices[rr];
				if (accurate_add_delete == 0){
					closeOutputsIndicesGlobal[rr].erase(remove(closeOutputsIndicesGlobal[rr].begin(), closeOutputsIndicesGlobal[rr].end(), i), closeOutputsIndicesGlobal[rr].end());
				}
				backgroundGrid.points[r].energy = new_energies[jj];
			}
			deleted_counter++;
		}

		// Tim Davison: 
		// - We can reach a point where the energy is really poor and we might end up deleting everything at once. 
		//   So, only delete at most a few samples, with a chance to fix them in the next optimizations step.
		// - Deleting samples is also going to change the energy and they may be coherent (or near) each other, therefore
		//   we would need to delete and then recalculate the neighborhoods. This is expensive, so defer deletion until
		//   we have cleaned up the output in the next optimization and maybe we don't need to delete.
		if(deleted_counter > 4)
			break;
	}

	// Remove all the elements in the list

	for (int u = 0; u < indices_to_delete.size(); u ++){
		
		outputCloud.deleteSample(indices_to_delete[u]);
		closeBrushIndicesGlobal.erase(closeBrushIndicesGlobal.begin() + indices_to_delete[u]);

		//update features indices
		for (int f = 0; f < outputCloud.features_indices.size(); f ++){
			for (int ff = 0; ff < outputCloud.features_indices[f].size(); ff ++){
				if ( outputCloud.features_indices[f][ff] > indices_to_delete[u] ){
					outputCloud.features_indices[f][ff] -= 1;
				}
			}
		}

	}
}
