#pragma once

#include "InputCloud.h"
#include "BackgroundGrid.h"
#include "OutputCloud.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif


class SimulationManager
{
public:
	SimulationManager();
	~SimulationManager(void);

	void addFeatureAtPosition( Eigen::Vector3f position, int feature_id);

	InputCloud inputCloud;
	BackgroundGrid backgroundGrid;
	OutputCloud outputCloud;

	// Nukes the inputCloud, backgroundGrid and outputCloud, brush_line_points, averageEnergies and discreteElementsAttributes.
	void clear();

	// Call init after configuring the inputCloud.
	void init();

	// Call init before calling run.
	void run();

	float disk_render_size;
	float output_canvas_size;
	float neighSize;
	float brushSize;
	float omega;

	bool is_2d;

	bool allow_tag_optimization; // i.e. attribute[1] optimization

	float assignStep_eps;
	float matchingStep_eps;
	float assignStepAttribute_eps;

	int draw_background_grid;
	int draw_background_energies;
	int draw_output_points;
	int draw_matching_points;

	/*
	Tim: we need to support other periodic border types, I'm changing this code.
	int use_periodic_borders;
	*/
	enum class PeriodicBorderType
	{
		None,
		Radius,
		Bounds
	};

	PeriodicBorderType periodic_border_type;

	Eigen::AlignedBox3f periodic_border_bounds;
	
	int allow_add_points;
	int allow_delete_points;

	bool allow_rotation;
	bool allow_scale;

	std::vector<Eigen::Vector3f> brush_line_points;

	int exp_approx_thres;

	int use_input_precomputation;
	int use_output_precomputation;

	std::vector<float> averageEnergies;

	int show_average_energy;

	int matching_position_new_seeds;

	int accurate_add_delete;

	int use_only_active_output_points;

	float assignStep_precision;
	float matchingStep_precision;
	float assignStepAttribute_precision;

	int assignStep_max_iterations;
	int matchingStep_max_iterations;
	int assignStepAttribute_max_iterations;

	int scene;

	int automatic_add_features;

	int consider_attributes_in_position_optimization;

	Eigen::Vector3f mousePos;

	std::vector<Eigen::Vector3f> discreteElementsAttributes;

	int show_mouse;

	bool is_paused;

	float factorTS;

	int automatic_multiscale_opt;


	int outputDomain;

protected:
	int SGPScene;

	Eigen::Vector3f first_color_tag;
	Eigen::Vector3f second_color_tag;

private:
	void computeNeighbors();
	void matchingStep(std::vector<int> indices_to_consider);
	void assignStep(std::vector<int> indices_to_consider);
	void addStep();
	void addFeaturesStep();
	void deleteFeaturesStep();
	void assignStepPosition(std::vector<int> indices_to_consider);
	void assignStepAttribute(std::vector<int> indices_to_consider);
	void assignStepAttributeTag(std::vector<int> indices_to_consider);
	void assignFeaturesStep();
	void computeEnergy();
	void setEnergyToPoint(int index);
	float computeEnergyForPoint(int index_brush_indices);
	void deleteBadSamples();
	float computeEnergyGivenFromOnePoint(int bg_point_index, int output_point_index);
	float estimateOutputPointScale(int index);
	void computeAverageEnergy();
	std::vector<int> findInputNeighbors(int ii);
	Eigen::Vector3f pickNewRandomMatchingPosition(int bg_point_index);
	std::vector<int> generateOutputIndicesBasedOnLowPointsDensity(std::vector<float> output_estimated_scales);
	std::vector<int> generateOutputIndicesBasedOnAllIndices();
	std::vector<int> generateOutputIndicesBasedOnZeroK();
	std::vector<int> generateOutputIndicesBasedOnNonZeroK();
	std::vector<int> generateMatchingIndicesBasedOnAllIndices();
	std::vector<int> generateMatchingIndicesBasedOnNonZeroK();
	std::vector<int> generateOutputIndicesBasedOnBrushPosition();
	std::vector<int> generateMatchingIndicesBasedOnBrushPosition();
	void autoExportCloud();
	void matchingDTE();
	void assignDTE();
	void computeNeighborsDTE();
	void preprocessInput();
	void multiScaleAdapt();

	std::vector<std::vector<int> > closeOutputsIndicesGlobal;
	std::vector<std::vector<int> > closeInputsIndicesGlobal;
	std::vector<std::vector<int> > closeBrushIndicesGlobal;


	 
	std::vector<Eigen::MatrixXf> dist_vec_second_part_global;

	std::vector<std::vector<int> > outputPreProcess_indices_global;
	std::vector<std::vector<int> > outputPreProcess_coefficients_global;

	int sim_step;
	bool autoExport;

	int frameExportStep;
	int sigmaExportStep;

	int automatic_delete_features;

	std::vector<std::vector<int> > neighborsDTE;
	std::vector<std::vector<int> > dteBigArray;
	std::vector<int> dteBigMatchings;

	std::vector<std::vector<int> > kco_neighbors;

	std::vector<float> times_iterations;
	

};

