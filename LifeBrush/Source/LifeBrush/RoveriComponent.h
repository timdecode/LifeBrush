// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "Components/ActorComponent.h"

#include "RegionGrowingComponent.h"

#include "Algorithm/Algorithm_Roveri.h"

#include "RoveriComponent.generated.h"




UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API URoveriComponent : public URegionGrowingComponent
{
	GENERATED_BODY()

public:
	virtual Algorithm* getOrCreateAlgorithm() override;
	virtual void loadParameters() override;
	virtual void initAlgorithms() override;

	virtual FString parametersString() override;

	virtual void didClear();

protected:
	virtual void _mainThreadGenerationWorkDone() override;


public:
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float intputCloud_target_average_distance = 8.5f;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float background_grid_overlap_factor = 1.65f;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool allow_tag_optimization = true; // TIM: this has to be true to use element types

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float disk_render_size = 0.65f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float output_canvas_size = 750.0f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float neighSize = 110.0f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float roveriBrushSize = 100.0f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float omega = 1.0f;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool is_2d = true;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float assignStep_eps = 0.031f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float matchingStep_eps = 0.031f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float assignStepAttribute_eps = 0.0015f;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int draw_background_grid = 1;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int draw_background_energies = 0;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int draw_output_points = 1;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int draw_matching_points = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int use_periodic_borders = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int allow_add_points = 1;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int allow_delete_points = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool allow_rotation = 0;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool allow_scale = 0;

	std::vector<Eigen::Vector3f> brush_line_points;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int exp_approx_thres = 10;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int use_input_precomputation = 0;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int use_output_precomputation = 0;

	std::vector<float> averageEnergies;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int show_average_energy = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int matching_position_new_seeds = 3;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int accurate_add_delete = 0;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int use_only_active_output_points = 0;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float assignStep_precision = 0.005f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float matchingStep_precision = 0.005f;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float assignStepAttribute_precision = 0.005f;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int assignStep_max_iterations = 10;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int matchingStep_max_iterations = 10;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int assignStepAttribute_max_iterations = 5;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int scene = 2;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int automatic_add_features = 0;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int consider_attributes_in_position_optimization  = 0;

	Eigen::Vector3f mousePos;

	std::vector<Eigen::Vector3f> discreteElementsAttributes;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int show_mouse = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool is_paused = false;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float factorTS = 1;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int automatic_multiscale_opt = 0;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) int outputDomain;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool debug_drawOutputPoints = false;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) bool debug_drawOutputMatchingPoints = false;

	UPROPERTY(EditAnywhere, Category = "Roveri-Simulation") bool debug_drawMatchingPointsHistory = false;
	UPROPERTY(EditAnywhere, Category = "Roveri-Simulation") int debug_drawOutputMatchingPoints_numMatchingPoints = 10;


	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) UStaticMesh * debug_outputPointMesh = nullptr;
	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) UMaterialInterface * debug_outputPointMaterial = nullptr;

	UPROPERTY( EditAnywhere, Category = "Roveri-Simulation" ) float debug_pointScale = 0.1f;

	UPROPERTY( EditAnywhere, Category = "Tim" ) ERoveriPeriodicBounds periodicBoundsMode = ERoveriPeriodicBounds::Radius;


protected:
	std::unique_ptr<Algorithm_Roveri> _algorithmRoveri;

	UInstancedStaticMeshComponent * _debug_outputPointsISMC = nullptr;
	UInstancedStaticMeshComponent * _debug_outputMatchingPointsISMC = nullptr;
	UInstancedStaticMeshComponent * _debug_outputMatchingPointsHistoryISMC = nullptr;


protected:
	UInstancedStaticMeshComponent& debug_outputPointsMesh();	// will always return an ISMC
	UInstancedStaticMeshComponent& debug_outputMatchingPointsMesh();
	UInstancedStaticMeshComponent& debug_outputMatchingPointsHistoryMesh();

	void updateDebugOutputPointsMesh();
	void updateDebugOutputMatchingPointsMesh();
	void updateDebugOutputMatchingPointsHistoryMesh();
};

