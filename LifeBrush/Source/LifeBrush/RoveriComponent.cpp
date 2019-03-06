// Fill out your copyright notice in the Description page of Project Settings.

#include "LifeBrush.h"
#include "RoveriComponent.h"

#include "Algorithm/Algorithm_Roveri.h"


Algorithm* URoveriComponent::getOrCreateAlgorithm()
{
	if (!_algorithmRoveri)
		_algorithmRoveri = std::make_unique<Algorithm_Roveri>(*_synthesisContext);

	return _algorithmRoveri.get();
}

void URoveriComponent::loadParameters()
{
	Super::loadParameters();


	_algorithmRoveri->simulation.inputCloud.target_average_distance = intputCloud_target_average_distance;

	_algorithmRoveri->simulation.backgroundGrid.overlap_factor = background_grid_overlap_factor;

	_algorithmRoveri->simulation.allow_tag_optimization = allow_tag_optimization;

	_algorithmRoveri->simulation.disk_render_size = disk_render_size;
	_algorithmRoveri->simulation.output_canvas_size = output_canvas_size;
	_algorithmRoveri->simulation.neighSize = neighSize;
	_algorithmRoveri->simulation.brushSize = brushSize;
	_algorithmRoveri->simulation.omega = omega;

	_algorithmRoveri->simulation.is_2d = is_2d;

	_algorithmRoveri->simulation.assignStep_eps = assignStep_eps;
	_algorithmRoveri->simulation.matchingStep_eps = matchingStep_eps;
	_algorithmRoveri->simulation.assignStepAttribute_eps = assignStepAttribute_eps;

	_algorithmRoveri->simulation.assignStep_max_iterations = assignStep_max_iterations;
	_algorithmRoveri->simulation.matchingStep_max_iterations = matchingStep_max_iterations;
	_algorithmRoveri->simulation.assignStepAttribute_max_iterations = assignStepAttribute_max_iterations;

	_algorithmRoveri->simulation.draw_background_grid = draw_background_grid;
	_algorithmRoveri->simulation.draw_background_energies = draw_background_energies;
	_algorithmRoveri->simulation.draw_output_points = draw_output_points;
	_algorithmRoveri->simulation.draw_matching_points = draw_matching_points;

	if (use_periodic_borders)
	{
		if (periodicBoundsMode == ERoveriPeriodicBounds::Radius)
			_algorithmRoveri->simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Radius;
		else if(periodicBoundsMode == ERoveriPeriodicBounds::BoundsInset)
			_algorithmRoveri->simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Bounds;
		else if(periodicBoundsMode == ERoveriPeriodicBounds::GenerativeLabel)
			_algorithmRoveri->simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Bounds;

		_algorithmRoveri->periodicBoundsMode = periodicBoundsMode;
	}
	else
		_algorithmRoveri->simulation.periodic_border_type = SimulationManager::PeriodicBorderType::None;

	_algorithmRoveri->simulation.allow_add_points = allow_add_points;
	_algorithmRoveri->simulation.allow_delete_points = allow_delete_points;

	_algorithmRoveri->simulation.allow_rotation = allow_rotation;
	_algorithmRoveri->simulation.allow_scale = allow_scale;

	_algorithmRoveri->simulation.brush_line_points = brush_line_points;

	_algorithmRoveri->simulation.exp_approx_thres = exp_approx_thres;

	_algorithmRoveri->simulation.use_input_precomputation = use_input_precomputation;
	_algorithmRoveri->simulation.use_output_precomputation = use_output_precomputation;

	_algorithmRoveri->simulation.averageEnergies =  averageEnergies;

	_algorithmRoveri->simulation.show_average_energy = show_average_energy;

	_algorithmRoveri->simulation.matching_position_new_seeds = matching_position_new_seeds;

	_algorithmRoveri->simulation.accurate_add_delete = accurate_add_delete;

	_algorithmRoveri->simulation.use_only_active_output_points = use_only_active_output_points;

	_algorithmRoveri->simulation.assignStep_precision = assignStep_precision;
	_algorithmRoveri->simulation.matchingStep_precision = matchingStep_precision;
	_algorithmRoveri->simulation.assignStepAttribute_precision = assignStepAttribute_precision;

	_algorithmRoveri->simulation.scene = scene;

	_algorithmRoveri->simulation.automatic_add_features = automatic_add_features;

	_algorithmRoveri->simulation.consider_attributes_in_position_optimization = consider_attributes_in_position_optimization;

	_algorithmRoveri->simulation.mousePos = mousePos;

	_algorithmRoveri->simulation.discreteElementsAttributes = discreteElementsAttributes;

	_algorithmRoveri->simulation.show_mouse = show_mouse;

	_algorithmRoveri->simulation.is_paused = is_paused;

	_algorithmRoveri->simulation.factorTS = factorTS;

	_algorithmRoveri->simulation.automatic_multiscale_opt = automatic_multiscale_opt;

	_algorithmRoveri->simulation.outputDomain = outputDomain;

	_algorithmRoveri->periodicBoundsMode = periodicBoundsMode;
}

void URoveriComponent::initAlgorithms()
{
	Super::initAlgorithms();

	_algorithmRoveri->init( _synthesisContext->threadPool );
}

FString URoveriComponent::parametersString()
{
	FString params = Super::parametersString();

	params += FString::Printf( TEXT( "inputCloud.target_average_distance %f\n" ), _algorithmRoveri->simulation.inputCloud.target_average_distance );

	params += FString::Printf( TEXT( "backgroundGrid.overlap_factor %f\n" ), _algorithmRoveri->simulation.backgroundGrid.overlap_factor );

	params += FString::Printf( TEXT( "allow_tag_optimization %d\n" ), _algorithmRoveri->simulation.allow_tag_optimization );

	params += FString::Printf( TEXT( "disk_render_size %f\n" ), _algorithmRoveri->simulation.disk_render_size );
	params += FString::Printf( TEXT( "output_canvas_size %f\n" ), _algorithmRoveri->simulation.output_canvas_size );
	params += FString::Printf( TEXT( "neighSize %f\n" ), _algorithmRoveri->simulation.neighSize );
	params += FString::Printf( TEXT( "brushSize %f\n" ), _algorithmRoveri->simulation.brushSize );
	params += FString::Printf( TEXT( "omega %f\n" ), _algorithmRoveri->simulation.omega );

	params += FString::Printf( TEXT( "is_2d %d\n" ), _algorithmRoveri->simulation.is_2d );

	params += FString::Printf( TEXT( "assignStep_eps %f\n" ), _algorithmRoveri->simulation.assignStep_eps );
	params += FString::Printf( TEXT( "matchingStep_eps %f\n" ), _algorithmRoveri->simulation.matchingStep_eps );
	params += FString::Printf( TEXT( "assignStepAttribute_eps %f\n" ), _algorithmRoveri->simulation.assignStepAttribute_eps );

	params += FString::Printf( TEXT( "assignStep_max_iterations %dn" ), _algorithmRoveri->simulation.assignStep_max_iterations );
	params += FString::Printf( TEXT( "matchingStep_max_iterations %d\n" ), _algorithmRoveri->simulation.matchingStep_max_iterations );
	params += FString::Printf( TEXT( "assignStepAttribute_max_iterations %\n" ), _algorithmRoveri->simulation.assignStepAttribute_max_iterations );

	params += FString::Printf( TEXT( "draw_background_grid %d\n" ), _algorithmRoveri->simulation.draw_background_grid );
	params += FString::Printf( TEXT( "draw_background_energies %d\n" ), _algorithmRoveri->simulation.draw_background_energies );
	params += FString::Printf( TEXT( "draw_output_points %d\n" ), _algorithmRoveri->simulation.draw_output_points );
	params += FString::Printf( TEXT( "draw_matching_points %d\n" ), _algorithmRoveri->simulation.draw_matching_points );

	//params += FString::Printf( TEXT( "use_periodic_borders %d\n" ), _algorithmRoveri->simulation.use_periodic_borders );
	if (_algorithmRoveri->simulation.periodic_border_type == SimulationManager::PeriodicBorderType::None)
		params += FString::Printf(TEXT("periodic_border_type None\n"));
	else if (_algorithmRoveri->simulation.periodic_border_type == SimulationManager::PeriodicBorderType::Radius)
		params += FString::Printf(TEXT("periodic_border_type Radius\n"));
	else if (_algorithmRoveri->simulation.periodic_border_type == SimulationManager::PeriodicBorderType::Bounds)
		params += FString::Printf(TEXT("periodic_border_type Bounds\n"));

	params += FString::Printf( TEXT( "allow_add_points %d\n" ), _algorithmRoveri->simulation.allow_add_points );
	params += FString::Printf( TEXT( "allow_delete_points %d\n" ), _algorithmRoveri->simulation.allow_delete_points );

	params += FString::Printf( TEXT( "allow_rotation %d\n" ), _algorithmRoveri->simulation.allow_rotation );
	params += FString::Printf( TEXT( "allow_scale %d\n" ), _algorithmRoveri->simulation.allow_scale );

	params += FString::Printf( TEXT( "exp_approx_thres %f\n" ), _algorithmRoveri->simulation.exp_approx_thres );

	params += FString::Printf( TEXT( "use_input_precomputation %d\n" ), _algorithmRoveri->simulation.use_input_precomputation );
	params += FString::Printf( TEXT( "use_output_precomputation %d\n" ), _algorithmRoveri->simulation.use_output_precomputation );

	params += FString::Printf( TEXT( "show_average_energy %d\n" ), _algorithmRoveri->simulation.show_average_energy );

	params += FString::Printf( TEXT( "matching_position_new_seeds %d\n" ), _algorithmRoveri->simulation.matching_position_new_seeds );

	params += FString::Printf( TEXT( "accurate_add_delete %d\n" ), _algorithmRoveri->simulation.accurate_add_delete );

	params += FString::Printf( TEXT( "use_only_active_output_points %d\n" ), _algorithmRoveri->simulation.use_only_active_output_points );

	params += FString::Printf( TEXT( "assignStep_precision %f\n" ), _algorithmRoveri->simulation.assignStep_precision );
	params += FString::Printf( TEXT( "matchingStep_precision %f\n" ), _algorithmRoveri->simulation.matchingStep_precision );
	params += FString::Printf( TEXT( "assignStepAttribute_precision %f\n" ), _algorithmRoveri->simulation.assignStepAttribute_precision );

	params += FString::Printf( TEXT( "scene %d\n" ), _algorithmRoveri->simulation.scene );

	params += FString::Printf( TEXT( "automatic_add_features %d\n" ), _algorithmRoveri->simulation.automatic_add_features );

	params += FString::Printf( TEXT( "consider_attributes_in_position_optimization %d\n" ), _algorithmRoveri->simulation.consider_attributes_in_position_optimization );



	params += FString::Printf( TEXT( "show_mouse %d\n" ), _algorithmRoveri->simulation.show_mouse );

	params += FString::Printf( TEXT( "is_paused %d\n" ), _algorithmRoveri->simulation.is_paused );

	params += FString::Printf( TEXT( "factorTS %f\n" ), _algorithmRoveri->simulation.factorTS );

	params += FString::Printf( TEXT( "automatic_multiscale_opt %d\n" ), _algorithmRoveri->simulation.automatic_multiscale_opt );

	params += FString::Printf( TEXT( "outputDomain %d\n" ), _algorithmRoveri->simulation.outputDomain );

	// periodicBoundsMode:
	{
		const UEnum * RoveriPeriodicBounds = FindObject<UEnum>(ANY_PACKAGE, TEXT("ERoveriPeriodicBounds"));
		int32 int32_periodicBoundsMode = (int32)generationMode;

		params += FString::Printf(TEXT("periodicBoundsMode: %s\n"), *(RoveriPeriodicBounds->GetDisplayNameTextByValue(int32_periodicBoundsMode).ToString()));
	}


	return params;
}

void URoveriComponent::_mainThreadGenerationWorkDone()
{
	updateDebugOutputPointsMesh();
	updateDebugOutputMatchingPointsMesh();
	updateDebugOutputMatchingPointsHistoryMesh();
}

void URoveriComponent::didClear()
{
	debug_outputMatchingPointsHistoryMesh().ClearInstances();
}

UInstancedStaticMeshComponent& URoveriComponent::debug_outputPointsMesh()
{
	if(_debug_outputPointsISMC == nullptr)
	{
		_debug_outputPointsISMC = NewObject<UInstancedStaticMeshComponent>( GetOwner() );

		_debug_outputPointsISMC->SetStaticMesh( debug_outputPointMesh );
		_debug_outputPointsISMC->SetMobility( EComponentMobility::Static );
		_debug_outputPointsISMC->SetCollisionProfileName( TEXT( "NoCollision" ) );
		_debug_outputPointsISMC->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
		_debug_outputPointsISMC->SetMaterial( 0, debug_outputPointMaterial );

		_debug_outputPointsISMC->RegisterComponent();
	}

	return *_debug_outputPointsISMC;
}

UInstancedStaticMeshComponent& URoveriComponent::debug_outputMatchingPointsMesh()
{
	if(_debug_outputMatchingPointsISMC == nullptr)
	{
		_debug_outputMatchingPointsISMC = NewObject<UInstancedStaticMeshComponent>( GetOwner() );

		_debug_outputMatchingPointsISMC->SetStaticMesh( debug_outputPointMesh );
		_debug_outputMatchingPointsISMC->SetMobility( EComponentMobility::Static );
		_debug_outputMatchingPointsISMC->SetCollisionProfileName( TEXT( "NoCollision" ) );
		_debug_outputMatchingPointsISMC->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
		_debug_outputMatchingPointsISMC->SetMaterial( 0, debug_outputPointMaterial );

		_debug_outputMatchingPointsISMC->RegisterComponent();
	}

	return *_debug_outputMatchingPointsISMC;
}

UInstancedStaticMeshComponent& URoveriComponent::debug_outputMatchingPointsHistoryMesh()
{
	if (_debug_outputMatchingPointsHistoryISMC == nullptr)
	{
		_debug_outputMatchingPointsHistoryISMC = NewObject<UInstancedStaticMeshComponent>(GetOwner());

		_debug_outputMatchingPointsHistoryISMC->SetStaticMesh(debug_outputPointMesh);
		_debug_outputMatchingPointsHistoryISMC->SetMobility(EComponentMobility::Static);
		_debug_outputMatchingPointsHistoryISMC->SetCollisionProfileName(TEXT("NoCollision"));
		_debug_outputMatchingPointsHistoryISMC->AttachToComponent(this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);
		_debug_outputMatchingPointsHistoryISMC->SetMaterial(0, debug_outputPointMaterial);

		_debug_outputMatchingPointsHistoryISMC->RegisterComponent();
	}

	return *_debug_outputMatchingPointsHistoryISMC;
}

void URoveriComponent::updateDebugOutputPointsMesh()
{
	UInstancedStaticMeshComponent& ismc = debug_outputPointsMesh();

	ismc.ClearInstances();

	if(debug_drawOutputPoints)
	{
		SimulationManager& sim = _algorithmRoveri->simulation;

		for(BackgroundPoint& point : sim.backgroundGrid.points)
		{
			FTransform transform( FQuat::Identity, unreal(point.position), FVector( debug_pointScale ) );

			ismc.AddInstanceWorldSpace( transform );
		}
	}
}

void URoveriComponent::updateDebugOutputMatchingPointsMesh()
{
	UInstancedStaticMeshComponent& ismc = debug_outputMatchingPointsMesh();

	ismc.ClearInstances();

	if(debug_drawOutputMatchingPoints)
	{
		SimulationManager& sim = _algorithmRoveri->simulation;
		Eigen::Vector3f offset = _algorithmRoveri->matchingPointsExemplarOffset();

		for(BackgroundPoint& point : sim.backgroundGrid.points)
		{
			FTransform transform( FQuat::Identity, unreal(point.matching_position + offset), FVector( debug_pointScale ) );

			ismc.AddInstanceWorldSpace( transform );
		}
	}
}

void URoveriComponent::updateDebugOutputMatchingPointsHistoryMesh()
{
	UInstancedStaticMeshComponent& ismc = debug_outputMatchingPointsHistoryMesh();

	if (debug_drawMatchingPointsHistory)
	{
		SimulationManager& sim = _algorithmRoveri->simulation;
		Eigen::Vector3f offset = _algorithmRoveri->matchingPointsExemplarOffset();

		// enumerate the centre points
		size_t size = std::sqrt(sim.backgroundGrid.points.size());

		size_t mid = size / 2;

		size_t numPointsMid = std::sqrt(debug_drawOutputMatchingPoints_numMatchingPoints) / 2;
		
		for (int x = mid - numPointsMid; x < mid + numPointsMid; ++x)
		{
			for (int y = mid - numPointsMid; y < mid + numPointsMid; ++y)
			{
				size_t i = x + y * size;

				BackgroundPoint& point = sim.backgroundGrid.points[i];

				FTransform transform(FQuat::Identity, unreal(point.matching_position + offset), FVector(debug_pointScale));

				ismc.AddInstanceWorldSpace(transform);
			}
		}


		//for( size_t i = 0; i < debug_drawOutputMatchingPoints_numMatchingPoints && i < sim.backgroundGrid.points.size(); ++i )
		//{
		//	BackgroundPoint& point = sim.backgroundGrid.points[i];

		//	FTransform transform(FQuat::Identity, unreal(point.matching_position + offset), FVector(debug_pointScale));

		//	ismc.AddInstanceWorldSpace(transform);
		//}
	}
}
