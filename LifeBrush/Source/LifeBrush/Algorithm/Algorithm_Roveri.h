//
//  Algorith_Roveri.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2018-05-24.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//
//  This is a wrapper around Roveri et al.'s SimulationManager from their 2015 Eurographics paper, Example-Based Repetitive Structure Synthesis

#pragma once

#include "Algorithm/Algorithm.h"
#include "Algorithm/SimpleUniformGrid.h"

#include "Algorithm/Roveri/SimulationManager.h"

#include <random>
#include <iostream>

#include "Algorithm_Roveri.generated.h"

UENUM(BlueprintType)
enum class ERoveriPeriodicBounds : uint8
{
	Radius UMETA(DisplayName = "Radius"),
	BoundsInset UMETA(DisplayName = "BoundsInset"),
	GenerativeLabel UMETA(DisplayName = "GenerativeLabel"),
};

class Algorithm_Roveri : public Algorithm
{
public:
	Algorithm_Roveri(SynthesisContext& context) : Algorithm(context) {}
	virtual ~Algorithm_Roveri() {}


	virtual AlgorithmResult generate( std::vector<PositionFace>& positions, float radius = -1.0f, AABB limits = AABB() );


	virtual void clear() override;

	virtual void loadExemplar() override;

	virtual void endRound( AlgorithmResult& result );

	Eigen::Vector3f matchingPointsExemplarOffset() { return _matchingPointsOffset; }


public:
	// Exposed publicly to set parameters on the simulation. However, don't mess with the inputCloud or backgroundGrid directly, that's
	// what this wrapper does internally.
	SimulationManager simulation;

	ERoveriPeriodicBounds periodicBoundsMode = ERoveriPeriodicBounds::Radius;

protected:
	virtual void _initialize() override;

	void _initializeMatchingPointExemplarOffset();

	int16 _readType(Sample& sample);
	void _writeType( Sample& sample, int16 type );
	void _calculateMaxType();

	void _seeding( std::vector<PositionFace>& positions, AABB limits );
	
	// sample helpers
	std::vector< Eigen::Vector3f > _dummyAttributes();

	// background grid helpers
	void _preComputeInputEnergiesCallBack();

	void _generatePeriodicBorders_radius();
	void _generatePeriodicBorders_boundsInset();

	void _generatePeriodicBorders_generativeLabel();
	void _generateBackgroundGridFromNeighSize();

	void _generatedBackgroundGridBrushIndices_volume();
	void _generateBackgroundGridOrientations_volume();

	void _generateBackgroundGridMatchingPositions();


	void _scaleInputAndMatchingPointWithBrushScale( float scale_value, int bg_point_index );
	void _alignInputAndMatchingPointWithBrushDirection3d( Eigen::Vector3f dir_vector, int bg_point_index );

protected:
	int16 _maxType;

	Eigen::Vector3f _matchingPointsOffset;
};