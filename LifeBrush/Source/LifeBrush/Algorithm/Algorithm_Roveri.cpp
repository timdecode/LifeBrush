//
//  Algorith_Roveri.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2018-05-24.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.


#include "LifeBrush.h"

#include "Algorithm/Algorithm_Roveri.h"

#include "Algorithm/Roveri/BackgroundPoint.h"

#define _USE_MATH_DEFINES // for C++  

#include <cmath>  
#include <vector>
#include "Utility.h"

using namespace Eigen;
using namespace std;

void Algorithm_Roveri::_initialize()
{
	Algorithm::_initialize();

	simulation.clear();

	_calculateMaxType();

	InputCloud& inputCloud = simulation.inputCloud;
	{
		_initializeMatchingPointExemplarOffset();

		for(auto& element : _exemplar)
		{
			auto& eNode = graph().node(element.nodeHandle());

			auto attributes = _dummyAttributes();

			Vector3f position = eigen(eNode.position);

			inputCloud.cloud.samples.emplace_back( position, attributes );

			auto& sample = inputCloud.cloud.samples.back();
			sample.hackScale = eNode.scale;
			_writeType( sample, element.type );
		}

		inputCloud.initCloud_automatic();

		if (periodicBoundsMode == ERoveriPeriodicBounds::Radius)
		{
			simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Radius;

			_generatePeriodicBorders_radius();
		}
		else if (periodicBoundsMode == ERoveriPeriodicBounds::BoundsInset)
		{
			simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Bounds;
			
			_generatePeriodicBorders_boundsInset();
		}
		else if (periodicBoundsMode == ERoveriPeriodicBounds::GenerativeLabel)
		{
			simulation.periodic_border_type = SimulationManager::PeriodicBorderType::Bounds;

			_generatePeriodicBorders_generativeLabel();
		}

		_generateBackgroundGridFromNeighSize();

		if(simulation.use_input_precomputation == 1) {
			_preComputeInputEnergiesCallBack();
		}
	}

	simulation.init();

	//_generateBackgroundGrid_volume();

	_generatedBackgroundGridBrushIndices_volume();

	_generateBackgroundGridOrientations_volume();
}

void Algorithm_Roveri::_initializeMatchingPointExemplarOffset()
{
	_matchingPointsOffset = Vector3f::Zero();

	if (_exemplar.size() == 0)
		return;

	FGraphNode& elementNode0 = _exemplar.elementAt(0).node(_exemplar.graph);
	AlignedBox3f box(eigen(elementNode0.position));

	for (auto& element : _exemplar)
	{
		FGraphNode& elementNode = element.node(_exemplar.graph);

		box.extend(eigen(elementNode.position));
	}

	_matchingPointsOffset = box.center();
}



void Algorithm_Roveri::_calculateMaxType()
{
	int16 maxType = 0;

	for(auto& ePtr : _exemplar)
	{
		if(ePtr.type > maxType)
			maxType = ePtr.type;
	}

	_maxType = maxType;
}

void Algorithm_Roveri::_preComputeInputEnergiesCallBack()
{
	simulation.inputCloud.precomputeInputEnergy( simulation.omega );
}

void Algorithm_Roveri::_generatePeriodicBorders_radius()
{
	float maxDistance = 0;
	for(int i = 0; i<simulation.inputCloud.cloud.samples.size(); i++) 
	{
		float dist = (simulation.inputCloud.cloud.samples[i].position - simulation.inputCloud.center_of_mass).norm();
		if(dist > maxDistance) {
			maxDistance = dist;
		}
	}
	simulation.inputCloud.max_distance_allowed = maxDistance - simulation.neighSize * 1.2f;

	simulation.inputCloud.insideBorders = std::vector<int> ( simulation.inputCloud.cloud.samples.size(), 1 );

	for(int i = 0; i < simulation.inputCloud.cloud.samples.size(); i++) 
	{
		float dist = (simulation.inputCloud.cloud.samples[i].position - simulation.inputCloud.center_of_mass).norm();

		if(dist > maxDistance - simulation.neighSize) {
			simulation.inputCloud.insideBorders[i] = 0;
		}
	}
}


void Algorithm_Roveri::_generatePeriodicBorders_boundsInset()
{
	AlignedBox3f box = simulation.inputCloud.samplesBounds();

	const float size = simulation.neighSize * .75;
	const Vector3f inset(size, size, simulation.is_2d ? 0.0f : size);

	box.min() += inset;
	box.max() -= inset;

	if (simulation.is_2d)
	{
		box.min().z() -= 1.0f;
		box.max().z() += 1.0f;
	}

	simulation.periodic_border_bounds = box;

	simulation.inputCloud.insideBorders = vector<int>(simulation.inputCloud.cloud.samples.size(), 1);

	for (int i = 0; i < simulation.inputCloud.cloud.samples.size(); i++)
	{
		float dist = (simulation.inputCloud.cloud.samples[i].position - simulation.inputCloud.center_of_mass).norm();

		simulation.inputCloud.insideBorders[i] = box.contains(simulation.inputCloud.cloud.samples[i].position);
	}
}

void Algorithm_Roveri::_generatePeriodicBorders_generativeLabel()
{
	simulation.inputCloud.insideBorders = vector<int>(simulation.inputCloud.cloud.samples.size());

	// hack, the labels were set on calling this
	AlignedBox3f box;

	size_t i = 0;
	for (auto& ePtr : _exemplar)
	{
		simulation.inputCloud.insideBorders[i] = ePtr.generative;

		if (ePtr.generative)
		{
			box.extend(eigen(ePtr.node(_exemplar.graph).position));
		}
	}

	const Vector3f inset(simulation.neighSize * 1.2f, simulation.neighSize * 1.2f, simulation.is_2d ? 0.0f : simulation.neighSize * 1.2f);

	box.min() += inset;
	box.max() -= inset;

	simulation.periodic_border_bounds = box;
}


void Algorithm_Roveri::_generateBackgroundGridFromNeighSize()
{
	simulation.backgroundGrid.generateBackgroundPoints( 
		simulation.output_canvas_size, 
		simulation.neighSize, 
		simulation.is_2d, 
		simulation.neighSize 
	);
}

void Algorithm_Roveri::_generatedBackgroundGridBrushIndices_volume()
{
	simulation.backgroundGrid.brush_indices.clear();

	for(int i = 0; i<simulation.backgroundGrid.points_large.size(); i++)
	{
		// we could check if the points are inside of a limit here
		simulation.backgroundGrid.brush_indices.push_back( i );
	}

	vector<BackgroundPoint> new_bg_points_optimized;
	vector<int> new_brush_indices_optimized;

	for(int ii = 0; ii<simulation.backgroundGrid.brush_indices.size(); ii++)
	{
		int i = simulation.backgroundGrid.brush_indices[ii];
		new_bg_points_optimized.push_back( simulation.backgroundGrid.points_large[i] );
		new_brush_indices_optimized.push_back( ii );
	}

	//optimized version
	simulation.backgroundGrid.points.clear();
	simulation.backgroundGrid.points = new_bg_points_optimized;
	simulation.backgroundGrid.brush_indices.clear();
	simulation.backgroundGrid.brush_indices = new_brush_indices_optimized;

	_generateBackgroundGridMatchingPositions();
}

void Algorithm_Roveri::_generateBackgroundGridOrientations_volume()
{
	for(int ii = 0; ii<simulation.backgroundGrid.brush_indices.size(); ii++) 
	{
		int i = simulation.backgroundGrid.brush_indices[ii];

		Vector3f current_orientation( 0, 0, 1 );
		_alignInputAndMatchingPointWithBrushDirection3d( current_orientation, i );
		float current_scale = 1;
		_scaleInputAndMatchingPointWithBrushScale( current_scale, i );
	}
}

void Algorithm_Roveri::_alignInputAndMatchingPointWithBrushDirection3d( Vector3f dir_vector, int bg_point_index )
{
	Vector3f rot_axis = (dir_vector.normalized().cross( Vector3f( 1, 0, 0 ) )).normalized();
	float d = dir_vector.normalized().dot( Vector3f( 1, 0, 0 ) );
	float angle = -acos( d );
	
	if(!simulation.allow_rotation) 
	{
		angle = 0;
		dir_vector = Vector3f( 0, 1, 0 );
		rot_axis = Vector3f( 0, 1, 0 );
	}

	simulation.backgroundGrid.points[bg_point_index].orientation = dir_vector;

	//cout << d << " " << angle<< endl;
	Matrix3f m;
	m = AngleAxisf( angle, rot_axis );
	if(abs( angle ) < 0.15f) 
	{
		m << 1, 0, 0,
			0, 1, 0,
			0, 0, 1;
	}

	for(int i = 0; i< simulation.inputCloud.cloud.samples.size(); i++)
	{
		Vector3f shiftedPosition = simulation.inputCloud.cloud.samples[i].position;// - simulationManager.inputCloud.center_of_mass;
		Vector3f rotatedPosition = m * shiftedPosition;// + simulationManager.inputCloud.center_of_mass;
		simulation.backgroundGrid.points[bg_point_index].inputCloud_positions.push_back( rotatedPosition );
		Vector3f rotatedAttribute = m * simulation.inputCloud.cloud.samples[i].attributes[0];
		simulation.backgroundGrid.points[bg_point_index].inputCloud_attributes.push_back( rotatedAttribute );
	}

	simulation.backgroundGrid.points[bg_point_index].matching_position = m * simulation.backgroundGrid.points[bg_point_index].matching_position;
}

void Algorithm_Roveri::_generateBackgroundGridMatchingPositions()
{
	BackgroundGrid& backgroundGrid = simulation.backgroundGrid;
	{
		for(int i = 0; i < backgroundGrid.points.size(); i++) {
			bool found = false;
			do 
			{
				int rand_input_sample_index = rand() % simulation.inputCloud.cloud.samples.size();

				if(simulation.inputCloud.insideBorders[rand_input_sample_index] == 1 || simulation.periodic_border_type == SimulationManager::PeriodicBorderType::None)
				{
					backgroundGrid.points[i].matching_position = simulation.inputCloud.cloud.samples[rand_input_sample_index].position;

					found = true;
				}
			} while (!found);
		}
	}
}

void Algorithm_Roveri::_scaleInputAndMatchingPointWithBrushScale( float scale_value, int bg_point_index )
{
	if(!simulation.allow_scale) 
	{
		scale_value = 1;
	}

	simulation.backgroundGrid.points[bg_point_index].scale = scale_value;

	for(int i = 0; i< simulation.inputCloud.cloud.samples.size(); i++) 
	{
		simulation.backgroundGrid.points[bg_point_index].inputCloud_positions[i] *= scale_value;
	}

	simulation.backgroundGrid.points[bg_point_index].matching_position *= scale_value;
}

std::vector< Eigen::Vector3f > Algorithm_Roveri::_dummyAttributes()
{
	std::vector<Vector3f> attributes;

	attributes.emplace_back( Vector3f( 0.0f, 0.0f, 1.0f ) ); // normal

	return attributes;
}

AlgorithmResult Algorithm_Roveri::generate( std::vector<PositionFace>& positions, float radius /*= -1.0f*/, AABB limits /*= AABB() */ )
{
	AlgorithmResult result;

	if(!_didInit)
		_initialize();

	_context.domain.emptyRecycleBin();

	_seeding( positions, limits );

	beginRound();

	simulation.run();

	endRound(result);


	Cloud& outputCloud = simulation.outputCloud.cloud;

	int sizeDifference = outputCloud.samples.size() - _context.domain.size();

	AlgorithmResult generationResult;
	// add elements
	if(sizeDifference > 0)
	{
		size_t numToAdd = sizeDifference;
		size_t cloudStart = outputCloud.samples.size() - numToAdd;

		for(size_t i = 0; i < numToAdd; ++i)
		{
			FElement templateElement;
			templateElement.position = outputCloud.samples[cloudStart + i].position;

			auto elementHandle = _context.domain.insert(unreal(outputCloud.samples[cloudStart + i].position));

			generationResult.generated.push_back(elementHandle);
		}
	}
	// remove elements
	else if( sizeDifference < 0 )
	{
		size_t numToRemove = -sizeDifference;

		std::vector<FGraphNodeHandle> toRemove;

		size_t endIndex = _context.domain.size() - 1;

		for(size_t i = 0; i < numToRemove; ++i)
		{
			FElementObject& element = _context.domain.elementAt( endIndex - i );
			toRemove.push_back( element.nodeHandle() );

			generationResult.removed.push_back( element.nodeHandle() );
		}

		_context.domain.erase( toRemove );
	}
	 
	// copy data to elements
	{
		size_t i = 0;
		for(auto& sample : outputCloud.samples)
		{
			FElementObject& elementObject = _context.domain.elementAt( i );
			FGraphNode& elementNode = _context.graph().node(elementObject.nodeHandle());

			elementNode.position = unreal(sample.position);
			elementNode.scale = sample.hackScale;
			elementObject.type = _readType( sample );

			++i;
		}
	}


	// everything is modified
	for(auto& e : _context.domain)
	{
		result.modified.push_back( e.nodeHandle() );
	}

	result.append( generationResult );


	return result;
}

static const size_t typeBase = 1;

int16 Algorithm_Roveri::_readType( Sample& sample )
{
	auto& attributes = sample.attributes;

	assert( attributes.size() >= typeBase + 1 );

	float max = 0;
	int max_i = typeBase;
	int max_c = 0;

	for(int i = typeBase; i < attributes.size(); ++i)
	{
		auto& v = attributes[i];

		for(int c = 0; c < 3; ++c)
		{
			if(v( c ) > max)
			{
				max = v( c );
				max_i = i;
				max_c = c;
			}
		}
	}

	return (max_i - typeBase) * 3 + max_c;
}

void Algorithm_Roveri::_writeType( Sample& sample, int16 type )
{
	auto& attributes = sample.attributes;

	int start = type / 3 + typeBase;
	int offset = type % 3;
	int max = _maxType / 3 + typeBase;

	if(attributes.size() < max + 1)
		attributes.resize( max + 1 );

	for(int i = typeBase; i < attributes.size(); ++i)
	{
		attributes[i] = Eigen::Vector3f::Zero();
	}

	attributes[start]( offset ) = 1.0f;
}

void Algorithm_Roveri::clear()
{
	Algorithm::clear();

	simulation.clear();
}

void Algorithm_Roveri::loadExemplar()
{

}

void Algorithm_Roveri::endRound( AlgorithmResult& result )
{
	RoundRecord& record = roundSummary().back();

	record.end = std::clock();

	record.elementsGenerated = result.generated.size();
	record.elementsRemoved = result.removed.size();
	record.elementsModified = result.modified.size();

	record.elementsTotal = output().size();

	if(_context.domain.size() && simulation.show_average_energy)
	{
		if( simulation.averageEnergies.size() )
			record.totalEnergy_bruteForce = simulation.averageEnergies.back();
	}
}

void Algorithm_Roveri::_seeding( std::vector<PositionFace>& positions, AABB limits )
{
	using namespace Eigen;

	if(positions.size() == 0 || _exemplar.size() == 0)
		return;

	InputCloud& inputCloud = simulation.inputCloud;
	OutputCloud& outputCloud = simulation.outputCloud;

	Vector3f center = inputCloud.center_of_mass;
	float distSqrd = simulation.neighSize * simulation.neighSize;

	for(PositionFace& pair : positions)
	{
		auto nearest = _context.domain.nearest( pair.position, simulation.neighSize );

		if(!nearest.element)
		{
			size_t i = 0;
			for(size_t i = 0; i < inputCloud.cloud.samples.size(); i++)
			{
				Sample& si = inputCloud.cloud.samples[i];

				if((si.position - center).squaredNorm() < distSqrd)
				{
					Vector3f dir = si.position - center;
					Vector3f position = pair.position + dir;

					simulation.outputCloud.addSample( position, si.attributes );

					Sample& so = outputCloud.cloud.samples.back();

					so.hackScale = si.hackScale;
				}
			}
		}
	}







}

#undef _USE_MATH_DEFINES


