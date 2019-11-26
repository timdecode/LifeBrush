//
//  Algorithm_RegionGrowing.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2018-01-04.
//  Copyright (c) Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "Algorithm/Algorithm.h"
#include "Algorithm/Algorithm_RegionGrowing.h"

#include "dbscan.hpp"
#include "SpaceMapping.hpp"
#include "ShipEditorSimulation/MeshSimulation.h"

#include "Eigen/Sparse"
#include "Eigen/Geometry"

#include <vector>

AlgorithmResult Algorithm_RegionGrowing::_reassignSourceExamples()
{
	AlgorithmResult result;

	std::set<FGraphNodeHandle> elements;

	for(auto& brushPoint : _brushPoints.pts)
	{
		auto generationArea = _context.domain.nearestInRadius( brushPoint.position, brushPoint.radius );

		for(FGraphNodeHandle& e : generationArea)
		{
			FElementObject& eObject = _context.graph().component<FElementObject>(e);

			if(!enableVolumeSurfaceInteraction && !eObject.surfaceIndex.isOnSurface() && generationMode == EGenerationMode::SurfacePainting)
				continue;

			elements.insert( e );
		}
	}

	// reassign source examples to elements from the example selection
	for (FGraphNodeHandle elementHandle : elements)
	{
		ElementTuple element(elementHandle, graph());
		
		auto& sourceExample = _sourceExampleMap.sourceExample(elementHandle);


		auto nearest = _nearestSimilarElements(
			elementHandle,
			_context.domain,
			_defaultSelection._generativeElements,
			_exemplar,
			element.element.generationParameters,
			false, // force space filling
			false // don't filter the candidate by type (include it)
		);

		if (!nearest.size())
			continue;

		FGraphNodeHandle nearestHandle;

		// filter the nearest set to only contain elements in the selection
		for (auto& near : nearest)
		{
			if (_defaultSelection.selection.find(near.element) != _defaultSelection.selection.end())
			{
				nearestHandle = near.element;
				break;
			}
		}

		if (!nearestHandle)
			continue;

		sourceExample.element = nearestHandle;
		sourceExample.weight = _defaultSelection.weight;

		ElementTuple nearestElement(nearestHandle, _exemplar.graph);


		// hack reassign graph objects
		// for now, we'll just copy all of the other components to this guy
		element.node.removeComponents(graph());

		for (auto& componentType : nearestElement.node.components)
		{
			FGraphCopyContext::copyComponent(componentType, nearestElement.node, element.node, _exemplar.graph, graph());
		}

		// mark the mesh as dirty
		if (FGraphMesh * mesh = graph().componentPtr<FGraphMesh>(elementHandle))
		{
			mesh->markDirty();
		}

		element.node.scale = nearestElement.node.scale;

		result.modified.push_back(elementHandle);
	}

		//// remove a source example, if it's not in our selections
		//for(auto it = sourceExamples.begin(); it != sourceExamples.end(); )
		//{
		//	FGraphNodeHandle sourceElement = it->element;

		//	bool found = false;
		//	for(auto& selection : _exampleSelections)
		//	{
		//		if(selection->selection.find( sourceElement ) != selection->selection.end())
		//		{
		//			found = true;
		//			break;
		//		}
		//	}

		//	if(found)
		//		++it;
		//	else // it's not in our selections, nuke it
		//		it = sourceExamples.erase( it );
		//}

		//// add new source examples for the element
		//for(auto& selection : _exampleSelections)
		//{
		//	bool found = false;
		//	for(auto& sourceExample : sourceExamples)
		//	{
		//		if(selection->selection.find( sourceExample.element ) != selection->selection.end())
		//		{
		//			found = true;
		//			break;
		//		}
		//	}

		//	if(found)
		//		continue;

		//	// the sourceExamples does not contain any elements in the current selection
		//	// so, find one

		//	auto nearest = _nearestSimilarElements(
		//		element,
		//		_output,
		//		selection->_generativeElementsByType[element->type],
		//		_exemplar,
		//		element->generationParameters
		//	);

		//	if(!nearest.size())
		//		continue;

		//	FGraphNodeHandle nearestElement = nearest[0].element;

		//	sourceExamples.emplace_back(
		//		nearestElement,
		//		selection->weight, 
		//		selection );

		//	// hack reassign graph objects
		//	element->graphObjects = nearestElement->graphObjects;

		//	result.modified.push_back( element );
		//}



	return result;
}

AlgorithmResult Algorithm_RegionGrowing::_generate( std::vector<PositionFace>& startingPositions, float radius /*= -1.0f*/, AABB limits /*= AABB() */ )
{
	using namespace Eigen;
	using namespace std;

	AlgorithmResult result;

	_rebuildFreespace();

	if(enableReassignment)
		result = _reassignSourceExamples();

	// we only care about those elements around the origin
	// The local horizon are all those elements within the synthesize radius of a starting point or within
	// the neighbourhood of a brush point.
	std::set<FGraphNodeHandle> localHorizon;
	if(radius > 0.0f)
	{
		auto buildLocalHorizon = [&]( PositionRadiusFace origin )
		{
			PositionFace positionFace( origin.position, origin.surfaceIndex );

			auto generationArea = _context.domain.nearestInRadius( origin.position, origin.radius );

			for(auto e : generationArea)
			{
				if(_horizon.find( e ) == _horizon.end())
					continue;

				// don't add elements to the problem if we are generating on a surface and the element isn't on a surface
				if (enableVolumeSurfaceInteraction)
				{
					localHorizon.insert(e);
				}
				else
				{
					ElementTuple element(e, graph());

					if (element.element.surfaceIndex.isOnSurface() == origin.surfaceIndex.isOnSurface())
						localHorizon.insert(e);
				}
			}

			// do we need to add any seed points?
			// yes, if there are no elements nearby or if there are no points in our generation radius
			if(generationArea.size() != 0)
				return;

			if(!_occlusionTester->isVisible( positionFace, 1.0f ))
				return;


			_initializeOutput( positionFace, limits, result );

			_expandFreeSpace( result.generated );
		};

		for(auto origin : startingPositions)
		{
			PositionRadiusFace positionFace( origin.position, generationParameters.radius, origin.surfaceIndex );
			buildLocalHorizon( positionFace );
		}

		for(auto brushPoint : _brushPoints.pts)
		{
			buildLocalHorizon( brushPoint );
		}

		localHorizon = _removeFrozen_usingFreespacePoints( localHorizon, _context.domain, result.frozen );

		// update the global horizon
		for(auto e : result.frozen)
		{
			_horizon.erase( e );
		}
	}
	else
	{
		_horizon = _removeFrozen_usingFreespacePoints( _horizon, _context.domain, result.frozen );

		localHorizon = _horizon;
	}

	EigenPointCloud cloud;
	EigenDynamicPointCloudIndexAdaptor startingPositions_kdTree( 3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams( 10 ) );

	{
		std::vector<Eigen::Vector3f> asPoints;
		for(auto& pair : startingPositions)
			asPoints.push_back( pair.position );

		cloud.pts = asPoints;

		if(cloud.pts.size() > 0)
			startingPositions_kdTree.addPoints( 0, cloud.pts.size() - 1 );
	}

	// seeds!!
	auto seeds = _seeds( localHorizon );

	if(seeds.size() > 0)
	{
		// prune the seeds
		if(radius > 0)
		{
			seeds = _constrainSeedsToBounds( radius, seeds, startingPositions_kdTree );
		}

		// precompute our nearest elements
		// build the cache of nearest elements
		std::vector<std::vector<SimilarElement>> similarElementsCache;

		std::vector<std::future<std::vector<SimilarElement>>> futures;

		for(int i = 0; i < seeds.size(); ++i)
		{
			// emplace back, not push_back, so that we call the future move constructor
			futures.emplace_back( _context.threadPool->push( [&, i]( int threadID ) mutable {
				FGraphNodeHandle seedHandle = seeds[i];
				FElementObject& seedElement = graph().component<FElementObject>(seedHandle);

				std::vector<Algorithm::SimilarElement> nearest;


				auto filteredGenerativeExemplars = _fastFilterGenerativeExemplars( seedHandle, _freeSpace );

				nearest = _generation_nearestSimilarElements( seedHandle, _context.domain, filteredGenerativeExemplars, _exemplar, seedElement.generationParameters );



				//// nearest will have the exemplar samples ordered first
				//// followed by the unbiased nearest exemplars
				//for(auto& s : unbiasedNearest)
				//{
				//	nearest.push_back( s );
				//}

				//// resort
				//sort( nearestSimilarElements.begin(), nearestSimilarElements.end(), []( const SimilarElement& a, const SimilarElement& b ) -> bool
				//{
				//	return a.cost < b.cost;
				//} );

				return nearest;
			} ) );
		}

		for(auto& f : futures)
			similarElementsCache.push_back( f.get() );

		auto& exampleSelection = _highestWeightExampleSelection();
		auto& selectionElements = exampleSelection.selection;

		for(int i = 0; i < seeds.size(); ++i)
		{
			FGraphNodeHandle seedHandle = seeds[i];
			ElementTuple seed(seedHandle, graph());


			// Warning, this is very hacky, because it doesn't consider newly synthesized neighbours
			vector<SimilarElement>& similarElements = similarElementsCache[i];

			vector<FGraphNodeHandle> generatedElements;

			float lastCost = 0.0f;
			if(similarElements.size())
				lastCost = similarElements[0].cost;

			int assigned = 0;


			for(SimilarElement& similar : similarElements)
			{
				// second pass at the seed
				if(ignoreBadSuggestions)
				{
					// check if the suggestion is terrible
					if(similar.cost > similar.fullPairings.size() * ignoreBadSuggestionsDistanceFactor * typeCost)
					{
						continue;
					}
				}

				lastCost = similar.cost;

				FGraphNodeHandle exampleHandle = similar.element;
				ElementTuple example(exampleHandle, _exemplar.graph);

				for(auto& pair : similar.rightPartialPairings)
				{
					FGraphNodeHandle exHandle_ = pair.second;
					ElementTuple ex_(exHandle_, _exemplar.graph);

					Vector3f prediction = ex_.position() - example.position();
					float predictionNorm = prediction.norm();
					if(predictionNorm > seed.element.generationParameters.radius || predictionNorm > seed.element.generationInnerRadius)
						continue;

					// now add all the siblings of the entity
					// we should follow the element connections, but for now, we just set this to the element
					std::vector<FGraphNodeHandle> toAdd = { exHandle_ };

					std::vector<Eigen::Vector3f> newPositions;
					std::vector<FSurfaceIndex> newSurfaces;

					auto mapping = _mapping(seed, example);

					for(FGraphNodeHandle& ex_siblingHandle : toAdd)
					{
						ElementTuple ex_sibling(ex_siblingHandle, _exemplar.graph);

						auto toSurface = mapping->toSurface(ex_sibling.position() );
						if(!toSurface.hit)
							continue;

						Vector3f& newPosition = toSurface.position;

						if(selectionElements.find(ex_siblingHandle) == selectionElements.end())
							continue;

						if(_overlaps( newPosition, ex_sibling.element.radius, ex_sibling.element.generationParameters ))
							continue;

						if(!_inBoundary( newPosition, ex_sibling.element.radius ))
							continue;

						if(!limits.contains( newPosition, ex_sibling.element.radius ))
							continue;

						PositionFace facePair( newPosition, toSurface.surfaceIndex );
						if(!_occlusionTester->isVisible( facePair, ex_sibling.element.radius ))
							continue;

						newPositions.push_back( newPosition );
						newSurfaces.push_back( toSurface.surfaceIndex );

						// remap for the next element
						mapping = _mapping( newPosition, toSurface.surfaceIndex, ex_sibling.position() );
					}

					if(newPositions.size() != toAdd.size())
						continue;

					assigned++;

					// we should also copy the connections, but we don't
					// instantiate the selected elements to add
					for(int j = 0; j < newPositions.size(); ++j)
					{
						auto& newPosition = newPositions[j];
						FGraphNodeHandle ex_siblingHandle = toAdd[j];
						ElementTuple ex_sibling(ex_siblingHandle, exemplarGraph());

						auto newSurfaceIndex = newSurfaces[j];

						int elementIndex = _context.domain.size();

						FQuat newRotation = ex_sibling.node.orientation;

						if (generationMode == EGenerationMode::SurfaceProjection)
						{
							auto rotationAndNormal = _context.meshInterface->rotationAndNormalAtSurfacePoint(unreal(newPosition));

							FQuat quat = rotationAndNormal.first;
							FQuat rotation = quat * ex_sibling.node.orientation;

							newRotation = rotation;
						}
						else if (!forceRotation && (generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting))
						{
							auto rotationAndNormal = _context.meshInterface->rotationAndNormalAtIndex(newSurfaceIndex);

							FQuat quat = rotationAndNormal.first;
							FQuat rotation = quat * ex_sibling.node.orientation;

							newRotation = rotation;
						}
						else if (forceRotation)
							newRotation = forcedRotation;


						FGraphNodeHandle newHandle = _copyExemplarToOutput( ex_siblingHandle, exemplarGraph(), newPosition, exampleSelection, newRotation);

						ElementTuple newElement(newHandle, graph());

						newElement.element.surfaceIndex = newSurfaceIndex;



						result.generated.push_back( newHandle );
						localHorizon.insert( newHandle );
					}
				}

				if(assigned > 0)
					break;
			}

			// place the free-space points in low priority
			if(assigned == 0 && _canRemoveFromHorizon( seedHandle ))
			{
				_horizon.erase( seedHandle );
				localHorizon.erase( seedHandle );

				result.frozen.push_back( seedHandle );

				// remove all the free-space points around the removed element
				Eigen::Vector3f seedPosition = eigen(seed.node.position);
				auto freespaceIndices = _freeSpace.outputIndicesInRadius( seedPosition, seed.element.generationParameters.radius );

				for(auto& pair : freespaceIndices)
					_freeSpace.removeOutputPoint( pair.first );
			}
		}
	}




	if(!disableOptimization)
	{
		AlgorithmResult optimizationResult;

		OptimizationProblem optimizationProblem;

		if(useGlobalOptimization)
		{
			for(auto& e : _context.domain)
			{
				optimizationProblem.activeElements.insert( e.nodeHandle() );
				optimizationProblem.elements.insert(e.nodeHandle());
				optimizationProblem.frozenElements.insert(e.nodeHandle());
			}
		}
		else
		{

			optimizationProblem.activeElements = localHorizon;

			for(auto& e : _context.domain)
				optimizationProblem.elements.insert(e.nodeHandle());

			optimizationProblem.frozenElements = _nearbyFrozenElements( optimizationProblem.activeElements );
		}

		for(int i = 0; i < optimizationRounds; ++i)
		{
			AlgorithmResult optimizationResult;

			optimizationProblem = _localOptimization( optimizationProblem, optimizationResult );

			result.append( optimizationResult );
		}
	}

	_expandFreeSpace( result.generated );

	return result;
}

std::vector<FGraphNodeHandle> Algorithm_RegionGrowing::_constrainSeedsToBounds( float radius, std::vector<FGraphNodeHandle>& seeds, EigenDynamicPointCloudIndexAdaptor &startingPositions_kdTree )
{
	std::vector<FGraphNodeHandle> toKeep;

	float radius_sq = radius * radius;

	for(FGraphNodeHandle seedHandle : seeds)
	{
		ElementTuple seedElement(seedHandle, graph());

		Eigen::Vector3f p = seedElement.position();


		PositionFace facePair( p, seedElement.element.surfaceIndex );

		if(!_occlusionTester->isVisible( facePair, seedElement.element.radius ))
			continue;

		// keep the point if we near a starting point
		{
			size_t index = 0;
			float distance = 0.0f;

			startingPositions_kdTree.knnSearch( &p( 0 ), 1, &index, &distance );

			bool nearStarting = distance < radius_sq;

			if(nearStarting)
			{
				toKeep.push_back(seedHandle);
				continue;
			}
		}

		// keep the point if we are near a brush point
		{
			size_t index = 0;
			float distance = 0.0f;

			_brushIndex->knnSearch( &p( 0 ), 1, &index, &distance );

			bool nearBrush = distance < radius_sq;

			if(nearBrush)
				toKeep.push_back(seedHandle);
		}
	}

	return toKeep;
}

AlgorithmResult Algorithm_RegionGrowing::generate( std::vector<PositionFace>& startingPositions, float radius /*= -1.0f*/, AABB limits /*= AABB() */ )
{
	// initialize (don't count towards round time)
	if(!_didInit)
	{
		AlgorithmResult result;

		_initialize();

		beginRound();

		for(auto start : startingPositions)
			_initializeOutput( start, limits, result );

		_expandFreeSpace( result.generated );

		endRound( result );

		return result;
	}
	else
	{
		beginRound();

		_context.domain.emptyRecycleBin();

		AlgorithmResult r = _generate( startingPositions, radius, limits );

		endRound( r );

		return r;
	}
}

void Algorithm_RegionGrowing::_initialize()
{
	Algorithm::_initialize();

	_initFreespacePoints();
}