//
//  Algorithm.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-07-20.
//  Copyright (c) 2015 EpicGames. All rights reserved.
//

#include "LifeBrush.h"

#include <unordered_set>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <memory>
#include <ctime>


#include "Algorithm.h"
#include "dbscan.hpp"
#include "SpaceMapping.hpp"

#include "Eigen/Sparse"
#include "Eigen/Geometry"
#include "ShipEditorSimulation/GraphSnapshot.h"
#include "Simulation/Aggregates.h"
#include "Simulation/FlexElements.h"

using namespace Eigen;
using namespace std;

class NoOcclusion : public OcclusionBase
{
public:
	virtual bool isVisible( PositionFace point, float radius )
	{
		return true;
	}
};

// -----------------------------------------------------------------------------
/// Initialization
// -----------------------------------------------------------------------------

void Algorithm::_conditionally_initializeThreadPool()
{
	if(_context.threadPool == nullptr)
	{
		//unsigned int n = std::thread::hardware_concurrency();

		_context.threadPool = std::make_shared<ctpl::thread_pool>();
		_context.threadPool->resize( nThreads );
	}
}

void Algorithm::_initialize()
{
    _maxRadius = _computeMaxRadius();

	_initElementParameters();
    
	_initDefaultSelection();
    _initExemplarTypeHistogram();

	_initFreespacePoints();

	_initBrushIndex();
    
    if( _occlusionTester == nullptr )
    {
        _occlusionTester = std::make_shared<NoOcclusion>();
    }
    
    _didInit = true;
}

float Algorithm::_computeMaxRadius()
{
    float maxRadius = 0.0f;
    for( auto& e : _exemplar )
    {
        const float radius = e.radius;
        if( radius > maxRadius )
            maxRadius = radius;
    }
    
    return maxRadius;
}

void Algorithm::_initElementParameters()
{
	// use the element parameters as is
	if(perElementParameters)
		return;

	// don't use per-element parameters? then override all of them with the algorithm's parameters
	for(auto& e : _exemplar)
	{
		e.minAssignmentDistance = minAssignmentDistance;
		e.freespaceRadius = freespaceRadius;
		e.generationParameters = generationParameters;
		e.optimizationParameters = optimizationParameters;
		e.generationInnerRadius = generationInnerRadius;
	}
}

void Algorithm::_initDefaultSelection()
{
	_defaultSelection.clear();

	// create the default selection
	for(auto& e : _exemplar)
	{
		_defaultSelection.selection.insert( e.nodeHandle() );
	}

	_defaultSelection.init( exemplar(), *this );
}

void Algorithm::_initExemplarTypeHistogram()
{
    _exemplarTypeHistogram.clear();
    
    for( auto& ex : _exemplar )
        _exemplarTypeHistogram[ex.type]++;
}

void Algorithm::_initFreespacePoints()
{
	_freeSpace.clear();



	const auto& exampleElements = _defaultSelection.selection;

	for(auto& elementNodeHandle : exampleElements)
	{
		auto elementNode = exemplarNode(elementNodeHandle);

		std::vector<FGraphNodeHandle>& coherentExamples = _defaultSelection._kCoherentElements[elementNodeHandle];

		std::vector<Eigen::Vector3f> points;
		for(auto& c : coherentExamples)
		{
			auto cTuple = exemplarTuple(c);

			auto position = cTuple.node.position;

			auto neighbours = _exemplar.neighbours(c, cTuple.element.generationParameters.radius );

			Element::eachNode(neighbours, exemplarGraph(), [&](FGraphNode& c_) {
				points.push_back(eigen(c_.position - position + elementNode.position));
			});
		}

		Point_kdTree<float, Eigen::Vector3f> kdTree(points);

		using DBScanElements = rg::DBScan<float, Eigen::Vector3f>;

		std::vector< std::vector<DBScanElements::PointAndIndex> > clusters_out;
		std::vector<DBScanElements::PointAndIndex> noise_out;

		DBScanElements::dbscan(kdTree, freespaceRadius, 2, clusters_out, noise_out );

		std::vector<Eigen::Vector3f>& freePoints = _freeSpace.freespacePoints(elementNodeHandle);

		auto setupCluster = [&]( std::vector<DBScanElements::PointAndIndex>& dbCluster ) mutable
		{
			// find the midpoint and add use that as the free-space point
			Eigen::Vector3f centre = Eigen::Vector3f::Zero();

			for(auto& c : dbCluster)
				centre += c.point;

			centre /= dbCluster.size();

			freePoints.push_back( centre );
		};

		// add noise as their own little clusters
		for(auto& noise : noise_out)
		{
			std::vector<DBScanElements::PointAndIndex> oneNoise;
			oneNoise.push_back( noise );

			setupCluster( oneNoise );
		}

		for(auto & dbCluster : clusters_out)
			setupCluster( dbCluster );
	}
}

void Algorithm::reloadContext()
{
	std::vector<FGraphNodeHandle> allExamples = _exemplar.elements();

	std::vector< std::future< std::pair<FGraphNodeHandle, FGraphNodeHandle > > > futures;

	// we need to recalculate the sourceExample for every output domain element
	for (auto& elementObject : _context.domain)
	{
		auto& found = _defaultSelection._generativeElementsByType[elementObject.type];

		futures.emplace_back(_context.threadPool->push([this, elementObject, &allExamples, &found](int threadID) mutable {

			std::vector<FGraphNodeHandle>& generativeExemplars = found.size() ? found : allExamples;

			auto nearest = _nearestSimilarElements(
				elementObject.nodeHandle(),
				_context.domain,
				generativeExemplars,
				_exemplar,
				elementObject.generationParameters
			);

			FGraphNodeHandle similarExample = nearest.size() ? nearest[0].element : FGraphNodeHandle::null;


			return std::make_pair(elementObject.nodeHandle(), similarExample);
		}));
	}

	for (auto& future : futures)
	{
		auto pair = future.get();

		if (pair.second)
		{
			auto& sourceExample = _sourceExampleMap.sourceExample(pair.first);

			sourceExample.element = pair.second;
			sourceExample.weight = 1.0f;
		}
	}
}

void Algorithm::loadOutputElements( FGraphSnapshot& toLoad )
{
	clearPainting();
	clearOutput();

	toLoad.restore(graph());

	reloadContext();
}

Algorithm::ExampleSelection& Algorithm::_highestWeightExampleSelection()
{
	return _defaultSelection;
}

void Algorithm::_initializeOutput(PositionFace& position, AABB& limits, AlgorithmResult& result)
{
    if( _exemplar.size() == 0 )
        return;
    
    Vector3f samplingPoint = Vector3f::Zero();
    
	// where to sample the exemplar from

	Algorithm::ExampleSelection& exampleSelection = _highestWeightExampleSelection();
	auto& selectionElements = exampleSelection.selection;

	// sample from the centroid
	std::vector<FGraphNodeHandle> elementsToSample;
	if( selectionElements.size() )
		elementsToSample = std::vector<FGraphNodeHandle>( selectionElements.begin(), selectionElements.end() );
	else
		elementsToSample = _exemplar.elements();

	// find the centre of the exemplar
	Vector3f centre = Vector3f::Zero();

	for(auto& elementHandle : elementsToSample)
		centre += eigen(exemplarNode(elementHandle).position);

	centre /= float( elementsToSample.size() );

	samplingPoint = centre;

	// which is the closest element to the sampling point?
	FGraphNodeHandle closestHandle;
	{
		float minDistanceSqrd = std::numeric_limits<float>::max();
		for (FGraphNodeHandle& eHandle : elementsToSample)
		{
			FGraphNode& eNode = exemplarNode(eHandle);

			float dist = (eigen(eNode.position) - samplingPoint).squaredNorm();

			if(dist < minDistanceSqrd)
			{
				closestHandle = eHandle;
				minDistanceSqrd = dist;
			}
		}
	}

	if(!closestHandle)
		return;

	FElementObject& closestElement = _exemplar.graph.component<FElementObject>(closestHandle);

	auto neighbourhood = _exemplar.nearestInRadius( samplingPoint, closestElement.generationParameters.radius );

	std::set<FGraphNodeHandle> selectedNeighbourhood;

	for(auto& neighbour : neighbourhood)
	{
		FElementObject& neighbourElement = _exemplar.graph.component<FElementObject>(neighbour);

		if(selectionElements.size() == 0 || selectionElements.find(neighbour) != selectionElements.end() )
			selectedNeighbourhood.insert(neighbour);
	}
    
    auto mapping = position.surfaceIndex.isOnSurface() ? _mapping(position.position, position.surfaceIndex, samplingPoint) : _mapping(position.position, samplingPoint);
    
	for (auto& elementHandle : selectedNeighbourhood)
	{
		ElementTuple element(elementHandle, _exemplar.graph);


		auto p = mapping->toSurface(element.position());

		if (!p.hit)
			continue;

		if (!_inBoundary(p.position, element.element.radius))
			continue;

		if (!limits.contains(p.position, element.element.radius))
			continue;

		Eigen::Vector3f p_ = p.position;

		PositionFace facePair(p.position, p.surfaceIndex);

		if (selectionElements.find(elementHandle) == selectionElements.end())
			continue;

		if (!_occlusionTester->isVisible(facePair, element.element.radius))
			continue;

		if (_overlaps(p_, element.element.radius, element.element.generationParameters))
			continue;

		FQuat q_ = FQuat::Identity;

		if (!forceRotation && (generationMode == EGenerationMode::SurfaceProjection || generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting))
		{
			auto rotationAndNormal = _context.meshInterface->rotationAndNormalAtIndex(p.surfaceIndex);


			FQuat quat = rotationAndNormal.first;
			FQuat rotation = quat * element.node.orientation;

			//added.node.orientation = rotation;

			q_ = rotation;
		}
		else if (forceRotation)
		{
			//added.node.orientation = forcedRotation;

			q_ = forcedRotation;
		}
		else
			q_ = element.node.orientation;

		FGraphNodeHandle addedHandle = _copyExemplarToOutput(elementHandle, _exemplar.graph, p_, exampleSelection, q_);

		ElementTuple added(addedHandle, graph());

		added.element.surfaceIndex = p.surfaceIndex;



		result.generated.push_back(addedHandle);
	}

	_context.domain.rebalance();
}

FGraphNodeHandle Algorithm::_copyExemplarToOutput(FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph, const Eigen::Vector3f& position, Algorithm::ExampleSelection& selection, FQuat rotation)
{
	FGraph& theGraph = _context.graph();

	FGraphNodeHandle newNodeHandle = FGraphCopyContext::copyAggregate(sourceNodeHandle, sourceGraph, theGraph, unreal(position), rotation);

	FGraphNode& newNode = theGraph.node(newNodeHandle);

	if( !newNode.hasComponent<FElementObject>() )	
		newNode.addComponent<FElementObject>(theGraph);

	_horizon.insert(newNodeHandle);

	_context.domain.didInsert(newNodeHandle);

	auto& sourceExample = _sourceExampleMap.sourceExample(newNodeHandle);

	sourceExample.element = sourceNodeHandle;
	sourceExample.weight = 1.0f;

	return newNodeHandle;

	//FMLAggregateNO * aggregate = sourceGraph.componentPtr<FMLAggregateNO>(sourceNodeHandle);

	//if (aggregate)
	//{
	//	TArray<FGraphNodeHandle> aggregateNodes;
	//	TArray<FGraphEdgeHandle> aggregateEdges;

	//	aggregate->edgesAndNodesInAggregate(sourceGraph, aggregateNodes, aggregateEdges);

	//	auto edgesToCopy = sourceGraph.edgesBetweenNodes(aggregateNodes);

	//	FGraphNode& sourceNode = sourceGraph.node(sourceNodeHandle);

	//	const FVector translation = unreal(position) - sourceNode.position;
	//	const FQuat rotation_fixed = rotation; // align the vector (otherwise unreal crashes on simd shit)
	//	const FQuat deltaRotation = rotation_fixed * sourceNode.orientation.Inverse();

	//	// Yea this is ugly
	//	const FTransform transform = FTransform(-sourceNode.position) * FTransform(deltaRotation) * FTransform(sourceNode.position) * FTransform(translation);


	//	auto result = FGraphCopyContext::copySubgraph(sourceGraph, _context.graph(), aggregateNodes, edgesToCopy, transform);

	//	return result.duplicatedNodes[0];
	//}
	//else
	//{
	//	FElementObject& sourceElement = sourceGraph.component<FElementObject>(sourceNodeHandle);

	//	auto newNodeHandle = _context.domain.insert(sourceElement, sourceGraph, unreal(position));

	//	_horizon.insert(newNodeHandle);

	//	_context.graph().node(newNodeHandle).orientation = rotation;

	//	auto& sourceExample = _sourceExampleMap.sourceExample(newNodeHandle);

	//	sourceExample.element = sourceNodeHandle;
	//	sourceExample.weight = 1.0f;

	//	return newNodeHandle;
	//}
}

void Algorithm::updateParticleBVH()
{
	FGraph& theGraph = _context.graph();

	auto& particles = theGraph.componentStorage<FFlexParticleObject>();

	const float hackRadius = 0.6;

	// update or insert into the BVH
	for (auto& particle : particles)
	{
		auto particleIndex = particle.nodeHandle().index;

		if (!particle.isValid())
		{
			if (_particleBVH.containsParticle(particleIndex))
				_particleBVH.removeParticle(particleIndex);
		}

		if (!particle.isValid()) continue;

		FVector position = theGraph.node(particle.nodeHandle()).position;

		if (_particleBVH.containsParticle(particleIndex))
			_particleBVH.updateParticle(particleIndex, position, hackRadius);
		else
			_particleBVH.insertParticle(particleIndex, position, hackRadius);
	}
}

void Algorithm::clearAlongPath(std::vector<PositionRadiusFace>& path)
{
	std::vector<FGraphNodeHandle> toRemove;

	for (auto& pathPoint : path)
	{
		Eigen::Vector3f& position = pathPoint.position;

		// remove elements
		auto elements = _context.domain.nearestInRadius(position, pathPoint.radius);

		for (FGraphNodeHandle& e : elements)
			toRemove.push_back(e);

		// remove freespace output points
		auto indices = _freeSpace.outputIndicesInRadius(position, pathPoint.radius);

		for (auto& pair : indices)
		{
			size_t i = pair.first;

			_freeSpace.removeOutputPoint(i);
		}
	}

	removeOutputElements(toRemove);
}

void Algorithm::removeOutputElements( std::vector<FGraphNodeHandle>& elements )
{
    // clear the horizon
    for(auto& e : elements)
    {
		_horizon.erase( e );
    }

    
    // clear the source examples
    for(auto& e : elements)
    {
		_sourceExampleMap.eraseElement( e );
    }

    
    // clear out the freespace
	// we should clear all output points attached to the elements, but we don't have that information, yet

	/*_freeSpace.removeOutputElements( elements );

    _rebuildFreespace();*/


    // clear the output
    _context.domain.erase( elements );
}

void Algorithm::addToHorizon( std::vector<FGraphNodeHandle>& elements )
{
    for(auto& e : elements)
        _horizon.insert( e );
}


// -----------------------------------------------------------------------------
/// LSA
// -----------------------------------------------------------------------------

void Algorithm::_sampleNeighbourhood(FGraphNodeHandle elementHandle, Domain& domain, Eigen::VectorXf& vector_out)
{
	ElementTuple element(elementHandle, domain.graph);

    auto neighbourhood = domain.neighbours(elementHandle, element.element.generationParameters.radius);
    neighbourhood.push_back(elementHandle);
    
    int n = 8;
    float cellSize = element.element.generationParameters.radius / float(n);
    
    // eigen is column major (column by column)
    vector_out = Eigen::VectorXf(n * n * n);
    typedef Eigen::Triplet<float> Triplet;
    
    for( auto& e : neighbourhood )
    {
        Vector3f p = eigen(e(domain.graph).position - element.node.position);
        p /= cellSize;
        
        vector_out(int(p.z()) * n * n + int(p.y()) * n + int(p.x())) += 1.0f;
    }
    
    
    // or
    
    for( int xi = 0; xi < n; ++xi )
    {
        for( int yi = 0; yi < n; ++yi )
        {
            for( int zi = 0; zi < n; ++zi )
            {
                
            }
        }
    }
    
}

auto Algorithm::_discretize(FGraphNodeHandle elementHandle, Domain& domain, const int n) -> Eigen::VectorXf
{
    // the bounding box extends from (p_e - r) to (p_e + r).
    const float r = generationParameters.radius;
    const float delta = (r * 2.0f) / float(n - 1);
    const float invDelta = 1.0f / delta;
    
    const Vector3f rvec = Vector3f(r, r, r);
    const Vector3f halfExtents = Vector3f(delta / 2.0f, delta / 2.0f, delta / 2.0f);
    
    const Vector3i toIndex = Vector3i(0, n, n * n);
    
    auto neighbourhood = domain.neighbours(elementHandle, generationParameters.radius);
    neighbourhood.push_back(elementHandle);

	auto elementPosition = eigen(domain.graph.node(elementHandle).position);
    
    VectorXf volume(n * n * n);
    for( int i = 0; i < n * n * n; ++i )
        volume(i) = 0.0f;
        
        for( auto e : neighbourhood )
        {
			FGraphNode& eNode = domain.graph.node(e);

            // shift to a bounding box from (0,0,0) to (2r, 2r, 2r)
            Vector3f p = eigen(eNode.position) - elementPosition + rvec;
            
            p = (p + halfExtents) * invDelta;
            
            Eigen::Vector3i ivec = p.cast<int>();
            
            int i = ivec.dot(toIndex);
            
            volume[i] += 1.0f;
        }
    
    return volume;
}

void Algorithm::_initLSA()
{
    const int n = 8;
    
    const int n3 = n * n * n;
    
    const int m = _exemplar.size();
    
    Eigen::MatrixXf A(n3, m);
    
    int j = 0;
    
    for( auto& e : _exemplar )
    {
        auto vec = _discretize(e.nodeHandle(), _exemplar, n);
        
        A.col(j) = vec;
        
        j++;
    }
    
    Eigen::JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
    
    const int k = 5;
    
    VectorXf singularValues = svd.singularValues();
    
    MatrixXf U_k = svd.matrixU().leftCols<k>();
    VectorXf singularValues_k = singularValues.head<k>();
    MatrixXf V_k = svd.matrixV().leftCols<k>();
    
    MatrixXf A_k = U_k * singularValues_k.asDiagonal() * V_k.transpose();
    
    VectorXf q = A.col(0);
    VectorXf answer = A.transpose() * q;
    
    VectorXf singularValuesInverse_k(k);
    for( int i = 0; i < k; ++i )
        singularValuesInverse_k(i) = singularValues_k(i) != 0.0f ? 1.0f / singularValues_k(i) : 0.0f;
    
    
    queryTransformMatrix = U_k * singularValuesInverse_k.asDiagonal();
    lowDimensionalNeighbourhoodMatrix = V_k;
    
    VectorXf q_k = q.transpose() * queryTransformMatrix;
    VectorXf answer_k = lowDimensionalNeighbourhoodMatrix * q_k;
    
    
    
    for( int i = 0; i < answer.size(); ++i )
        UE_LOG(LogTemp, Warning, TEXT("%d %f"), i, answer(i));
    
    for( int i = 0; i < answer_k.size(); ++i )
        UE_LOG(LogTemp, Warning, TEXT("%d %f"), i, answer_k(i));
    
    std::stringstream out;
    
    out << "U_k\n" << U_k << "\nV_k---\n" << V_k << "\nS_k---\n" << singularValues_k << "\nA_k---\n" << A_k;
    
    FString fString(out.str().c_str());
    
    UE_LOG(LogTemp, Warning, TEXT("%s"), *fString);
    
    return;
}

void Algorithm::setOcclusionTester( std::shared_ptr<OcclusionBase> const& occlusionTester )
{
	if( occlusionTester == nullptr )
		_occlusionTester = std::make_shared<NoOcclusion>();
	else
		_occlusionTester = occlusionTester;
}

void Algorithm::init( std::shared_ptr<ctpl::thread_pool> threadPool )
{
	if(threadPool)
		setThreadPool( threadPool );
	else
		_conditionally_initializeThreadPool();
}


void Algorithm::loadExemplar()
{

}

void Algorithm::setMeshInterface( std::shared_ptr<tcodsMeshInterface> meshInterface )
{
	_context.meshInterface = meshInterface;
}

void Algorithm::setToMeshTransform( const Eigen::Affine3f& transform )
{
	_toMeshTransform = transform;
}

// -----------------------------------------------------------------------------
/// Generation
// -----------------------------------------------------------------------------

AlgorithmResult Algorithm::generate( std::vector<PositionFace>& startingPositions, float radius, AABB limits)
{
	return AlgorithmResult();
}

void Algorithm::_expandFreeSpace(std::vector<FGraphNodeHandle>& generated)
{
    // update the freespace index
	Element::eachTuple(generated, graph(), [&](ElementTuple element)
	{
		Eigen::Vector3f position = element.position();

		auto& sourceExample = _sourceExampleMap.sourceExample(element.handle());

		if (!sourceExample.element)
			return;

		auto sourceElement = exemplarTuple(sourceExample.element);

		auto mapping = _mapping(element, sourceElement);

		float searchRadiusSqrd = freespaceRadius * freespaceRadius;

		auto points = _freeSpace.freespacePoints(sourceElement.handle());
		for (Eigen::Vector3f& point : points)
		{
			// we need to project to the plane and offset the positions
			auto onSurface = mapping->toSurface(point);

			if (!onSurface.hit)
				continue;

			Eigen::Vector3f& onSurfacePoint = onSurface.position;

			// check for overlap
			auto nearest = _context.domain.nearest(onSurfacePoint, freespaceRadius);

			if (!nearest.element)
			{
				_freeSpace.addOutputPoint(onSurfacePoint);
			}
		}

		// remove any free-space points that overlap with the generated element
		auto freespaceIndices = _freeSpace.outputIndicesInRadius(position, freespaceRadius);

		for (auto& indicesPair : freespaceIndices)
			_freeSpace.removeOutputPoint(indicesPair.first);
	});
}

bool Algorithm::_overlaps(const Eigen::Vector3f& position, const float radius, FNeighbourhoodParameters& parameters, const FGraphNodeHandle elementToIgnore /*= nullptr*/)
{
    // we have to use _maxRadius for the query, since our radius and neighbours radius may be different
    //    float maxAllowedRadius = (radius + _maxRadius) / relaxation;
    float searchRadius = parameters.radius;
    
    auto neighbours = _context.domain.nearestInRadius(position, searchRadius);
    
    for( FGraphNodeHandle& neighbourHandle : neighbours )
    {
        if( elementToIgnore == neighbourHandle)
            continue;
        
		ElementTuple neighbour = outputTuple(neighbourHandle);

        float allowedRadius = (radius + neighbour.element.radius) * relaxation;
        
        float distance = (neighbour.position() - position).norm();
        if( distance < allowedRadius  )
            return true;
    }

	// check for particle overlap
	{
		FVector unrealPosition = unreal(position);

		const float hackParticleRadius = 0.6;

		unrealAABB::AABB query(unrealPosition, hackParticleRadius);

		bool hasHit = false;

		// check if the rule application overlaps with any occluders
		// - we don't overlap with the target of the rule
		// - we might need to exclude all of the handles in the rule
		_particleBVH.query(query, [&](unsigned int particleIndex) {
			hasHit = true;

			// don't keep querying, we are done
			return false;
		});

		if (hasHit)
			return true;
	}

    
    return false;
}

bool Algorithm::_inBoundary(const Vector3f& position, const float radius)
{
    // check boundary conditions
    if( generationMode == EGenerationMode::SpaceFilling )
    {
		if(seedsIgnoreMeshInBoundary)
			return true;

        auto nearest = _context.meshInterface->nearestPointOnMesh(unreal(position));
        
        float radiusSquared = radius * radius;
        Eigen::Vector3f nearestPoint = eigen(nearest.point);
        if( (nearestPoint - position).squaredNorm() < radiusSquared )
            return false;
        
        // look up the face and see if its normal is pointing away from us
        auto face = _context.meshInterface->mesh(nearest.surfaceIndex.sectionIndex).face(nearest.surfaceIndex.faceIndex);
        auto normal = eigen( face->normal() );

        Eigen::Vector3f direction = eigen(nearest.point) - position;
        
        direction.normalize();
        
        if( flipSurfaceNormals )
            direction = direction * -1.0f;
        
        if( normal.dot(direction) >= 0.0f )
            return false;
    }
    else if( generationMode == EGenerationMode::SurfacePainting || generationMode == EGenerationMode::SpacePainting )
    {
		if(_brushPoints.pts.size() == 0)
			return false;
       
		size_t foundIndex;
		float dSqrd;
		
		_brushIndex->knnSearch( &position.x(), 1, &foundIndex, &dSqrd );

		auto brushPoint = _brushPoints.pts[foundIndex];
		float brushSizeSqrd = brushPoint.radius * brushPoint.radius;

		if(dSqrd >= brushSizeSqrd)
			return false;
    }

    return true;
}

/**
 Marks elements as frozen if they do not have any right-partial-pairings.
 */
std::set<FGraphNodeHandle> Algorithm::_removeFrozen(const std::set<FGraphNodeHandle>& elements, Domain& domain, std::vector<FGraphNodeHandle>& frozen_out)
{
    std::set<FGraphNodeHandle> result;
    
    for( FGraphNodeHandle eHandle : elements )
    {
		ElementTuple e(eHandle, domain.graph);
        auto eNeighbours = domain.neighbours(eHandle, e.element.generationParameters.radius, e.element.generationParameters.kNearest);
        
        unsigned int rightPairingsCount = 0;
        
		auto& sourceExample = _sourceExampleMap.sourceExample( eHandle );

		bool foundRightPairings = false;


		auto& coherentNeighbours = sourceExample.kCoherentNeighbours( _defaultSelection );

		auto mapping = _mapping( e );

		for(FGraphNodeHandle candidateHandle : coherentNeighbours)
		{
			ElementTuple candidate(candidateHandle, exemplarGraph());
			if(e.element.type != candidate.element.type)
				continue;

			mapping->setExemplarOrigin( candidate.position() );

			auto candidateNeighbours = _exemplar.neighbours( candidateHandle, candidate.element.generationParameters.radius, candidate.element.generationParameters.kNearest );

			SimilarElement similar = distance( e.node, eNeighbours, domain.graph, candidate.node, candidateNeighbours, _exemplar.graph, mapping.get() );

			if(similar.rightPartialPairings.size())
			{
				result.insert( eHandle );
				foundRightPairings = true;
				break;
			}
		}
		
		if(!foundRightPairings)
			frozen_out.push_back( eHandle );
    }
    
    return result;
}

std::set<FGraphNodeHandle> Algorithm::_removeFrozen_usingFreespacePoints(const std::set<FGraphNodeHandle>& elements, Domain& domain, std::vector<FGraphNodeHandle>& frozen_out)
{
    std::set<FGraphNodeHandle> result;
    
    for( FGraphNodeHandle eHandle : elements )
    {
		ElementTuple e(eHandle, domain.graph);
		Eigen::Vector3f ePosition = e.position();

		bool hasNearest = _freeSpace.nearestInRadius(ePosition, e.element.generationParameters.radius) >= 0;

        if(!hasNearest && _canRemoveFromHorizon(eHandle) )
            frozen_out.push_back(eHandle);
        else
            result.insert(eHandle);
    }
    
    return result;
}

std::vector<FGraphNodeHandle> Algorithm::_seeds(const std::set<FGraphNodeHandle>& elements_in)
{
    vector<FGraphNodeHandle> result;
    
    if( elements_in.size() == 0 )
        return result;
    
    std::vector<FGraphNodeHandle> elements;


	std::copy( elements_in.begin(), elements_in.end(), std::back_inserter( elements ) );
    
    // build a local kd-tree for the elements
    EigenPointCloud cloud;

	EigenStaticPointCloudAdaptor index( 3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams( 10 ) );

    
    for(FGraphNodeHandle elementHandle : elements )
    {
        const auto p = elementHandle(graph()).position;
        cloud.pts.emplace_back(eigen(p));
    }

	index.buildIndex();
    
	FGraphNodeHandle elementHandle = FGraphNodeHandle::null;

    unordered_set<FGraphNodeHandle> available;
	for (FGraphNodeHandle e : elements)
	{
		ElementTuple element(e, graph());

		if( element.element.generationParameters.kNearest == 0 )
			continue;

		if (!elementHandle)
			elementHandle = e;

		available.insert(e);
	}
    
	if(elementHandle)
	{
		int k = graph().component<FElementObject>(elementHandle).generationParameters.kNearest;

		std::vector<size_t> indices(k);
		std::vector<float> distances(k);

		while (available.size() > 0)
		{
			ElementTuple element(elementHandle, graph());

			result.push_back(elementHandle);
			available.erase(elementHandle);

			Vector3f p = element.position();

			k = element.element.generationParameters.kNearest;

			// the std spec says that capacity will not be reduced, but the size will (we want this behavior for performance avoid too many allocations
			indices.resize(k);
			distances.resize(k);

			nanoflann::KNNResultSet<float, size_t> resultSet(k);
			resultSet.init(&indices[0], &distances[0]);
			index.findNeighbors(resultSet, &p(0), nanoflann::SearchParams());

			FGraphNodeHandle next;

			// remove all the elements in the interval between this and the seedSeparation'th element
			int neighbourIndex = 1; // skip e
			for (int i = 0; i < seedSeparation && neighbourIndex < resultSet.size(); ++i, ++neighbourIndex)
			{
				assert(neighbourIndex < k);

				int indicesIndexOfNeighbour = indices[neighbourIndex];

				FGraphNodeHandle neighbour = elements.at(indicesIndexOfNeighbour);

				if (available.find(neighbour) != available.end())
				{
					available.erase(neighbour);
					next = neighbour;
				}
			}

			if (next)
				elementHandle = next;
			else if (available.size() > 0)
				elementHandle = *available.begin();
		}



	}

    if( result.size() )
    {
		Eigen::Vector3f p0 = eigen(result[0](graph()).position);
        _orderNearest(result, graph(), p0);
        return result;
    }
    else
        return result;
}

void Algorithm::_orderNearest(std::vector<FGraphNodeHandle>& elements, FGraph& graph, const Vector3f referencePoint)
{
    sort(elements.begin(), elements.end(), [referencePoint,&graph](FGraphNodeHandle a, FGraphNodeHandle b)
    {
		Eigen::Vector3f pa = eigen(a(graph).position);
		Eigen::Vector3f pb = eigen(b(graph).position);
             return (pa - referencePoint).squaredNorm() < (pb - referencePoint).squaredNorm();
    });
}

// -----------------------------------------------------------------------------
/// Distance
// -----------------------------------------------------------------------------

Algorithm::SimilarElement Algorithm::distance_densityBased
( 
	FGraphNode& element, std::vector<FGraphNodeHandle>& elementNeighbours, FGraph& elementGraph,
	FGraphNode& example, std::vector<FGraphNodeHandle>& exampleNeighbours, FGraph& exampleGraph,
	rg::Mapping* mapping 
)
{
	int elementsSize = elementNeighbours.size();
	int examplesSize = exampleNeighbours.size();

	vector<Vector3f> elementPredictions( elementsSize );
	vector<Vector3f> examplePredictions( examplesSize );

	const Vector3f examplePosition = eigen(example.position);
	const Vector3f elementPosition = eigen(element.position);

	bool workInExemplarSpace = mapping->compareInExemplarSpace();

	if(workInExemplarSpace)
	{
		for (int i = 0; i < elementsSize; ++i)
		{
			FGraphNode& node = elementNeighbours[i](elementGraph);
			elementPredictions[i] = mapping->toExemplarVector(eigen(node.position) - elementPosition);
		}

		for (int i = 0; i < examplesSize; ++i)
		{
			FGraphNode& node = exampleNeighbours[i](exampleGraph);
			examplePredictions[i] = eigen(node.position) - examplePosition;
		}
	}
	else
	{
		for (int i = 0; i < elementsSize; ++i)
		{
			FGraphNode& node = elementNeighbours[i](elementGraph);
			elementPredictions[i] = eigen(node.position);
		}

		for(int i = 0; i < examplesSize; ++i)
		{
			FGraphNode& node = exampleNeighbours[i](exampleGraph);

			auto surface = mapping->toSurface( eigen(node.position) );
			examplePredictions[i] = surface.position;
		}
	}

	vector<pair<int, int>> indexPairs;
	float cost = 0.0f;

	for(int i = 0; i < elementsSize; ++i)
	{
		const Vector3f& p = elementPredictions[i];

		for(int j = 0; j < examplesSize; ++j)
		{
			const Vector3f& p_ = examplePredictions[j];

			// https://en.wikipedia.org/wiki/Multivariate_normal_distribution ?
		}
	}




	SimilarElement result;
	result.element = example.handle();
	result.cost = cost;

	// convert pairs of indices to element pairings
	// --------------------------------------------
	auto& fullPairings = result.fullPairings;
	auto& leftPartialPairings = result.leftPartialPairings;
	auto& rightPartialPairings = result.rightPartialPairings;


	vector<bool> assignedElements( elementsSize );
	vector<bool> assignedExamples( examplesSize );

	float minDistanceSquared = minAssignmentDistance * minAssignmentDistance;

	for(auto pair : indexPairs)
	{
		if(pair.first < 0 || pair.second < 0)
			break;

		FGraphNodeHandle e = elementNeighbours[pair.first];
		FGraphNodeHandle ex = exampleNeighbours[pair.second];

		Vector3f p_e = elementPredictions[pair.first];
		Vector3f p_ex = examplePredictions[pair.second];

		if((p_e - p_ex).squaredNorm() < minDistanceSquared)
		{
			fullPairings.emplace_back( e, ex );

			assignedElements[pair.first] = true;
			assignedExamples[pair.second] = true;
		}
	}

	for(int i = 0; i < elementsSize; ++i)
	{
		if(assignedElements[i])
			continue;

		leftPartialPairings.emplace_back( elementNeighbours[i], FGraphNodeHandle::null);
	}

	for(int i = 0; i < examplesSize; ++i)
	{
		if(assignedExamples[i])
			continue;

		rightPartialPairings.emplace_back(FGraphNodeHandle::null, exampleNeighbours[i] );
	}



	return result;
}

void Algorithm::beginRound()
{
	RoundRecord& record = roundSummary().emplace();

	record.start = std::clock();
}

void Algorithm::endRound( AlgorithmResult& result )
{
	RoundRecord& record = roundSummary().back();

	record.end = std::clock();

	record.elementsGenerated = result.generated.size();
	record.elementsRemoved = result.removed.size();
	record.elementsModified = result.modified.size();

	record.elementsTotal = output().size();

	if(enableRoundSummaries)
	{
		if(_context.domain.size() && calculate_kCoherenceEnergy)
		{
			record.totalEnergy_kCoherence = totalEnergy_kCoherence( optimizationParameters );
			record.averageEnergy_kCoherence = record.totalEnergy_kCoherence / float( _context.domain.size() );
		}

		if(_context.domain.size() && calculate_bruteForceEnergy)
		{
			record.totalEnergy_bruteForce = totalEnergy_bruteForce( optimizationParameters );
			record.averageEnergy_bruteForce = record.totalEnergy_bruteForce / float( _context.domain.size() );
		}
	}
}

auto Algorithm::totalEnergy_kCoherence( FNeighbourhoodParameters parameters ) -> double
{
	std::vector< std::future<float> > futures;

	for(auto& element : _context.domain)
	{
		futures.emplace_back( _context.threadPool->push( [this, parameters, &element]( int threadID ) mutable {
			float cost = 0.0f;
			
			auto& sourceExample = _sourceExampleMap.sourceExample( element.nodeHandle() );


			auto& coherentNeighbours = sourceExample.kCoherentNeighbours( _defaultSelection );

			std::vector<FGraphNodeHandle> generative;
			for(FGraphNodeHandle e : coherentNeighbours)
			{
				FElementObject& eObject = exemplarGraph().component<FElementObject>(e);

				if(eObject.generative)
					generative.push_back( e );
			}

			if(useCinpactEnergy)
			{
				auto nearestVec = _nearestElements_discreteCinpactSimilarity( element.nodeHandle(), _context.domain, generative, _exemplar );

				if(nearestVec.size())
					cost += nearestVec[0].cost;
				else
					cost += 0.0f;
			}
			else
			{
				auto nearestVec = _nearestSimilarElements(element.nodeHandle(), _context.domain, generative, _exemplar, parameters );

				if(nearestVec.size())
					cost += nearestVec[0].cost;
				else
					cost += 0.0f;
			}
			
			return cost;
		} ) );
	}

	double energy = 0.0;

	for(auto& f : futures)
		energy += f.get();

	return energy;
}

auto Algorithm::totalEnergy_bruteForce( FNeighbourhoodParameters parameters ) -> double
{
	std::vector< std::future<float> > futures;

	for(auto& element : _context.domain)
	{
		futures.emplace_back( _context.threadPool->push( [this, parameters, &element]( int threadID ) mutable {
			if(useCinpactEnergy)
			{
				auto nearestVec = _nearestElements_discreteCinpactSimilarity( element.nodeHandle(), _context.domain, _defaultSelection._generativeElements, _exemplar );

				if(nearestVec.size())
					return nearestVec[0].cost;
				else
					return 0.0f;
			}
			else
			{
				auto nearestVec  = _nearestSimilarElements( element.nodeHandle(), _context.domain, _defaultSelection._generativeElements, _exemplar, parameters );

				if(nearestVec.size())
					return nearestVec[0].cost;
				else
					return 0.0f;
			}
		} ) );
	}

	double energy = 0.0;

	for(auto& f : futures)
		energy += f.get();

	return energy;
}

auto Algorithm::sourceExamples(FGraphNodeHandle element) -> std::vector<FGraphNodeHandle>
{
	std::vector<FGraphNodeHandle> result;

	auto& sourceExample = _sourceExampleMap.sourceExample(element);

	if(sourceExample.element)
		result.push_back(sourceExample.element);

	return result;
}

Algorithm::SimilarElement Algorithm::distance(FGraphNode& elementNode, std::vector<FGraphNodeHandle>& elementNeighbours, FGraph& elementGraph,
                                              FGraphNode& exampleNode, std::vector<FGraphNodeHandle>& exampleNeighbours, FGraph& exampleGraph,
                                              rg::Mapping* mapping)
{
	ElementTuple element(elementNode, elementGraph);
	ElementTuple example(exampleNode, exampleGraph);

    int elementsSize = elementNeighbours.size();
    int examplesSize = exampleNeighbours.size();
    
    vector<Vector3f> elementPredictions(elementsSize);
    vector<Vector3f> examplePredictions(examplesSize);
    
    const Vector3f examplePosition = example.position();
    const Vector3f elementPosition = element.position();
    
	bool workInExemplarSpace = mapping->compareInExemplarSpace();

	if(workInExemplarSpace)
	{
		for (int i = 0; i < elementsSize; ++i)
		{
			Eigen::Vector3f p_i = eigen(elementNeighbours[i](elementGraph).position);
			elementPredictions[i] = mapping->toExemplarVector(p_i - elementPosition);
		}

		for (int i = 0; i < examplesSize; ++i)
		{
			Eigen::Vector3f p_i = eigen(exampleNeighbours[i](exampleGraph).position);
			examplePredictions[i] = p_i - examplePosition;
		}
	}
	else
	{
		for (int i = 0; i < elementsSize; ++i)
		{
			Eigen::Vector3f p_i = eigen(elementNeighbours[i](elementGraph).position);
			elementPredictions[i] = p_i;
		}

		for(int i = 0; i < examplesSize; ++i)
		{
			Eigen::Vector3f p_i = eigen(exampleNeighbours[i](exampleGraph).position);

			auto surface = mapping->toSurface(p_i);
			examplePredictions[i] = surface.position;
		}
	}
    
    vector<pair<int, int>> indexPairs;
    float cost = 0.0f;

	if(element.element.type != example.element.type)
		cost += typeCost;

    // this is where we would normally do ICP
    cost = _closestPoint(elementPredictions, elementNeighbours, elementGraph, examplePredictions, exampleNeighbours, exampleGraph, indexPairs);
    
    SimilarElement result;
    result.element = example.handle();
    result.cost = cost;
    
    // convert pairs of indices to element pairings
    // --------------------------------------------
    auto& fullPairings = result.fullPairings;
    auto& leftPartialPairings = result.leftPartialPairings;
    auto& rightPartialPairings = result.rightPartialPairings;
    
    
    vector<bool> assignedElements(elementsSize);
    vector<bool> assignedExamples(examplesSize);
    
    float minDistanceSquared = minAssignmentDistance * minAssignmentDistance;
    
    for( auto pair : indexPairs )
    {
        if( pair.first < 0 || pair.second < 0 )
            break;
        
        FGraphNodeHandle e = elementNeighbours[pair.first];
        FGraphNodeHandle ex = exampleNeighbours[pair.second];
        
        Vector3f p_e = elementPredictions[pair.first];
        Vector3f p_ex = examplePredictions[pair.second];
        
        if( (p_e - p_ex).squaredNorm() < minDistanceSquared )
        {
            fullPairings.emplace_back(e, ex);
            
            assignedElements[pair.first] = true;
            assignedExamples[pair.second] = true;
        }
    }
    
    for( int i = 0; i < elementsSize; ++i )
    {
        if( assignedElements[i] )
            continue;
        
        leftPartialPairings.emplace_back(elementNeighbours[i], FGraphNodeHandle::null);
    }
    
    for( int i = 0; i < examplesSize; ++i )
    {
        if( assignedExamples[i] )
            continue;
        
        rightPartialPairings.emplace_back(FGraphNodeHandle::null, exampleNeighbours[i]);
    }
    
    return result;
}

float Algorithm::_closestPoint(const std::vector<Vector3f>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
                               const std::vector<Vector3f>& bPositions, const std::vector<FGraphNodeHandle>& bNeighbours, FGraph& bGraph,
                               std::vector<std::pair<int,int>>& pairings)
{
    float cost = 0.0f;
    
    int aSize = aPositions.size();
    int bSize = bPositions.size();
    
    pairings.clear();
    pairings.reserve(aSize);
    
    float thresholdDistanceSquared = minAssignmentDistance * minAssignmentDistance;
    
    struct FastIndexSetIndex
    {
        int index;      // the index at this position
        int runStart;   // the start of a run of indices
    };
    std::vector<FastIndexSetIndex> indices(bSize + 1);
    
    for( int i =0; i < bSize + 1; ++i )
        indices[i] = {i, i};
    
    for( int aIndex = 0; aIndex < aSize; ++aIndex )
    {
        const Vector3f& aPosition = aPositions[aIndex];
        FGraphNodeHandle aHandle = aNeighbours[aIndex];
		FElementObject& aElement = aGraph.component<FElementObject>(aHandle);
        
        FGraphNodeHandle bNearestHandle;
		FElementObject * bNearest = nullptr;
        int bNearestIndex = -1;
        float minDSquared = numeric_limits<float>::max();
        
        
        int bIndex = indices[0].index;
        while( bIndex < bSize )
        {
            FGraphNodeHandle bHandle = bNeighbours[bIndex];
            const Vector3f& bPosition = bPositions[bIndex];

			FElementObject& bElement = bGraph.component<FElementObject>(bHandle);
            
            float dSquared = (aPosition - bPosition).squaredNorm();
            
            // add in the cost of unmatched types
            if( aElement.type != bElement.type )
                dSquared += typeCost * typeCost;
            
            if( dSquared < minDSquared )
            {
                minDSquared = dSquared;
                bNearestHandle = bHandle;
                bNearestIndex = bIndex;
				bNearest = &bElement;
            }
            
            bIndex = indices[bIndex + 1].index;
        }
        
        if( bNearestHandle )
        {
            // erase an index
            static FastIndexSetIndex nullIndex = {-1, -1};
            
            auto& bPrevious = bNearestIndex - 1 >= 0 ? indices[bNearestIndex - 1] : nullIndex;
            auto& bNext = indices[bNearestIndex + 1];
            auto& bCur = indices[bNearestIndex];
            
            bCur.index = bNext.index;
            
            if( bPrevious.index == bNearestIndex )
            {
                indices[bCur.runStart].index = bNext.index;
                indices[bNext.index].runStart = bCur.runStart;
            }
            else
                indices[bCur.index].runStart = bNearestIndex;
            
            pairings.emplace_back(aIndex, bNearestIndex);
            
            cost += sqrt(minDSquared);
            
            if( aElement.type != bNearest->type )
                cost += typeCost;
        }
        else
            pairings.emplace_back(aIndex, -1);
    }
    
    return cost;
}

std::vector<Algorithm::SimilarElement> Algorithm::_nearestSimilarElements
(
	FGraphNodeHandle elementHandle, Domain& elementDomain,
	std::vector<FGraphNodeHandle> candidates, Domain& candidateDomain,
	const FNeighbourhoodParameters& searchParameters,
	bool forceSpaceFilling,
	bool filterCandidateByType
)
{
    vector<Algorithm::SimilarElement> results;
    
	ElementTuple element(elementHandle, elementDomain.graph);

    auto elementNeighbours = elementDomain.neighbours(elementHandle, searchParameters.radius, searchParameters.kNearest, true, nullptr, enableVolumeSurfaceInteraction );

    auto mapping = _mapping(element, forceSpaceFilling);
	    
    for( FGraphNodeHandle candidateHandle : candidates )
    {
		ElementTuple candidate(candidateHandle, candidateDomain.graph);

        if(element.element.type != candidate.element.type && filterCandidateByType)
            continue; 
        
        mapping->setExemplarOrigin(candidate.position());
        
        auto candidateNeighbours = candidateDomain.neighbours(candidateHandle, searchParameters.radius, searchParameters.kNearest);
        
        SimilarElement similar = distance(
			element.node, elementNeighbours, elementDomain.graph, 
			candidate.node, candidateNeighbours, candidateDomain.graph, 
			mapping.get()
		);
        
        results.push_back(similar);
    }
    
    sort(results.begin(), results.end(), [](const SimilarElement& a, const SimilarElement& b) -> bool
    {
        return a.cost < b.cost;
    });
    
    return results;
}

std::vector<Algorithm::SimilarElement> Algorithm::_generation_nearestSimilarElements(FGraphNodeHandle elementHandle, Domain& elementDomain,
                                                                                     std::vector<FGraphNodeHandle> candidates, Domain& candidateDomain,
                                                                                     const FNeighbourhoodParameters& searchParameters)
{
    vector<Algorithm::SimilarElement> results;

	ElementTuple element(elementHandle, elementDomain.graph);
    
    float radius = std::max(searchParameters.radius, sourceHistogramRadius);
    
    auto neighbours = elementDomain.neighbours(elementHandle, radius, -1);
    
    auto histogram = _histogram(elementHandle, neighbours);
    
    // we only need to run one neighbourhood query, we can extract the element neighbourhood from
    // this query
    vector<FGraphNodeHandle> elementNeighbours;
    elementNeighbours.reserve(searchParameters.kNearest);
    
    float searchRadiusSqrd = searchParameters.radius * searchParameters.radius;
	Eigen::Vector3f elementPosition = element.position();
    
    for( int i = 0; i < searchParameters.kNearest && i < neighbours.size(); i++ )
    {
		FGraphNode& neighbour_i = neighbours[i](elementDomain.graph);
        float d = (eigen(neighbour_i.position) - elementPosition).squaredNorm();
        
        if( d < searchRadiusSqrd )
            elementNeighbours.push_back(neighbours[i]);
    }
    
    auto mapping = _mapping(element);
    
    for( FGraphNodeHandle candidateHandle : candidates )
    {
		ElementTuple candidate(candidateHandle, candidateDomain.graph);

		// ignore element's that don't match on id if we are painting
        if( element.element.type != candidate.element.type && generationMode != EGenerationMode::SurfacePainting )
            continue;
        
        mapping->setExemplarOrigin(candidate.position());
        
        auto candidateNeighbours = candidateDomain.neighbours(candidateHandle, searchParameters.radius, searchParameters.kNearest);
        
        SimilarElement similar = distance(
			element.node, elementNeighbours, elementDomain.graph,
			candidate.node, candidateNeighbours, candidateDomain.graph,
			mapping.get()
		);
        
        // -----------------
        // Sampling distance
        float samplingDistance = 0.0f;
        
        for( FGraphNodeHandle eHandle : candidateNeighbours )
        {
            if( histogram.find(eHandle) == histogram.end() )
                continue;

			FGraphNode& eNode = eHandle(candidateDomain.graph);
            
            float dist = (eigen(eNode.position) - candidate.position()).norm();
            
            if( dist > sourceHistogramRadius )
                continue;
            
            dist += sourceHistogramRadius - dist;
            
            samplingDistance += dist * histogram[eHandle];
        }
        
        // ---------------
        // Local Histogram
        std::map<int,int> localHistogram;
        
		for (FGraphNodeHandle eHandle : elementNeighbours)
		{
			FElementObject& eElement = elementDomain.graph.component<FElementObject>(eHandle);
			localHistogram[eElement.type]++;
		}
        
        localHistogram[element.element.type]++;
        int localHistogramCount = elementNeighbours.size() + 1 + similar.rightPartialPairings.size();
        
		for (auto& pair : similar.rightPartialPairings)
		{
			FElementObject& rightElement = candidateDomain.graph.component<FElementObject>(pair.second);
			localHistogram[rightElement.type]++;
		}
        
        float typeHistogramTerm = _compareHistograms(_exemplarTypeHistogram, _exemplar.size(),
                                                     localHistogram, localHistogramCount);
        
        similar.cost += samplingDistance * samplingDistanceWeight + typeHistogramTerm * sourceHistogramWeight;
        
        results.push_back(similar);
    }
    
    sort(results.begin(), results.end(), [](const SimilarElement& a, const SimilarElement& b) -> bool
         {
             return a.cost < b.cost;
         });
    
    return results;
}


std::vector<Algorithm::SimilarElement> Algorithm::_lsa_nearestSimilarElements(FGraphNodeHandle element, Domain& elementDomain, const FNeighbourhoodParameters& searchParameters)
{
    vector<Algorithm::SimilarElement> results;
    
    //    Eigen::VectorXf q = _discretize(element, elementDomain, 8);
    //
    //    Eigen::VectorXf q_k = q.transpose() * queryTransformMatrix;
    //    Eigen::VectorXf answer_k = lowDimensionalNeighbourhoodMatrix * q_k;
    //
    //    std::vector<size_t> indices(answer_k.size());
    //    for( int i = 0; i < indices.size(); ++i )
    //        indices[i] = i;
    //
    //    std::sort(indices.begin(), indices.end(),
    //              [&answer_k](size_t a, size_t b) {
    //                  return answer_k[a] > answer_k[b];
    //              });
    //
    //    for( size_t i : indices )
    //    {
    //        results.emplace_back();
    //
    //        SimilarElement& result = results.back();
    //
    //        result.element = _exemplar.elementAt(i);
    //        result.cost
    //
    //    }
    //
    //    vector<Algorithm::SimilarElement> results;
    
    return results;
}


// -----------------------------------------------------------------------------
/// Freespace
// -----------------------------------------------------------------------------

void Algorithm::_rebuildFreespace()
{
	_freeSpace.rebuild();
}

std::vector<FGraphNodeHandle> Algorithm::_fastFilterGenerativeExemplars(FGraphNodeHandle seedHandle, Freespace& freespace)
{
    std::vector<FGraphNodeHandle> result;

	if (!_sourceExampleMap.contains(seedHandle))
		return result;

	ElementTuple seed = outputTuple(seedHandle);

    // find the nearby freespace points
	std::vector< Eigen::Vector3f > freePoints = freespace.outputPointsInRadius(seed.position(), seed.element.generationParameters.radius );
    
    // mapping
    auto mapping = _mapping(seed);

	if(mapping->compareInExemplarSpace())
	{
		// remap the free points
		for( Eigen::Vector3f& point : freePoints )
			point = mapping->toExemplar( point );



		auto& sourceExample = _sourceExampleMap.sourceExample( seedHandle );

		auto& candidates = sourceExample.kCoherentNeighbours( _defaultSelection );

		for(FGraphNodeHandle exemplarHandle : candidates)
		{
			ElementTuple exemplarElement = exemplarTuple(exemplarHandle);

			Eigen::Vector3f exemplarElementPosition = exemplarElement.position();

			int generativeCount = 0;

			// check for overlaps with freepoints, if there are, then this exemplar will generate elements
			for(auto& freePoint : freePoints)
			{
				//auto exemplarElementPosition = exemplarElement.position();

				Eigen::Vector3f point = freePoint + exemplarElementPosition;

				auto nearest = _exemplar.nearest( point, seed.element.generationParameters.radius );

				if(nearest.element)
					generativeCount++;
			}

			if(generativeCount > 0)
				result.push_back( exemplarHandle );
		}
	}
	else
	{
		// get the relevant freespace points in the cache
		const float generationRadiusSqrd = seed.element.generationParameters.radius * seed.element.generationParameters.radius;

		auto& sourceExample = _sourceExampleMap.sourceExample( seedHandle );

		auto& candidates = sourceExample.kCoherentNeighbours( _defaultSelection );

		// now check if any free space points overlap with an exemplar
		for(FGraphNodeHandle exemplarHandle : candidates)
		{
			ElementTuple exemplar = exemplarTuple(exemplarHandle);

			mapping->setExemplarOrigin(exemplar.position());

			auto neighbours = _exemplar.neighbours( exemplarHandle, seed.element.generationParameters.radius, -1, false );

			int generativeCount = 0;

			for(auto ex_ : neighbours)
			{
				auto ex_position = exemplarGraph().node(ex_).position;
				auto exemplarPoint = mapping->toSurface( eigen(ex_position) );

				if(!exemplarPoint.hit)
					continue;

				// check for overlaps with freepoints, if there are, then this exemplar will generate elements
				for(auto& freePoint : freePoints)
				{
					float d = (freePoint - exemplarPoint.position).squaredNorm();

					if(d < generationRadiusSqrd)
					{
						generativeCount++;
						break;
					}
				}

				if(generativeCount)
					break;
			}

			if(generativeCount > 0)
				result.push_back(exemplarHandle);
		}
	}
    
    return result;
}

// -----------------------------------------------------------------------------
/// Histogram
// -----------------------------------------------------------------------------

float Algorithm::_sourceHistogramCost(const std::vector<FGraphNodeHandle>& outputElements, const std::vector<FGraphNodeHandle>& exemplarElements)
{
    map<FGraphNodeHandle, int> histogram;
     
    for( FGraphNodeHandle e : outputElements )
    {
		auto& sourceExample = _sourceExampleMap.sourceExample( e );

		FGraphNodeHandle source = sourceExample.element;

		if(!source)
			continue;

		if(histogram.find( source ) == histogram.end())
			histogram[source] = 1;
		else
			histogram[source]++;
    }
    
    float cost = 0.0f;
    
    for( FGraphNodeHandle e_ : exemplarElements )
    {
        if( histogram.find(e_) != histogram.end() )
            cost += histogram[e_];
    }
    
    return cost;
}

std::map<FGraphNodeHandle, int> Algorithm::_histogram(FGraphNodeHandle element, std::vector<FGraphNodeHandle>& neighbours)
{
    std::map<FGraphNodeHandle, int> histogram;
    
    for( FGraphNodeHandle e_ : neighbours )
    {
		if( !_sourceExampleMap.contains(e_ ) ) continue;

		auto& sourceExample = _sourceExampleMap.sourceExample( e_ );

		if(sourceExample.element)
			histogram[sourceExample.element]++;
    }
    
    return histogram;
}

float Algorithm::_compareHistograms(std::map<int,int> histogramA, int countA_in,
                                    std::map<int,int> histogramB, int countB_in)
{
    float distance = 0.0f;
    
    for( auto& pair : histogramA )
    {
        int countA = pair.second;
        int countB = histogramB.find(pair.first) == histogramB.end() ? 0 : histogramB[pair.first];
        
        float ratioA = float(countA) / float(countA_in);
        float ratioB = float(countB) / float(countB_in);
        
        distance += std::abs(ratioA - ratioB);
    }
    
    return distance;
}

// -----------------------------------------------------------------------------
/// Optimization
// -----------------------------------------------------------------------------

AlgorithmResult Algorithm::horizonOptimization()
{
    AlgorithmResult result;
    
	OptimizationProblem problem;
	problem.activeElements = _horizon;

	problem.frozenElements = _nearbyFrozenElements( problem.activeElements );
	
	{
		for(auto& e : _context.domain)
			problem.elements.insert( e.nodeHandle() );
	}
    _localOptimization( problem, result);
    
    _context.domain.rebalance();
    
    return result;
}

AlgorithmResult Algorithm::globalOptimization()
{
    AlgorithmResult result;
    
	OptimizationProblem problem;

	{
		for(auto& e : _context.domain)
			problem.elements.insert(e.nodeHandle());
	}

	problem.activeElements = problem.elements;
	problem.frozenElements = problem.elements;

	_localOptimization( problem, result );
    
    _context.domain.rebalance();
    
    return result;
}

void Algorithm::clear()
{
	clearOutput();
	clearPainting();

	_exemplar.clear();

	_defaultSelection.clear();

	_didInit = false;
}

void Algorithm::clearOutput()
{
	roundSummary().clear();

	_horizon.clear();
	_context.domain.clear();
	_sourceExampleMap.clear();

	_initFreespacePoints();
}

void Algorithm::clearPainting()
{
	_brushPoints.pts.clear();
	_initBrushIndex();
}



std::vector<FGraphNodeHandle> _filter( std::vector<FGraphNodeHandle> original, std::unordered_set<FGraphNodeHandle> toRemove )
{
	std::vector<FGraphNodeHandle> result;

	for(FGraphNodeHandle e : original)
	{
		if(toRemove.find( e ) == toRemove.end())
			result.push_back( e );
	}

	return result;
}

std::set<FGraphNodeHandle> Algorithm::_nearbyFrozenElements(const std::set<FGraphNodeHandle>& activeElements)
{
	std::set<FGraphNodeHandle> frozen;

	for(auto e : activeElements)
	{
		auto sphere = _sphereOfInfluence( e, frozenElementRadius );

		for(FGraphNodeHandle inSphere : sphere)
		{
			if(activeElements.find( inSphere ) == activeElements.end())
				frozen.insert( inSphere );
		}
	}

	return frozen;
}


std::vector<Algorithm::PredictionNeighbourhood> Algorithm::_predictionNeighbourhoods( std::set<FGraphNodeHandle>& elements )
{
	// reserve space in the results vector
	size_t n = 0;

	for(FGraphNodeHandle e : elements)
	{
		auto& sourceExample = _sourceExampleMap.sourceExample( e );

		if( sourceExample.element )
			n += 1;
	}

	std::vector<PredictionNeighbourhood> results( n );

	vector< future< void > > futures;

	size_t i = 0;
	Element::eachElement(elements, graph(), [&](FElementObject& element)
	{
		auto& sourceExample = _sourceExampleMap.sourceExample(element.nodeHandle());

		if (!sourceExample.element)
			return;

		// get the nearest element
		// Careful, this can reallocate _defaultSelection._kCoherentElements 
		// Don't do that when we are on multiple threads
		auto& coherentNeighbours = sourceExample.kCoherentNeighbours(_defaultSelection);

		futures.emplace_back(_context.threadPool->push([this, &element, &sourceExample, &coherentNeighbours, i, &results](int threadID) {
			auto nearestVec = _nearestSimilarElements(
				element.nodeHandle(), _context.domain,
				coherentNeighbours, _exemplar,
				element.optimizationParameters,
				false, // forceSpaceFilling
				false // filterCandidatesByType
			);

			PredictionNeighbourhood result;

			if (nearestVec.size() == 0)
			{
				result.element = element.nodeHandle();

				return;
			}

			result.element = element.nodeHandle();
			result.similarNeighbourhood = nearestVec[0];
			result.weight = sourceExample.weight;

			results[i] = result;
		}));

		++i;
	
	});


	for(auto& f : futures)
	{
		f.wait();
	}

	return results;
}

void Algorithm::_buildPredictions(
	const std::vector<PredictionNeighbourhood>& predictionNeighbourhoods,
	const std::set<FGraphNodeHandle>& limitPredictionsTo, 
	std::map<FGraphNodeHandle, unsigned int>& problem, 
	std::vector< std::vector< PredictedPosition > >& predictions_in_out
)
{
	for(const PredictionNeighbourhood& f : predictionNeighbourhoods)
	{
		if(!f.element)
			continue;

		if(!f.similarNeighbourhood.element)
			continue;

		FGraphNodeHandle elementHandle = f.element;
		FGraphNodeHandle exampleHandle = f.similarNeighbourhood.element;

		ElementTuple element(elementHandle, graph());
		ElementTuple example(exampleHandle, exemplarGraph());

		const size_t i = problem[elementHandle];

		auto& pairings = f.similarNeighbourhood.fullPairings;

		auto mapping = _mapping(element, example);

		for(auto& pair : pairings)
		{
			FGraphNodeHandle elementHandle_ = pair.first;
			FGraphNodeHandle exampleHandle_ = pair.second;

			ElementTuple element_(elementHandle_, graph());
			ElementTuple example_(exampleHandle_, exemplarGraph());

			// only create predictions for the active elements
			if(limitPredictionsTo.find( elementHandle_ ) == limitPredictionsTo.end())
				continue;

			if(problem.find( elementHandle_ ) == problem.end())
				continue;

			unsigned int j = problem[elementHandle_];

			// this is where we have to walk backwards on the integral curve
			// this is going to involve a shortest-path calculation :/
			auto mappedPosition = mapping->toSurface( example_.position() );
			if(!mappedPosition.hit)
				continue;

			if(!enableVolumeSurfaceInteraction)
			{
				// if one element is on a surface and the other in a volume, don't predict
				if (element.element.surfaceIndex.isOnSurface() != element_.element.surfaceIndex.isOnSurface())
					continue;
			}

			predictions_in_out[j].emplace_back( 
				i, // .elementIndex
				mappedPosition.position, // .position
				f.weight // .weight
			);
		}
	}
}



void Algorithm::_removeNonOptimizingElements( std::set<FGraphNodeHandle>& localElements )
{
	std::vector<FGraphNodeHandle> toRemove;

	for(FGraphNodeHandle fHandle : localElements)
	{
		FElementObject& f = graph().component<FElementObject>(fHandle);
		if(f.optimizationParameters.radius == 0.0f || f.optimizationParameters.kNearest == 0)
			toRemove.push_back( fHandle );
	}

	for(FGraphNodeHandle f : toRemove)
	{
		localElements.erase( f );
	}
}

std::set<FGraphNodeHandle> Algorithm::_filterSurfaceVolumeInteractions( std::set<FGraphNodeHandle>& elements )
{
	if(enableVolumeSurfaceInteraction)
		return elements;

	std::set<FGraphNodeHandle> filtered;

	if(generationMode == EGenerationMode::SurfacePainting)
	{
		for(FGraphNodeHandle eHandle : elements)
		{
			FElementObject& e = graph().component<FElementObject>(eHandle);
			if( e.surfaceIndex.isOnSurface() )
				filtered.insert( eHandle );
		}
	}
	else // volume elements
	{
		for(FGraphNodeHandle eHandle : elements)
		{
			FElementObject& e = graph().component<FElementObject>(eHandle);

			if ( !e.surfaceIndex.isOnSurface() )
				filtered.insert(eHandle);
		}
	}

	return filtered;
}

Algorithm::OptimizationProblem Algorithm::_localOptimization( OptimizationProblem problem_in, AlgorithmResult& result )
 {
	using namespace Eigen;
	using namespace std;

	typedef SparseMatrix<float> Matrix;
	typedef Triplet<float> Triplet;

	std::set<FGraphNodeHandle>& activeElements = problem_in.activeElements;
	std::set<FGraphNodeHandle>& elements = problem_in.elements;
	std::set<FGraphNodeHandle>& frozenElements = problem_in.frozenElements;

	// add the elements to the problem
	std::set<FGraphNodeHandle> localElements = activeElements;


	// add the frozen guys to the problem
	for(FGraphNodeHandle f : frozenElements)
		localElements.insert( f );

	// add remove any elements that don't have an optimization radius
	_removeNonOptimizingElements( localElements );

	// filter local elements for volume-surface interactions
	localElements = _filterSurfaceVolumeInteractions( localElements );
    
    // setup the problem map
    // we associate with each element and index in the A matrix (the problem matrix)
    std::map<FGraphNodeHandle, unsigned int> problem;
    std::vector<FGraphNodeHandle> indexToElement(localElements.size());
    
    unsigned int n = 0;
    for(auto e : localElements )
    {
        problem[e] = n;
        indexToElement[n] = e;
        ++n;
    } 

    if( n == 0 )
        return problem_in;
    
    std::vector<std::vector< PredictedPosition > > predictions(n);
    for( unsigned int i = 0; i < n; ++i )
        predictions[i] = std::vector< PredictedPosition >();
    
	for(auto keyValue : problem)
	{
		// this first prediction is a dummy for the covariance calculation
		// for some reason we do this twice, is this a mistake?
		FGraphNodeHandle handle = keyValue.first;
		unsigned int i = keyValue.second;

		auto element = outputTuple(handle);

		predictions[i].emplace_back( i, element.position() );
		predictions[i].emplace_back( i, element.position() );
	}

	// build the cache of nearest elements
	// collect the predicted positions from each frozen neighbourhood
	std::vector<PredictionNeighbourhood> predictionNeighbourhoods = _predictionNeighbourhoods( localElements );

	_buildPredictions( predictionNeighbourhoods, localElements /* limit to */, problem, predictions );


    // allocate the linear system
    // allocate Ax = b
	VectorXf b[3];
	for(int c = 0; c < 3; ++c)
	{
		b[c] = VectorXf( n );
		b[c].setZero();
	}

	std::vector<Triplet> a_t[3];
   
    // add the initial positions
	for(auto keyValue : problem)
	{
		FGraphNode& e = outputNode(keyValue.first);
		const unsigned int i = keyValue.second;

		const float weight = 0.5;

		for(int c = 0; c < 3; ++c)
		{
			a_t[c].emplace_back( i, i, weight );
			b[c]( i ) += weight * e.position[ c ];
		}
	}

    // build a variance for each updated element position
    // Compute the eigen vectors and values for the prediction clusters
    for( unsigned int i = 0; i < n; ++i )
    {
        auto& positions = predictions[i];
        
        int m = positions.size();
        
        Matrix3f cinv;
        Vector3f mean;
        
        if( m <= 1 )
        {
            cinv = Matrix3f::Identity();
            mean = Vector3f::Zero();
        }
        else
        {
            MatrixXf x(m, 3);
            
            for( int j = 0; j < m; ++j )
                x.row(j) = positions[j].position;
            
            // it's a weird eigen type
            auto mean_ = x.colwise().mean();
            MatrixXf a = x.rowwise() - mean_;
            Matrix3f c = a.transpose() * a * (1.0f / float(m - 1));
            
            cinv = pinv<Matrix3f>(c);
            mean = mean_;
        }
        
        // ignore the first item, it's a dummy for the actual position used in the covariance calculation
        float normalizeTerm = 1.0f / (positions.size() - 1);
        for( int index = 1; index < positions.size(); ++index )
        {
            auto& prediction = positions[index];
            
            int j = prediction.elementIndex;
			FGraphNode& element = outputNode(indexToElement[j]);
            
            const Vector3f predictedPosition = prediction.position;
            
            // find how the prediction position varies wrt to the covariance of
            // the other predictions
            Vector3f p = predictedPosition - mean;
            
            float w = std::exp( -.5f * p.transpose() * cinv * p ) / normalizeTerm;

			w *= prediction.weight;
            
            if( w != w || w == std::numeric_limits<float>::infinity() )
                w = 0.0f;
            
            const Vector3f predictedDirection = predictedPosition - eigen(element.position);
            
            for( int c = 0; c < 3; ++c )
            {
                a_t[c].emplace_back(i, i, w);
                a_t[c].emplace_back(j, j, w);
                a_t[c].emplace_back(i, j, -w);
                a_t[c].emplace_back(j, i, -w);
                
                b[c](j) -= w * predictedDirection(c);
                b[c](i) += w * predictedDirection(c);
            }
        }
        
        
    }
    
    // finally, we can build A
    Matrix a[3] = {Matrix(n,n), Matrix(n,n), Matrix(n,n)};
    for( int c = 0; c < 3; ++c )
    {
        a[c].setFromTriplets(a_t[c].begin(), a_t[c].end());
        a[c].finalize();
        a[c].makeCompressed();
    }
    
    VectorXf x[3] = {VectorXf(n), VectorXf(n), VectorXf(n)};
    
    for( int c = 0; c < 3; ++c )
    {
        Eigen::SimplicialCholesky<Matrix> cholesky(a[c]);
        x[c] = cholesky.solve(b[c]);
    }
    
    unordered_set<FGraphNodeHandle> modified;
    
	{
		for(auto& prediction : predictionNeighbourhoods)
		{
			if(prediction.element && prediction.similarNeighbourhood.element)
			{
				auto& sourceExample = _sourceExampleMap.sourceExample( prediction.element );

				sourceExample.element = prediction.similarNeighbourhood.element;
				sourceExample.weight = prediction.weight;
			}
		}
	}


    // update the positions
    for( auto keyValue : problem )
    {
        unsigned int i = keyValue.second;
        
        Vector3f p;
        for( int c = 0; c < 3; ++c )
            p(c) = x[c](i);
        
        auto e = outputTuple(keyValue.first);

		// update the position, p
        if( generationMode == EGenerationMode::SurfaceProjection )
        {
			auto& sourceExample = _sourceExampleMap.sourceExample( e.handle() );

			if (sourceExample.element)
			{
				auto& ex = exemplarNode(sourceExample.element);

				auto nearest = _context.meshInterface->nearestPointOnMesh( unreal( p ) );
				auto rotationAndNormal = _context.meshInterface->rotationAndNormalAtIndex( nearest.surfaceIndex );

				p = eigen( nearest.point ) + eigen( rotationAndNormal.second ) * ex.position.Z;

				FQuat quat = rotationAndNormal.first;
				FQuat rotation = quat * ex.orientation;

				e.node.orientation = rotation;
			}
        }
		else if((generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting))
		{
			auto& sourceExample = _sourceExampleMap.sourceExample(e.handle());

			if (sourceExample.element)
			{
				auto& ex = exemplarNode(sourceExample.element);

				auto nearest = _context.meshInterface->nearestPointOnMesh( unreal( p ) );
				auto rotationAndNormal = _context.meshInterface->rotationAndNormalAtIndex( nearest.surfaceIndex );

				p = eigen( nearest.point );// +eigen( rotationAndNormal.second ) * ex->position.z();

				FQuat quat = rotationAndNormal.first;
				FQuat rotation = quat * ex.orientation;

				e.node.orientation = rotation;

				e.element.surfaceIndex = nearest.surfaceIndex;
			}
		}

		if(forceRotation)
			e.node.orientation = forcedRotation;
        
        e.node.position = unreal(p);
        


        modified.insert(e.handle());
    }
    
    
    // votes and update types
	if(useTypeVoting)
	{
		// perform voting
		unordered_map<FGraphNodeHandle, unordered_map<int, int>> typeVotes;

		for(PredictionNeighbourhood& f : predictionNeighbourhoods)
		{
			auto element = outputTuple(f.element);
			
			auto& pairings = f.similarNeighbourhood.fullPairings;

			for(auto pair : pairings)
			{
				auto& other = outputElement(pair.first);

				if(!enableVolumeSurfaceInteraction)
				{
					// if one element is on a surface and the other in a volume, don't predict
					if( element.element.surfaceIndex.isOnSurface() != other.surfaceIndex.isOnSurface() )
						continue;
				}


				auto& votes = typeVotes[pair.first];
				votes[exemplarElement(pair.second).type]++;
			}
		}

		// analyze votes
		for(auto& elementVotes : typeVotes)
		{
			auto& element = outputElement(elementVotes.first);
			auto& votes = elementVotes.second;

			int maxType = 0;
			int maxCount = 0;
			int originalCount = 0;

			// determine max vote
			for(auto pair : votes)
			{
				if(pair.second > maxCount)
				{
					maxType = pair.first;
					maxCount = pair.second;
				}

				if(pair.first == element.type)
					originalCount++;
			}

			// reassign
			if(element.type != maxType && maxCount > originalCount)
			{
				element.type = maxType;
				
				modified.insert( element.nodeHandle() );

				// we changed the type, so we'll have to find a new source exemplar
				auto& sourceExample = _sourceExampleMap.sourceExample( element.nodeHandle() );

				auto& generativeByType = sourceExample.generativeElementsByType( maxType, _defaultSelection );

				auto nearest = _nearestSimilarElements( element.nodeHandle(), _context.domain, generativeByType, _exemplar, element.generationParameters );

				if(!nearest.size())
					continue;

				
				auto& sourceNode = exemplarNode(nearest[0].element);

				// assign graph objects
				// hack reassign graph objects
				// for now, we'll just copy all of the other components to this guy
				auto& elementNode = outputNode(element.nodeHandle());
				for (auto& componentType : sourceNode.components)
				{
					FGraphCopyContext::copyComponent(componentType, sourceNode, elementNode, exemplarGraph(), graph());
				}

				sourceExample.element = sourceNode.handle();
			}
		}
	}

	// rebalance the output kd-tree
	{
		std::clock_t start = std::clock();

		_context.domain.rebalance();

		double deltaTime = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		UE_LOG( LogTemp, Warning, TEXT( "Rebalanced in %f" ), (float)(deltaTime) );
	}

	// remove overlaps
	if(removeOverlaps)
	{
		std::unordered_set<FGraphNodeHandle> toUpdateFreespace;
		std::unordered_set<FGraphNodeHandle> toRemove;

		for(auto keyValue : problem)
		{
			auto element = outputTuple(keyValue.first);

			// element has been removed, don't remove his neighbours or we will have also removed the element that removed him (creating a void)
			if(toRemove.find( element.handle() ) != toRemove.end())
				continue;

			auto neighbours = _context.domain.nearestInRadius( element.position(), element.element.generationParameters.radius );

			auto overlaps = _context.domain.nearestInRadius( element.position(), element.element.radius * relaxation );

			// remove the other guys
			bool didRemove = false;

			Element::eachTuple(neighbours, graph(), [&](ElementTuple e)
			{
				if (e.handle() == element.handle())
					return;

				float allowedRadius = (element.element.radius + e.element.radius) * relaxation;

				float distance = (e.position() - element.position()).norm();
				if (distance < allowedRadius)
				{
					toRemove.insert(e.handle());
					didRemove = true;
				}
			});


			if( didRemove )
				toUpdateFreespace.insert( element.handle() );
		}

		std::vector<FGraphNodeHandle> toRemoveVec( toRemove.begin(), toRemove.end() );
		removeOutputElements( toRemoveVec );

		for(FGraphNodeHandle r : toRemoveVec)
		{
			modified.erase( r );
			activeElements.erase( r );
			elements.erase( r );
			frozenElements.erase( r );
		}

		// wake up the guys that weren't removed (but caused a removal)
		{
			std::vector<FGraphNodeHandle> asVector( toUpdateFreespace.begin(), toUpdateFreespace.end() );

			_expandFreeSpace( asVector );
		}

		// check if result contains any of the removed elements
		result.modified = _filter( result.modified, toRemove );
		result.generated = _filter( result.generated, toRemove );
		result.frozen = _filter( result.frozen, toRemove );
		result.removed = toRemoveVec;

		_context.domain.rebalance();
	}
    
    std::copy(modified.begin(), modified.end(), std::back_inserter(result.modified));

	return problem_in;
}





std::vector<FGraphNodeHandle> Algorithm::_sphereOfInfluence(FGraphNodeHandle e, float r)
{
	if(r < 0)
		r = generationInnerRadius;

    return _context.domain.neighbours(e, r, -1, false);
}

void Algorithm::_initBrushIndex()
{
	_brushPoints.pts.clear();
    
    _brushIndex.reset( new PositionFacePairCloudAdaptor( 3, _brushPoints, nanoflann::KDTreeSingleIndexAdaptorParams( 10 /* max leaf */ ) ) );
}



void Algorithm::addBrushPoint( PositionRadiusFace& point )
{
	size_t foundIndex;
	float dSqrd; 

	const float separationDistSqrd = std::pow( point.radius * 0.5f, 2.0f );

	if(_brushPoints.kdtree_get_point_count())
	{
		_brushIndex->knnSearch( &point.position.x(), 1, &foundIndex, &dSqrd );

		// don't add a brush point if it's very near another one
		if(dSqrd < separationDistSqrd)
			return;
	}


	_brushPoints.pts.push_back( point );
	size_t index = _brushPoints.kdtree_get_point_count() - 1;
	_brushIndex->addPoints( index, index );

	// wake up elements near the point
	std::vector<FGraphNodeHandle> nearest = _context.domain.nearestInRadius( point.position, point.radius * 1.5f );

	_expandFreeSpace( nearest );
	
	// add them to the horizon
	for(auto e : nearest)
		_horizon.insert( e );
}


void Algorithm::setCurrentBrushPoint( PositionRadiusFace& point )
{
	_currentBrushPoint = point;
	_hasCurrentBrushPoint = true;
}

void Algorithm::clearCurrentBrushPoint()
{
	_hasCurrentBrushPoint = false;
}


void Algorithm::clearBrushPoints()
{
	_initBrushIndex();
}

std::vector<FGraphNodeHandle> Algorithm::eraseAt( PositionRadiusFace brushPoint )
{
	_removeBrushPoint( brushPoint );

	std::vector<FGraphNodeHandle> inRadius = _context.domain.nearestInRadius( brushPoint.position, brushPoint.radius );

	if( inRadius.size() )
		this->removeOutputElements( inRadius );

	auto indices = _freeSpace.outputIndicesInRadius(brushPoint.position, brushPoint.radius);

	for (auto& pair : indices)
	{
		_freeSpace.removeOutputPoint(pair.first);
	}

	return inRadius;
}

Algorithm::ExampleSelectionPtr Algorithm::getExampleSelection()
{
	return (ExampleSelectionPtr)&_defaultSelection;
}


std::vector<FGraphNodeHandle> Algorithm::getExampleSelectionVector()
{
	return std::vector<FGraphNodeHandle>(_defaultSelection.selection.begin(), _defaultSelection.selection.end());
}

void Algorithm::updateExampleSelection( std::vector<FGraphNodeHandle> selection, float weight /*= 1.0f */ )
{
	_defaultSelection.selection = std::unordered_set<FGraphNodeHandle>( selection.begin(), selection.end() );
	_defaultSelection.weight = weight;

	_defaultSelection.init( _exemplar, *this );

	_initFreespacePoints();
}



void Algorithm::_removeBrushPoint( PositionRadiusFace& point )
{
	if(_brushPoints.kdtree_get_point_count() == 0)
		return;

	std::vector<std::pair<size_t, float>> points;
	nanoflann::SearchParams searchParams;

	_brushIndex->radiusSearch( &point.position.x(), point.radius, points, searchParams );

	if(points.size() == 0)
		return;

	std::unordered_set<size_t> toRemove;

	for(auto& pair : points)
		toRemove.insert( pair.first );

	// hack, we'll just rebuild the index
	auto pts = _brushPoints.pts;

	_brushPoints.pts.clear();
	_brushIndex.reset( new PositionFacePairCloudAdaptor( 3, _brushPoints, nanoflann::KDTreeSingleIndexAdaptorParams( 10 /* max leaf */ ) ) );

	for(size_t i = 0; i < pts.size(); ++i )
	{
		if(toRemove.find( i ) == toRemove.end())
			_brushPoints.pts.push_back( pts[i] );
	}

	if(_brushPoints.pts.size() > 0)
		_brushIndex->addPoints( 0, _brushPoints.pts.size() - 1 );
}







// Idea: keep points in the horizon that are overlapping the current brush point
//bool Algorithm::_canRemoveFromHorizon(FGraphNodeHandle handle)
//{
//	if (!_hasCurrentBrushPoint || _brushPoints.pts.size() == 0)
//		return true;
//
//	FGraphNode& element = outputNode(handle);
//	auto position = eigen(element.position);
//
//	// find nearest brush point
//	size_t foundIndex;
//	float dSqrd;
//
//	_brushIndex->knnSearch(&position.x(), 1, &foundIndex, &dSqrd);
//
//	auto& brushPoint = _brushPoints.pts[foundIndex];
//
//	float brushSqrd = std::pow(brushPoint.radius, 2.0f);
//
//	return dSqrd > brushSqrd;
//}

// Idea: keep points in the horizon that are overlapping the current brush point
bool Algorithm::_canRemoveFromHorizon( FGraphNodeHandle handle )
{
	if(!_hasCurrentBrushPoint)
		return true;

	FGraphNode& element = outputNode(handle);

	float brushSqrd = _currentBrushPoint.radius * _currentBrushPoint.radius;

	float distSqrd = (eigen(element.position) - _currentBrushPoint.position).squaredNorm();

	return distSqrd > brushSqrd;
}

auto Algorithm::_cinpact( const float u, const float k, const float c ) -> float
{
	if( u >= c || u <= -c )
		return 0.0f;

	auto result = std::exp(
		(-k * (u * u)) /
		((c * c) - (u * u)) 
	);

	return result;
}

auto Algorithm::_cinpactSum( const Eigen::Vector3f point, const std::vector<FGraphNodeHandle>& elements, const int16_t elementType, FGraph& graph ) -> float
{
	float theSum = 0.0f;

	const float k = 3.0f;

	Element::eachTuple(elements, graph, [&](ElementTuple e)
	{
		if (elementType != e.element.type)
			return;

		const float c = cellExtents;// element->radius * 2.0f;

		const float u = (e.position() - point).norm();

		theSum += _cinpact(u, k, c);
	});

	return theSum;
}

auto Algorithm::_discreteCinpactSimiliarity( const FGraphNodeHandle aHandle, Domain& aDomain, const FGraphNodeHandle bHandle, Domain& bDomain ) -> float
{
	using namespace Eigen;

	int numIndices_half = std::floor( generationParameters.radius / cellExtents );

	const Vector3f start = -Vector3f( cellExtents * numIndices_half, cellExtents * numIndices_half, cellExtents * numIndices_half );
	const Vector3f end = -start;

	Vector3f offset = start;

	float cost = 0.0f;

	FGraphNode& a = aDomain.graph.node(aHandle);
	FGraphNode& b = bDomain.graph.node(bHandle);

	// pull the types out of the generative-exemplar-by-type map
	std::vector<int16_t> types;
	for(auto& pair : _defaultSelection._generativeElementsByType)
	{
		types.push_back( pair.first );
	}

	// calculate the cost as the cinpact difference between elements of the same type
	for(offset.x() = start.x(); offset.x() <= end.x(); offset.x() += cellExtents)
	{
		for(offset.y() = start.y(); offset.y() <= end.y(); offset.y() += cellExtents)
		{
			for(offset.z() = start.z(); offset.z() <= end.z(); offset.z() += cellExtents)
			{
				Vector3f p_a = eigen(a.position) + offset;
				Vector3f p_b = eigen(b.position) + offset;

				auto aElements = aDomain.nearestInRadius( p_a, cellExtents );
				auto bElements = bDomain.nearestInRadius( p_b, cellExtents );

				for(auto t : types)
				{
					cost += std::abs( _cinpactSum( p_a, aElements, t, aDomain.graph ) - _cinpactSum( p_b, bElements, t, bDomain.graph ) );
				}
			}
		}
	}

	return cost;
}

auto Algorithm::_nearestElements_discreteCinpactSimilarity( const FGraphNodeHandle element, Domain& domain, std::vector<FGraphNodeHandle> toConsiderInExemplar, Domain& exemplar ) 
-> std::vector<SimilarElement_NoPairings>
{
	std::vector<SimilarElement_NoPairings> result;

	for(FGraphNodeHandle exemplarElement : toConsiderInExemplar)
	{
		float cost = _discreteCinpactSimiliarity( element, domain, exemplarElement, exemplar );

		SimilarElement_NoPairings similar;
		similar.element = exemplarElement;
		similar.cost = cost;

		result.push_back( similar );
	}

	sort( result.begin(), result.end(), []( const SimilarElement_NoPairings& a, const SimilarElement_NoPairings& b ) -> bool
	{
		return a.cost < b.cost;
	} );

	return result;
}

void Algorithm::ExampleSelection::clear()
{
	selection.clear();

	_kCoherentElements.clear();

	_generativeElements.clear();
	_generativeElementsByType.clear();
}

void Algorithm::ExampleSelection::init( Domain& exemplar, Algorithm& algorithm )
{
	_initGenerativeExamples( exemplar );
	_initKCoherentNeighbours( exemplar, algorithm );
}

void Algorithm::ExampleSelection::_initGenerativeExamples( Domain& exemplar )
{
	_generativeElements.clear();
	_generativeElementsByType.clear();

	Element::eachElement(selection, exemplar.graph, [&](FElementObject& element) {
		if (!element.generative)
			return;

		_generativeElements.push_back(element.nodeHandle());

		auto& byTypeVector = _generativeElementsByType[element.type];
		byTypeVector.push_back(element.nodeHandle());
	});
}

void Algorithm::ExampleSelection::_initKCoherentNeighbours( Domain& exemplar, Algorithm& algorithm )
{
	_kCoherentElements.clear();

	const float clusteringRadiusSqrd = algorithm.kCoherenceClusteringPenaltyRadius * algorithm.kCoherenceClusteringPenaltyRadius;

	std::vector< std::future< std::pair<FGraphNodeHandle, vector<FGraphNodeHandle> > > > futures;

	for(FElementObject& ex : exemplar)
	{
		// limit the operation to the selection
		if( selection.find(ex.nodeHandle()) == selection.end() ) continue;

		futures.emplace_back( algorithm._context.threadPool->push( [this, &ex, &exemplar, &algorithm, clusteringRadiusSqrd]( int threadID ) mutable {
			vector<SimilarElement> nearestExamples = algorithm._nearestSimilarElements(
				ex.nodeHandle(), exemplar, 
				_generativeElements, exemplar,
				ex.generationParameters, 
				true, // forceSpaceFilling
				false // filterCandidatesByType
			);

			int k = min<float>( nearestExamples.size(), algorithm.kCoherence );

			vector<FGraphNodeHandle> kNeighbours( k );

			for(int i = 0; i < k; ++i)
				kNeighbours[i] = nearestExamples[i].element;

			return std::make_pair( ex.nodeHandle(), kNeighbours );
		} ) );
	}

	for(auto& f : futures)
	{
		auto pair = f.get();

		_kCoherentElements[pair.first] = pair.second;
	}
}
