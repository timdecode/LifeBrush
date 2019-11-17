//
//  Algorithm.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-07-20.
//  Copyright (c) 2016 Timothy Davison. All rights reserved.
//

#pragma once

#include <map>
#include <memory>
#include <set>

#include "Domain.h"
#include "SynthesisContext.h"

#include "tcods/HalfEdge.h"

#include "tcodsMeshInterface.h"
#include "Algorithm/SpaceMapping.hpp"
#include "Eigen/Geometry"
#include "Algorithm/RoundRecord.h"


#include "ElementKDTree.h"
#include "AlgorithmEnums.h"

#include "ctpl_stl.h"

struct FGraphSnapshot;

struct AlgorithmResult
{
    std::vector<FGraphNodeHandle> generated;
    std::vector<FGraphNodeHandle> modified;
    
    std::vector<FGraphNodeHandle> frozen;

	std::vector<FGraphNodeHandle> removed;
    
    void append(AlgorithmResult& other)
    {
    	// using set for faster enumeration
        std::set<FGraphNodeHandle> generatedSet(generated.begin(), generated.end());
        std::set<FGraphNodeHandle> modifiedSet(modified.begin(), modified.end());
        std::set<FGraphNodeHandle> frozenSet(frozen.begin(), frozen.end());
		std::set<FGraphNodeHandle> removedSet(removed.begin(), removed.end());
        
        for( auto e : other.generated )
            generatedSet.insert(e);
        for( auto e : other.modified )
            modifiedSet.insert(e);
        for( auto e : other.frozen )
			frozenSet.insert(e);
		for(auto e : other.removed)
			removedSet.insert( e );

		for(auto e : other.removed)
		{
			generatedSet.erase( e );
			modifiedSet.erase( e );
			frozenSet.erase( e );
		}
        
        generated = std::vector<FGraphNodeHandle>(generatedSet.begin(), generatedSet.end());
        modified = std::vector<FGraphNodeHandle>(modifiedSet.begin(), modifiedSet.end());
        frozen = std::vector<FGraphNodeHandle>(frozenSet.begin(), frozenSet.end());
		removed = std::vector<FGraphNodeHandle>( removedSet.begin(), removedSet.end() );
    }
};

struct PositionFace
{
	PositionFace() { surfaceIndex.setOffSurface(); }
	PositionFace( Eigen::Vector3f position, FSurfaceIndex surfaceIndex ) : position(position), surfaceIndex(surfaceIndex) {}

	Eigen::Vector3f position;

	FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface;
};

struct PositionRadiusFace
{
	PositionRadiusFace() {}
	PositionRadiusFace( Eigen::Vector3f position, float radius, FSurfaceIndex surfaceIndex) : position(position), radius(radius), surfaceIndex(surfaceIndex) {}

	Eigen::Vector3f position;

	FSurfaceIndex surfaceIndex = FSurfaceIndex::OffSurface;


	float& operator[](int32 i) { return position[i]; }

	float radius = 0.0f;
};

struct PositionFacePairCloud
{
	std::vector<PositionRadiusFace>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the distance between the vector "p1[0:size-1]" and the data point with index "idx_p2" stored in the class:
	inline float kdtree_distance( const float *p1, const size_t idx_p2, size_t /*size*/ ) const
	{
		const auto& p2_3f = pts[idx_p2].position;
		const auto& p1_3f = *reinterpret_cast<const Eigen::Vector3f*>(p1);

		return (p2_3f - p1_3f).squaredNorm();
	}

	inline float kdtree_get_pt( const size_t idx, int dim ) const
	{
		return pts[idx].position( dim );
	}

	template <class BBOX>
	bool kdtree_get_bbox( BBOX& /*bb*/ ) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexDynamicAdaptor<
	nanoflann::L2_Simple_Adaptor<float, PositionFacePairCloud>,
	PositionFacePairCloud,
	3 /* dim */
> PositionFacePairCloudAdaptor;



class OcclusionBase
{
public:
	virtual bool isVisible( PositionFace point, float radius )
	{
		return true;
	}
};

class Algorithm
{
public:
	struct OptimizationProblem
	{
		std::set<FGraphNodeHandle> activeElements;
		std::set<FGraphNodeHandle> frozenElements;
		std::set<FGraphNodeHandle> elements;
	};

    struct SimilarElement
    {
        FGraphNodeHandle element;
        
        std::vector<std::pair<FGraphNodeHandle, FGraphNodeHandle> > fullPairings;
        std::vector<std::pair<FGraphNodeHandle, FGraphNodeHandle> > leftPartialPairings;
        std::vector<std::pair<FGraphNodeHandle, FGraphNodeHandle> > rightPartialPairings;
        
        float cost;
    };

	struct SimilarElement_NoPairings
	{
		FGraphNodeHandle element;
		float cost = 0.0f;
	};

	struct PredictionNeighbourhood
	{
		FGraphNodeHandle element;
		SimilarElement similarNeighbourhood;

		float weight = 1.0f;
	};

	struct PredictedPosition
	{
		PredictedPosition(size_t elementIndex, Eigen::Vector3f position, float weight) : elementIndex(elementIndex), position(position), weight(weight) {}

		PredictedPosition(size_t elementIndex, Eigen::Vector3f position) : elementIndex(elementIndex), position(position), weight(1.0f) {}

		int32_t elementIndex;
		Eigen::Vector3f position;

		float weight = 1.0f;
	};
    
    struct AABB
    {
        Eigen::AlignedBox3f aabb;
        
        AABB()
        {
            const Eigen::Vector3f one = Eigen::Vector3f(1.0f, 1.0f, 1.0f);
            aabb = Eigen::AlignedBox3f(-1000000.0f * one, 1000000.0f * one);
        }
        
        bool contains(const Eigen::Vector3f& position)
        {
            return aabb.contains(position);
        }

		bool contains( const Eigen::Vector3f& position, float radius )
		{
			const Eigen::Vector3f halfExtents( radius, radius, radius );

			const Eigen::AlignedBox3f radiusBox( position - halfExtents, position + halfExtents );

			return aabb.contains( radiusBox );
		}
        
        struct Hit
        {
            float tMin;
            float tMax;
            bool hit;
        };
        
        
        Hit intersection(Eigen::Vector3f origin, Eigen::Vector3f direction)
        {
            // from http://psgraphics.blogspot.ca/2016/02/new-simple-ray-box-test-from-andrew.html
            
            Eigen::Vector3f dirInv = direction.cwiseInverse();
            
            Eigen::Vector3f t0 = (aabb.min() - origin).cwiseProduct(dirInv);
            Eigen::Vector3f t1 = (aabb.max() - origin).cwiseProduct(dirInv);
            
            Hit result;
            
            result.tMin = t0.minCoeff();
            result.tMax = t1.maxCoeff();
            result.hit = false;
            
            for( int c = 0; c < 3; c++ )
            {
                float t0_ = t0(c);
                float t1_ = t1(c);
                
                if( dirInv(c) < 0.0f )
                    std::swap(t0_, t1_);
                
                result.tMin = t0_ > result.tMin ? t0_ : result.tMin;
                result.tMax = t1_ < result.tMax ? t1_ : result.tMax;
                
                if( result.tMax <= result.tMin )
                    return result;
            }
            
            result.hit = true;
            return result;
        }
    };

	// Freespace Clustering
	struct Freespace
	{
		std::unordered_map<FGraphNodeHandle, std::vector<Eigen::Vector3f>> _freespacePoints;
		EigenPointCloud _freespaceCloud;
		EigenDynamicPointCloudIndexAdaptor* _freespaceIndex = nullptr;
		std::set<size_t> _freespaceRemovedIndices;

		size_t _pointCount = 0;

		void clear()
		{
			_freespacePoints.clear();
			_freespaceRemovedIndices.clear();
			_freespaceCloud.pts.clear();

			_pointCount = 0;

			if(_freespaceIndex)
				delete _freespaceIndex;

			_freespaceIndex = new EigenDynamicPointCloudIndexAdaptor( 3, _freespaceCloud, nanoflann::KDTreeSingleIndexAdaptorParams( 10 /* max leaf */ ) );
		}
		
		void rebuild()
		{
			std::vector<Eigen::Vector3f> points;

			points.reserve( _freespaceCloud.pts.size() - _freespaceRemovedIndices.size() );

			for(int i = 0; i < _freespaceCloud.pts.size(); ++i)
			{
				if(_freespaceRemovedIndices.find( i ) == _freespaceRemovedIndices.end())
					points.push_back( _freespaceCloud.pts[i] );
			}

			_freespaceCloud.pts.clear();

			if(_freespaceIndex)
				delete _freespaceIndex;

			_freespaceIndex = new EigenDynamicPointCloudIndexAdaptor( 3, _freespaceCloud, nanoflann::KDTreeSingleIndexAdaptorParams( 10 /* max leaf */ ) );

			_freespaceCloud.pts = std::move( points );

			_freespaceRemovedIndices.clear();

			if(_freespaceCloud.pts.size() > 0)
				_freespaceIndex->addPoints( 0, _freespaceCloud.pts.size() - 1 );
		}

		auto outputSize() -> size_t
		{
			return _pointCount;
		}

		void addOutputPoint( Eigen::Vector3f point )
		{
			size_t newIndex = _freespaceCloud.kdtree_get_point_count();

			_freespaceCloud.pts.emplace_back( point );
			_freespaceIndex->addPoints( newIndex, newIndex );

			_pointCount++;
		}

		void removeOutputPoint( size_t index )
		{
			_freespaceIndex->removePoint( index );

			_freespaceRemovedIndices.insert( index );

			_pointCount--;
		}

		auto outputPoint( size_t index ) -> Eigen::Vector3f
		{
			return _freespaceCloud.pts[index];
		}

		// Returns -1 if there isn't a nearest, otherwise the index of the nearest output freespace point.
		auto nearestInRadius(Eigen::Vector3f atPosition, float radius) -> int32
		{
			nanoflann::KNNResultSet<float, size_t> resultSet(1);

			float distSqrd = std::numeric_limits<float>::max();
			size_t index = -1;
			resultSet.init(&index, &distSqrd);

			nanoflann::SearchParams searchParams;
			searchParams.sorted = false;

			bool found = _freespaceIndex->findNeighbors(resultSet, &atPosition(0), searchParams);

			if( found && distSqrd < std::pow(radius, 2.0f) )
				return index;
			else
				return -1;
		}

		auto outputIndicesInRadius( Eigen::Vector3f atPosition, float inRadius ) -> std::vector< std::pair<size_t, float> >
		{
			std::vector<std::pair<size_t, float> > dynoResult;

			nanoflann::SearchParams searchParams;
			searchParams.sorted = false;
			float searchRadiusSqrd = inRadius * inRadius;

			_freespaceIndex->radiusSearch( &atPosition( 0 ), searchRadiusSqrd, dynoResult, searchParams );

			return dynoResult;
		}

		auto outputPointsInRadius( Eigen::Vector3f atPosition, float inRadius ) -> std::vector < Eigen::Vector3f >
		{
			auto indices = outputIndicesInRadius( atPosition, inRadius );

			std::vector < Eigen::Vector3f > points;
			points.reserve( indices.size() );

			for(auto& pair : indices)
				points.push_back( outputPoint( pair.first ) );

			return points;
		}

		auto freespacePoints( FGraphNodeHandle exemplarElement ) -> std::vector<Eigen::Vector3f>&
		{
			return _freespacePoints[exemplarElement];
		}
	};

	



public:
	Algorithm(SynthesisContext& context) : _context(context), _exemplar(_exemplarGraph)
	{
		_exemplarGraph.init();
	}

	virtual ~Algorithm() {}

	virtual void init( std::shared_ptr<ctpl::thread_pool> threadPool = nullptr );

	virtual AlgorithmResult generate( std::vector<PositionFace>& positions, float radius = -1.0f, AABB limits = AABB() );
	virtual AlgorithmResult horizonOptimization(); // local optimization (of the horizon)
	virtual AlgorithmResult globalOptimization();	

	virtual void clear();

	virtual void loadExemplar();

	// reload the internal data structures with the current state of the context
	void reloadContext();
	
	void setThreadPool(std::shared_ptr<ctpl::thread_pool> threadPool) { _context.threadPool = threadPool; }
    void setOcclusionTester(std::shared_ptr<OcclusionBase> const& occlusionTester);
	void setMeshInterface(std::shared_ptr<tcodsMeshInterface> meshInterface);
	void setToMeshTransform(const Eigen::Affine3f& transform);
	
	Domain& exemplar() { return _exemplar; }
    Domain& output() { return _context.domain; }

	std::set<FGraphNodeHandle> horizon() { return _horizon; }

	inline FGraph& graph() { return _context.graph();  }
	inline FGraph& exemplarGraph() { return _exemplar.graph; }

    
    

    

	void clearOutput();
	void clearPainting();

	void loadOutputElements(FGraphSnapshot& toLoad );

	void removeOutputElements( std::vector<FGraphNodeHandle>& elements );
	void clearAlongPath(std::vector<PositionRadiusFace>& path);

    void addToHorizon( std::vector<FGraphNodeHandle>& elements );

	void addBrushPoint( PositionRadiusFace& point );    
	void clearBrushPoints();

	std::vector<FGraphNodeHandle> eraseAt( PositionRadiusFace brushPoint );

	void setCurrentBrushPoint( PositionRadiusFace& point );
	void clearCurrentBrushPoint();

	typedef struct ExampleSelection *ExampleSelectionPtr;
	ExampleSelectionPtr getExampleSelection();
	std::vector<FGraphNodeHandle> getExampleSelectionVector();

	void updateExampleSelection( std::vector<FGraphNodeHandle> selection, float weight = 1.0f );

	SimilarElement distance( FGraphNode& element, std::vector<FGraphNodeHandle>& elementNeighbours, FGraph& elementGraph,
                             FGraphNode& example, std::vector<FGraphNodeHandle>& exampleNeighbours, FGraph& exampleGraph,
                             rg::Mapping* mapping );

	SimilarElement distance_densityBased( FGraphNode& element, std::vector<FGraphNodeHandle>& elementNeighbours, FGraph& elementGraph,
		FGraphNode& example, std::vector<FGraphNodeHandle>& exampleNeighbours, FGraph& exampleGraph,
		rg::Mapping* mapping );
	
	std::unique_ptr<RecordSummary> _totalSummary;
	virtual auto roundSummary() -> RecordSummary&
	{
		if(!_totalSummary)
			_totalSummary = std::make_unique<RecordSummary>();

		return *(_totalSummary.get());
	}

	virtual void beginRound();
	virtual void endRound(AlgorithmResult& result);

	auto totalEnergy_kCoherence( FNeighbourhoodParameters parameters ) -> double;
	auto totalEnergy_bruteForce( FNeighbourhoodParameters parameters ) -> double;

	auto sourceExamples(FGraphNodeHandle element) -> std::vector<FGraphNodeHandle>;

	inline FElementObject& outputElement(FGraphNodeHandle handle) { return graph().component<FElementObject>(handle); }
	inline FElementObject& exemplarElement(FGraphNodeHandle handle) { return exemplarGraph().component<FElementObject>(handle); }

	inline FGraphNode& outputNode(FGraphNodeHandle handle) { return graph().node(handle); }
	inline FGraphNode& exemplarNode(FGraphNodeHandle handle) { return exemplarGraph().node(handle); }

	inline ElementTuple outputTuple(FGraphNodeHandle handle) { return ElementTuple(handle, graph()); }
	inline ElementTuple exemplarTuple(FGraphNodeHandle handle) { return ElementTuple(handle, exemplarGraph()); }

public:
	template<typename MatrixType>
	static MatrixType pinv(MatrixType& inMatrix)
	{
		using namespace Eigen;
		Eigen::JacobiSVD<MatrixType> svd(inMatrix, ComputeFullU | ComputeFullV);

		auto singularValues_inv = svd.singularValues();

		double tolerance = 1.e-4;

		for (int i = 0; i < inMatrix.cols(); ++i)
			singularValues_inv(i) = singularValues_inv(i) > tolerance ? 1.0 / singularValues_inv(i) : 0.0;

		MatrixType returnValue = (svd.matrixV() * singularValues_inv.asDiagonal() * svd.matrixU().transpose());

		return returnValue;
	}

public:
    EGenerationMode generationMode = EGenerationMode::SpaceFilling;

	bool perElementParameters = false;

	bool sketchBasedForceTangentPlane = false;

	bool enableRoundSummaries = false;
	bool calculate_kCoherenceEnergy = false;
	bool calculate_bruteForceEnergy = false;
	bool useCinpactEnergy = false;
    
    float typeCost = 300.0f;

	bool useTypeVoting = false;

	bool enableReassignment = false;

	bool usePatchSelection = false;

	float gradientTerm = 1.0f;

    float minAssignmentDistance = 2.0f;
    
    float freespaceRadius = 0.0f;
    
    uint32_t kCoherence = 5;
    float kCoherenceClusteringPenalty = 0.0f;
    float kCoherenceClusteringPenaltyRadius = 5.0f;
    
    float relaxation = 0.8f;

	bool enableVolumeSurfaceInteraction = true;
    
    bool ignoreBadSuggestions = false;
	float ignoreBadSuggestionsDistanceFactor = 10.0f;

	bool forceRotation = false;
	FQuat forcedRotation = FQuat::Identity;

	bool disableOptimization = false;
    
    FNeighbourhoodParameters generationParameters = {10.0f, -1};
    FNeighbourhoodParameters optimizationParameters = {10.0f, -1};
    float generationInnerRadius = 50000.0f;
	float frozenElementRadius = 0.0f;

    int seedSeparation = 6;

	bool seedsIgnoreMeshInBoundary = false;

	int optimizationRounds = 1;
	bool useGlobalOptimization = false;
       
    float sourceHistogramWeight = 1.0f;
    float sourceHistogramRadius = 5.0f;
    
    float samplingDistanceWeight = 1.0f;
    
    bool flipSurfaceNormals = false;

	bool removeOverlaps = false;

	FIntVector patchInitializationExemplarDivisions = FIntVector( 1, 1, 1 ); // a huge polymorphism hack, i know
	float patchMarginPercentage = 0.0f;

	int32 nThreads = 2;

	float voidSize = 10.0f;

	float cellExtents = 1.0f;

protected:
	void _conditionally_initializeThreadPool();
	virtual void _initialize();
	void _initElementParameters();

	void _initDefaultSelection();

	void _initExemplarTypeHistogram();
    void _initFreespacePoints();
	void _initLSA();

	void _removeBrushPoint( PositionRadiusFace& point );

    
	// forceSpaceFilling is an option so that we can setup the KCoherentNeighbours with a normal mapping
    auto _mapping(Eigen::Vector3f const & surfaceOrigin, Eigen::Vector3f const & exemplarOrigin, bool forceSpaceFilling = false ) -> std::unique_ptr<rg::Mapping>
    {
		if(generationMode == EGenerationMode::SpaceFilling || generationMode == EGenerationMode::SpacePainting || forceSpaceFilling)
			return std::unique_ptr<rg::Mapping>( new rg::SpaceMapping( _context.meshInterface.get(), surfaceOrigin, exemplarOrigin ) );
		else if(generationMode == EGenerationMode::SurfaceProjection || sketchBasedForceTangentPlane)
			return std::unique_ptr<rg::Mapping>( new rg::SurfaceMapping( _context.meshInterface.get(), surfaceOrigin, exemplarOrigin ) );
		else
			return std::unique_ptr<rg::SurfaceWalk>( new rg::SurfaceWalk( _context.meshInterface.get(), surfaceOrigin, exemplarOrigin ) );
    }

	auto _mapping(ElementTuple element, ElementTuple example, bool forceSpaceFilling = false ) -> std::unique_ptr<rg::Mapping>
	{
		if((generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting) 
			&& !(forceSpaceFilling || sketchBasedForceTangentPlane) )
			return std::unique_ptr<rg::SurfaceWalk>( new rg::SurfaceWalk( _context.meshInterface.get(), element.position(), element.element.surfaceIndex, example.position() ) );
		else
			return _mapping( element.position(), example.position() );
	}

	auto _mapping( ElementTuple element, bool forceSpaceFilling = false ) -> std::unique_ptr<rg::Mapping>
	{
		if((generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting) && !(forceSpaceFilling|| sketchBasedForceTangentPlane))
			return std::unique_ptr<rg::SurfaceWalk>( new rg::SurfaceWalk( _context.meshInterface.get(), element.position(), element.element.surfaceIndex, Eigen::Vector3f::Zero() ) );
		else
			return _mapping(element.position(), Eigen::Vector3f::Zero(), forceSpaceFilling );
	}

	auto _mapping( Eigen::Vector3f & surfaceOrigin, FSurfaceIndex surfaceIndex, Eigen::Vector3f exemplarOrigin = Eigen::Vector3f::Zero() )->std::unique_ptr<rg::Mapping>
	{
		if((generationMode == EGenerationMode::SurfaceWalking || generationMode == EGenerationMode::SurfacePainting) && !sketchBasedForceTangentPlane)
			return std::unique_ptr<rg::SurfaceWalk>( new rg::SurfaceWalk( _context.meshInterface.get(), surfaceOrigin, surfaceIndex, exemplarOrigin ) );
		else
			return _mapping( surfaceOrigin, exemplarOrigin );
	}


    // Initializes the output domain at the given position by copying elements from
    // the exemplar to the output domain.
    // not called by _initialize()
    void _initializeOutput(PositionFace& position, AABB& limits, AlgorithmResult& result);
    
    /*
     Returns the elements that can still generate new elements. The elements that can no longer generate new elements
     are placed in the frozen_out vector.
     \param elements The elements to filter.
     \param domain The domain to search for neighbours in.
     \param frozen_out The elements removed from elements.
     \return The remaining elements that can predict new elements.
     */
    std::set<FGraphNodeHandle> _removeFrozen(const std::set<FGraphNodeHandle>& elements, Domain& domain, std::vector<FGraphNodeHandle>& frozen_out);
    std::set<FGraphNodeHandle> _removeFrozen_usingFreespacePoints(const std::set<FGraphNodeHandle>& elements, Domain& domain, std::vector<FGraphNodeHandle>& frozen_out);

    std::vector<FGraphNodeHandle> _seeds(const std::set<FGraphNodeHandle>& elements);

    float _closestPoint(const std::vector<Eigen::Vector3f>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
                        const std::vector<Eigen::Vector3f>& bPositions, const std::vector<FGraphNodeHandle>& bNeighbours, FGraph& bGraph,
                        std::vector<std::pair<int,int>>& pairings_out);
    
    float _sourceHistogramCost(const std::vector<FGraphNodeHandle>& outputElements, const std::vector<FGraphNodeHandle>& exemplarElements);
    std::map<FGraphNodeHandle, int> _histogram(FGraphNodeHandle element, std::vector<FGraphNodeHandle>& neighbours);
    float _compareHistograms(std::map<int,int> histogramA, int countA,
                             std::map<int,int> histogramB, int countB);

	// forceSpaceFilling is an option for setting up the k-nearest neighbours in the exemplar
    std::vector<Algorithm::SimilarElement> _nearestSimilarElements(
		FGraphNodeHandle element, Domain& elementDomain,
        std::vector<FGraphNodeHandle> candidates, Domain& candidateDomain,
        const FNeighbourhoodParameters& searchParameters,
		bool forceSpaceFilling = false,
		bool filterCandidatesByType = true
	);
    
    std::vector<Algorithm::SimilarElement> _generation_nearestSimilarElements(FGraphNodeHandle element, Domain& elementDomain,
                                                                              std::vector<FGraphNodeHandle> candidates, Domain& candidateDomain,
                                                                              const FNeighbourhoodParameters& searchParameters);


    void _orderNearest(std::vector<FGraphNodeHandle>& elements, FGraph& graph, const Eigen::Vector3f referencePoint);
    
    

	enum PredictionType : uint8_t 
	{
		FrozenElements,
		ActiveElements
	};

	OptimizationProblem _localOptimization( OptimizationProblem problem, AlgorithmResult& result);

	void _removeNonOptimizingElements( std::set<FGraphNodeHandle>& localElements );

	std::set<FGraphNodeHandle> _filterSurfaceVolumeInteractions( std::set<FGraphNodeHandle>& elements );

	
	std::set<FGraphNodeHandle> _nearbyFrozenElements( const std::set<FGraphNodeHandle>& activeElements );
    
    bool _overlaps(const Eigen::Vector3f& position, const float radius, FNeighbourhoodParameters& parameters, const FGraphNodeHandle elementToIgnore = FGraphNodeHandle());
    bool _inBoundary(const Eigen::Vector3f& position, const float radius);
    
    std::vector<FGraphNodeHandle> _sphereOfInfluence(FGraphNodeHandle e, float r = -1.0);
    
    float _maxRadius = 0.0f;
    float _computeMaxRadius();
    

	// Latent Semantic Analysis Idea
    void _sampleNeighbourhood(FGraphNodeHandle element, Domain& domain, Eigen::VectorXf& vector_out);
    auto _discretize(FGraphNodeHandle element, Domain& domain, const int n) -> Eigen::VectorXf;
    
	// Freespace
    void _rebuildFreespace();
	std::vector<FGraphNodeHandle> _fastFilterGenerativeExemplars( FGraphNodeHandle seed, Freespace& freespace );
    
    void _expandFreeSpace(std::vector<FGraphNodeHandle>& generated);


	void _initBrushIndex();

protected:
	struct ExampleSelection
	{
	public:
		std::unordered_set<FGraphNodeHandle> selection;

		// k-coherence search
		std::unordered_map<FGraphNodeHandle, std::vector<FGraphNodeHandle> > _kCoherentElements;

		std::unordered_map<int16, std::vector<FGraphNodeHandle>> _generativeElementsByType;

		std::vector<FGraphNodeHandle> _generativeElements;

		float weight = 1.0f;

		void clear();

		void init( Domain& exemplar, Algorithm& algorithm );

	protected:
		void _initGenerativeExamples( Domain& exemplar );

		void _initKCoherentNeighbours( 
			Domain& exemplar, 
			Algorithm& algorithm
		);
	};

	friend ExampleSelection;

	Algorithm::ExampleSelection _defaultSelection;

	Algorithm::ExampleSelection& _highestWeightExampleSelection();

protected:
	struct SourceExampleMap
	{
	public:
		struct SourceExample
		{
			SourceExample() {}
			SourceExample( FGraphNodeHandle elementHandle, float weight ) :
				element(elementHandle), weight( weight ) {}

			FGraphNodeHandle element;
			float weight = 0.0f;

			std::vector<FGraphNodeHandle>& kCoherentNeighbours( ExampleSelection& defaultSelection )
			{
				return defaultSelection._kCoherentElements[element];
			}

			std::vector<FGraphNodeHandle>& generativeElementsByType( int16 type, ExampleSelection& defaultSelection )
			{
				return defaultSelection._generativeElementsByType[type];
			}
		};

		bool contains(FGraphNodeHandle outputElement)
		{
			return _sourceExample.find(outputElement) != _sourceExample.end();
		}

		SourceExample& sourceExample(FGraphNodeHandle outputElement )
		{
			return _sourceExample[outputElement];
		}

		void eraseElement(FGraphNodeHandle outputElement )
		{
			auto found = _sourceExample.find( outputElement );

			if(found != _sourceExample.end())
				_sourceExample.erase( found );
		}

		void clear()
		{
			_sourceExample.clear();
		}

	protected:
		std::unordered_map<FGraphNodeHandle, SourceExample > _sourceExample; // key = output element, value = the example element that create the output element
	};

protected:


	std::vector<PredictionNeighbourhood> _predictionNeighbourhoods( std::set<FGraphNodeHandle>& elements );

	void _buildPredictions(
		const std::vector<PredictionNeighbourhood>& predictionNeighbourhoods,
		const std::set<FGraphNodeHandle>& excludeFromPredictions,
		std::map<FGraphNodeHandle, unsigned int>& problem,
		std::vector< std::vector< PredictedPosition > >& predictions_in_out
	);

protected:
	FGraphNodeHandle _copyExemplarToOutput( FGraphNodeHandle sourceElement, FGraph& sourceGraph, const Eigen::Vector3f& position, Algorithm::ExampleSelection& selection, FQuat q = FQuat::Identity );


protected:
    bool _didInit = false;
    
	SynthesisContext& _context;

	FGraph _exemplarGraph;
    Domain _exemplar;
    
    std::map<int,int> _exemplarTypeHistogram; // a histogram of elementID counts in the exemplar
    
    Eigen::Affine3f _toMeshTransform = Eigen::Affine3f(Eigen::Affine3f::Identity());
    
    std::set<FGraphNodeHandle> _horizon;

	SourceExampleMap _sourceExampleMap;

    // LSA
    Eigen::MatrixXf queryTransformMatrix;
    Eigen::MatrixXf lowDimensionalNeighbourhoodMatrix;
    std::vector<SimilarElement> _lsa_nearestSimilarElements(FGraphNodeHandle element, Domain& elementDomain, const FNeighbourhoodParameters& searchParameters);
    
    

	Freespace _freeSpace;

    std::shared_ptr<OcclusionBase> _occlusionTester;

	// Painting
	PositionFacePairCloud _brushPoints;
	std::unique_ptr<PositionFacePairCloudAdaptor> _brushIndex;

	bool _hasCurrentBrushPoint = false;
	PositionRadiusFace _currentBrushPoint;

	bool _canRemoveFromHorizon( FGraphNodeHandle element );

protected:
	// CINPACT energy
	// --------------
	auto _cinpact( const float u, const float k, const float c ) -> float;

	auto _cinpactSum( const Eigen::Vector3f point, const std::vector<FGraphNodeHandle>& elements, const int16_t elementType, FGraph& graph ) -> float;
	auto _discreteCinpactSimiliarity( const FGraphNodeHandle a, Domain& aDomain, const FGraphNodeHandle b, Domain& bDomain ) -> float;

	auto _nearestElements_discreteCinpactSimilarity( const FGraphNodeHandle element, Domain& domain, std::vector<FGraphNodeHandle> toConsiderInExemplar, Domain& exemplar )->std::vector<SimilarElement_NoPairings>;
};


