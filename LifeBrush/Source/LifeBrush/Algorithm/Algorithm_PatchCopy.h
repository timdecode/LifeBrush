//
//  Algorithm_PatchCopy.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2017-11-22.
//  Copyright (c) 2017 Timothy Davison. All rights reserved.
//

#pragma once

#include "Algorithm/Algorithm.h"

#include <random>

class Algorithm_PatchCopy : public Algorithm
{
protected:
	FRandomStream _randomStream;
	AABB exemplarLimits = AABB();

public:
	Algorithm_PatchCopy(SynthesisContext& context) : Algorithm(context) {}
	virtual ~Algorithm_PatchCopy() {}


	virtual AlgorithmResult generate( std::vector<PositionFace>& positions, float radius /* = -1.0f */, AABB limits /* = AABB() */ )
	{

		AlgorithmResult result;

		if(!_didInit)
		{
			_initialize();

			// don't time initialization
			beginRound();

			result = _patchInitialization( limits );
		}
		else
		{
			beginRound();

			result = globalOptimization();
		}

		endRound( result );

		return result;
	}

	virtual auto roundSummary() -> RecordSummary&
	{
		if(!_totalSummary)
			_totalSummary = std::unique_ptr<RecordSummary>( new PatchRecordSummary );

		return *(_totalSummary.get());
	}

protected:
	struct CubicPatch
	{
	public:
		Eigen::Vector3f centroid = Eigen::Vector3f::Zero();

		std::vector<FGraphNodeHandle> elements;
	};

	virtual void _initialize()
	{
		Algorithm::_initialize();

		_randomStream.Initialize( 42 );

		_updateExemplarLimits();
	}

private:
	// Based on code from [Ma et al. 2011]:
	/**
		void CParticleSystemSyn::InitializeOutputTrajectoriesViaPatchCopy3( int numOfParticles , CBoundaryConstraint* ptrBoundaryConstraint )
		{
			m_outputGroup.ClearParticleSystem();
			int inputIdx = ptrBoundaryConstraint->GetInputIdx();
			vector<CubicPatch> patches = GetCubicPatches( m_vecInputExemplar[inputIdx] ); //GetCubicPatches(m_inputSequence); //
			const int numOfPatches = int( patches.size() );
			Vec3f patchSize = CParticleSystemConfig::m_patchSize;
			Vec3f cubicMin = ptrBoundaryConstraint->GetPosMin();
			Vec3f cubicMax = ptrBoundaryConstraint->GetPosMax();
			Vec3i patchNum;
			for(int i = 0; i < 3; i++)
			{
				patchNum[i] = floor( (cubicMax[i] - cubicMin[i]) / patchSize[i] + 0.5f );
			}
			patchNum[2] = 1; // For 2D synthesis
			if(patchNum[2] == 1) cubicMin[2] = 0.0f;
			Vec3f patchCen;
			//const int numOfInputFrames = int(m_inputSequence.size());
			vector<Vec3f> vecPatchCen;
			for(int i = 0; i < patchNum[0]; i++)
			{
				patchCen[0] = cubicMin[0] + (i + 0.5f) * patchSize[0];
				for(int j = 0; j < patchNum[1]; j++)
				{
					patchCen[1] = cubicMin[1] + (j + 0.5f) * patchSize[1];
					for(int k = 0; k < patchNum[2]; k++)
					{
						patchCen[2] = cubicMin[2] + (k + 0.5f) * patchSize[2];
						vecPatchCen.push_back( patchCen );
					}
				}
			}
			random_shuffle( vecPatchCen.begin(), vecPatchCen.end() );
			for(int i = 0; i<int( vecPatchCen.size() ); i++)
			{
				patchCen = vecPatchCen[i];
				int ri = int( rand() / Flt( RAND_MAX ) * numOfPatches );
				ri = ri % numOfPatches;
				CubicPatch& patch = patches[ri];
				vector<int> indices = patch.m_patchIndices;
				int patchFrameIdx = patch.m_patchFrameIdx;
				Vec3f patchCenSrc = patch.m_patchCen;
				for(int n = 0; n<int( indices.size() ); n++)
				{
					int idxSrc = indices[n];
					//CParticleData softBodyData = m_inputSequence[patchFrameIdx].GetParticleData(idxSrc);
					CParticleData softBodyData = m_vecInputExemplar[inputIdx].GetParticleData( idxSrc );
					Vec3f center = softBodyData.GetPos() - patchCenSrc + patchCen;
					softBodyData.SetPos( center );
					softBodyData.TranslateSoftBody( patchCen - patchCenSrc );
					Vec3f pos = softBodyData.GetVecSamplePos()[0];
					if(numOfParticles > 0 && m_outputGroup.GetNumOfSoftBodies() >= numOfParticles)
					{
						break;
					}
					if(ptrBoundaryConstraint->InsideBoundaryNew( pos ) == false)
					{
						continue;
					}
					m_outputGroup.AddSoftBody( softBodyData );
				} // End-For-n
			}
			if(m_outputGroup.GetNumOfSoftBodies() < numOfParticles)
			{
				cout << "Try to initialize again...\n";
				InitializeOutputTrajectoriesViaPatchCopy3( numOfParticles, ptrBoundaryConstraint );
			}
			cout << "Have initialized " << m_outputGroup.GetNumOfSoftBodies() << " soft bodies via patch copy with a boundary constraint!\n";
		}
	*/ 
	auto _patchInitialization( AABB& limits ) -> AlgorithmResult
	{
		AlgorithmResult result;

		// we could use the generationParameters instead, but we use optimizationParameters as that is in line with Ma et al.
		auto halfCellSize = _cellSize() * 0.5f;

		// 1. Compute a vector patch center, shuffle them
		auto margin = halfCellSize * patchMarginPercentage;
		std::vector<CubicPatch> exemplarPatches = _cubicPatches_exemplar( halfCellSize, margin );

		if(exemplarPatches.size() == 0)
			return result;

		std::random_device rd;
		std::mt19937 g( rd() );

		std::shuffle( exemplarPatches.begin(), exemplarPatches.end(), g );

		// 2. Prepare the destination patches
		std::vector<Eigen::Vector3f> outputCenters = _cubicPatchCenters( halfCellSize, limits );

		// 3. Copy
		for(Eigen::Vector3f& center : outputCenters)
		{
			int32 chosen = _randomStream.RandRange( 0, exemplarPatches.size() - 1 );

			auto& patch = exemplarPatches[chosen];

			auto newElements = _cubicPatch_copyToOutput( patch, center );

			result.generated.insert( result.generated.end(), newElements.begin(), newElements.end() );
		}

		return result;
	}

	auto _cubicPatch_copyToOutput( CubicPatch& patch, Eigen::Vector3f targetCentroid ) -> std::vector<FGraphNodeHandle>
	{
		std::vector<FGraphNodeHandle> result;

		Eigen::Vector3f offset = targetCentroid - patch.centroid;

		for(FGraphNodeHandle exampleElement : patch.elements)
		{
			const auto position = exampleElement(graph()).position;

			Eigen::Vector3f newPosition = eigen(position) + offset;

			//if( _overlaps( newPosition, exampleElement->radius * relaxation ) )
			//	continue;

			FGraphNodeHandle newElement = _copyExemplarToOutput(exampleElement, exemplar().graph, newPosition, _defaultSelection);

			result.push_back( newElement );
		}

		return result;
	}

	auto _cubicPatchCenters( Eigen::Vector3f halfCellSize, const AABB limits_in ) -> std::vector<Eigen::Vector3f>
	{
		using namespace Eigen;

		auto cellSize = halfCellSize * 2.0f;

		const auto limits = _expandedLimits( limits_in, cellSize );


		std::vector<Vector3f> centers;

		Vector3f start = limits.aabb.min() + halfCellSize;
		Vector3f end = limits.aabb.max();

		Vector3f cellMin = Vector3f::Zero();


		for(cellMin.x() = start.x(); cellMin.x() < end.x(); cellMin.x() += cellSize.x())
		{
			for(cellMin.y() = start.y(); cellMin.y() < end.y(); cellMin.y() += cellSize.y())
			{
				for(cellMin.z() = start.z(); cellMin.z() < end.z(); cellMin.z() += cellSize.z())
				{
					centers.emplace_back( cellMin );
				}
			}
		}

		return centers;
	}

	auto _cubicPatches_exemplar( Eigen::Vector3f halfCellSize, Eigen::Vector3f margin ) -> std::vector<Algorithm_PatchCopy::CubicPatch>
	{
		auto centers = _cubicPatchCenters( halfCellSize, exemplarLimits );

		std::vector<CubicPatch> patches;

		for( auto& center : centers )
		{
			CubicPatch patch;
			patch.centroid = center;

			auto elementsInBox = exemplar().nearestInBox( center, halfCellSize + margin );
			patch.elements = elementsInBox;

			patches.push_back( patch );
		}

		return patches;
	}

	auto _updateExemplarLimits() -> void
	{
		Eigen::AlignedBox3f box;

		Domain& domain = this->exemplar();

		for( auto elementHandle : _defaultSelection._generativeElements)
		{
			FElementObject& elementObject = graph().component<FElementObject>(elementHandle);
			FGraphNode& elementNode = elementHandle(graph());

			Eigen::Vector3f elementHalfExtents(elementObject.radius, elementObject.radius, elementObject.radius );

			const auto p = eigen(elementNode.position);
			Eigen::AlignedBox3f elementBox(p - elementHalfExtents, p + elementHalfExtents );

			box.extend( elementBox );
		}

		exemplarLimits.aabb = box;
	}

	auto _cellSize() -> Eigen::Vector3f
	{
		Eigen::Vector3f exemplarSize = exemplarLimits.aabb.sizes();

		for(int c = 0; c < 3; c++)
			exemplarSize( c ) = exemplarSize( c ) / (patchInitializationExemplarDivisions[c] > 0 ? patchInitializationExemplarDivisions[c] : 0);

		return exemplarSize;
	}

	auto _expandedLimits(AABB limits, Eigen::Vector3f cellSize) -> AABB
	{
		auto min = limits.aabb.min();
		auto max = limits.aabb.max();
		auto size = max - min;

		auto halfSize = cellSize * 0.5f;

		for(int c = 0; c < 3; c++)
		{
			if(size( c ) < cellSize(c))
			{
				min( c ) -= halfSize(c);
				max( c ) += halfSize(c);
			}
		}

		limits.aabb = Eigen::AlignedBox3f( min, max );

		return limits;
	}
};
