// Copyright (c) 2019 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "ShipEditorSimulation/MeshSimulation.h"

#include "Simulation/FlexElements.h"
#include "Simulation/Aggregates.h"
#include "ElementEditor/RegionGrowingGenerator.h"
#include "ElementEditor/DiscreteElementEditorComponent.h"
#include "MolecularLego/MolecularLego_Relaxation.h"

#include "MolecularLegoGenerator.h"
#include "Simulation/MeshFilamentSimulation.h"

std::vector<UClass *> UMolecularLegoGenerator::dependencies()
{
	std::vector<UClass*> result = { URegionGrowingGenerator::StaticClass() };

	return result;
}

void UMolecularLegoGenerator::init(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::init(context, simulationManager);
}

void UMolecularLegoGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::attach(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;

	_simulationManager->detachAllSimulations();

	_simulationManager->registerSimulation<UMeshSimulation>();
	_simulationManager->attachSimulation<UMeshSimulation>();

	_simulationManager->registerSimulation<UMLElementSimulation>();
	_simulationManager->attachSimulation<UMLElementSimulation>();

	_simulationManager->registerSimulation<UMeshFilamentSimulation>();
	_simulationManager->attachSimulation<UMeshFilamentSimulation>();

	checkf(exampleSimulationManager(), TEXT("We need an example manager, perhaps the dependencies were not satisfied."));

	FGraph * exampleGraph = exampleSimulationManager()->graph();

	checkf(exampleGraph, TEXT("We need a graph."));

	_buildElementBVH();

	_initPath();
}

void UMolecularLegoGenerator::detach()
{
	Super::detach();

	flexSimulation->pause();

	elementBVH.removeAll();

	_initPath();
}

void UMolecularLegoGenerator::tick(float deltaT)
{
	if (!_context)
		return;

	_tick(deltaT);
}

void UMolecularLegoGenerator::_tick(float deltaT)
{
	size_t n = _brushPath.brushPoints().size();

	if (n < 2)
		return;

	UMLElementSimulation * mlSimulation = _simulationManager->simulation<UMLElementSimulation>();

	FGraph& graph = _context->graph();

	graph.beginTransaction();

	auto& elements = graph.componentStorage<FMLElement>();

	auto selection = exampleSelection();
	if (selection.Num() == 0)
		return;

	FGraphNodeHandle anyHandle = selection[0];

	const float minAssignmentDistance = mlSimulation->minAssignmentDistance;

	auto toEnumerate = _activeElements.Array();

	for (FGraphNodeHandle nodeHandle : toEnumerate)
	{
		FGraphNode& node = graph.node(nodeHandle);

		FQuat q_node = node.orientation;
		FVector p_node = node.position;
		auto i_node = elements.elementIndex(nodeHandle);

		FMLElement * element = elements.componentPtrForNode(nodeHandle);
		if( !element ) continue;

		if( !segmentForHandle.Contains(nodeHandle) ) continue;

		uint32_t segmentIndex = segmentForHandle[nodeHandle];

		int numGenerated = 0;
		// if we didn't generate anything because of being outside of the brush path, then this could still be an active element
		bool isBrushLimitied = false;

		const FVector p_a = _brushPath.brushPoints()[segmentIndex].position;
		const FVector p_b = _brushPath.brushPoints()[segmentIndex + 1].position;

		const FVector pOnPath_node = FMath::ClosestPointOnSegment(p_node, p_a, p_b);

		for (auto& rule : mlSimulation->cachedRules)
		{
			if (rule.positions.empty()) continue;

			// find the min pairings
			for (int rulePosition_i = 0; rulePosition_i < rule.positions.size(); ++rulePosition_i)
			{
				const FVector p_base = rule.positions[rulePosition_i].position;
				const FQuat qInv_base = rule.positions[rulePosition_i].orientation.Inverse();

				// only match the correct types
				if( rule.positions[rulePosition_i].type != element->type ) continue;

				for (int rulePosition_j = 0; rulePosition_j < rule.positions.size(); ++rulePosition_j)
				{
					// apply the rule into node's coordinate frame
					const FVector pred = qInv_base.RotateVector(rule.positions[rulePosition_j].position - p_base);
					const FVector p = node.orientation.RotateVector(pred) + node.position;

					const FQuat newRotation = node.orientation * (qInv_base * rule.positions[rulePosition_j].orientation);

					const FVector pOnPath = FMath::ClosestPointOnInfiniteLine(p_a, p_b, p);

					float distance = FVector::Dist(pOnPath, pOnPath_node);
					// find the sign (which direction to walk along the path)
					distance = FVector::DotProduct(pOnPath - pOnPath_node, p_b - p_a) >= 0.0f ? distance : -distance;

					LinearPath::PositionRotation frame;
	
					FVector p_;
					FQuat q_;

					if (mode == EMolecularLegotGeneratorMode::AlongPath)
					{
						if (!_brushPath.frameAlongPath(pOnPath_node, segmentIndex, distance, frame))
							continue;

						FVector offset = p - pOnPath;

						// move to identity space
						offset = _brushPath.segmentRotation(segmentIndex).UnrotateVector(offset);
						offset = frame.rotation.RotateVector(offset);

						p_ = frame.position + offset;

						// to identity
						FQuat q = _brushPath.segmentRotation(segmentIndex).Inverse() * newRotation;

						// to frame
						q_ = frame.rotation * q;
					}
					else
					{
						p_ = p;
						q_ = newRotation;
					}


					// elements must be within the bounds
					if( !_insetBounds().IsInside(p_) ) 
						continue;

					// check if we are inside the brush path
					size_t nearestSegment_a = 0;
					size_t nearestSegment_b = 0;

					const float searchRadius = maxBrushRadius;//  elementRadius + segmentLength * 2.0f;

					if (!_hasNearestBrushSegment(p_, searchRadius, nearestSegment_a, nearestSegment_b))
					{
						isBrushLimitied = true;
						continue;
					}

					// check for overlap with other elements
					{
						unrealAABB::AABB query(p_, elementRadius);

						bool hit = false;

						elementBVH.query(query, [&](unsigned int particleIndex) {
							hit = true;

							// exit the query on the first hit
							return false;
						});

						// spawn an element here
						if (!hit)
						{
							// copy one node to the end
							FGraph * exampleGraph = exampleSimulationManager()->graph();
							FGraph * targetGraph = _simulationManager->graph();


							FGraphNodeHandle newElement = _copyElement(anyHandle, *exampleGraph, *targetGraph, p_, q_);
								
							_insertElement(newElement, p_, elementRadius, frame.segment);

							numGenerated++;
						}
					}
				}
			}
		}
	}

	graph.endTransaction();
}

void UMolecularLegoGenerator::setSelection(TArray<AElementActor*> selection)
{
	URegionGrowingGenerator*  regionGrowingGenerator = elementEditorComponent->generator<URegionGrowingGenerator>();

	_selection.Empty();

	if (!regionGrowingGenerator) return;

	for (AElementActor* elementActor : selection)
	{
		FGraphNodeHandle handle = regionGrowingGenerator->handleForExampleActor(elementActor);

		_selection.Add(handle);
	}
}

bool UMolecularLegoGenerator::_hasNearestBrushSegment(FVector point, float radius, size_t& segment_a, size_t& segment_b)
{
	if (_brushPath.brushPoints().empty())
		return false;

	auto points = _brushPath.radiusSearch(point, radius);

	if (points.empty())
		return false;

	// now find the nearest brush segments
	return true;
}

void UMolecularLegoGenerator::_buildElementBVH()
{
	FGraph& graph = _context->graph();

	auto& elements = graph.componentStorage<FMLElement>();

	elementBVH.removeAll();

	for (auto& e : elements)
	{
		if (!e.isValid())
			continue;

		FGraphNode& node = graph.node(e.nodeHandle());

		if (!node.isValid())
			continue;

		elementBVH.insertParticle(node.handle().index, node.position, elementRadius);
	}
}


UGraphSimulationManager * UMolecularLegoGenerator::exampleSimulationManager()
{
	if (!_exampleSimulationManager)
	{
		URegionGrowingGenerator*  regionGrowingGenerator = elementEditorComponent->generator<URegionGrowingGenerator>();
		
		if (!regionGrowingGenerator)
			return nullptr;
		
		_exampleSimulationManager = regionGrowingGenerator->exampleSimulationManager();
	}

	return _exampleSimulationManager;
}


TArray<FGraphNodeHandle> UMolecularLegoGenerator::exampleSelection()
{
	return _selection;
}

FGraphNodeHandle UMolecularLegoGenerator::_copyElement(
	FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph, 
	FGraph& targetGraph,
	const FVector position,
	const FQuat rotation_in)
{
	FMLAggregateNO * aggregate = sourceGraph.componentPtr<FMLAggregateNO>(sourceNodeHandle);

	if (aggregate)
	{
		TArray<FGraphNodeHandle> aggregateNodes;
		TArray<FGraphEdgeHandle> aggregateEdges;

		aggregate->edgesAndNodesInAggregate(sourceGraph, aggregateNodes, aggregateEdges);

		auto edgesToCopy = sourceGraph.edgesBetweenNodes(aggregateNodes);

		FGraphNode& node = sourceGraph.node(sourceNodeHandle);

		// Grr, this BUG again:
		// We can't acces rotation_in from this constructor because unreal tries an aligned vector align!!
		// so we have to copy to a local variable to get the alignment
		FQuat vectorAlignedRotation = rotation_in;

		// position - theZeroedPosition
		const FVector translation = position - vectorAlignedRotation.RotateVector(node.position);
		// And this is where the crash happens (FTransform's constructor does an aligned vector load on the rotation)
		const FTransform transform(vectorAlignedRotation, translation);

		auto result = FGraphCopyContext::copySubgraph(sourceGraph, targetGraph, aggregateNodes, edgesToCopy, transform);

		return result.duplicatedNodes[0];
	}
	else
	{
		TArray<FGraphNodeHandle> oneNode;
		oneNode.Add(sourceNodeHandle);

		TArray<FGraphEdgeHandle> noEdges;

		const FVector translation = position - sourceGraph.node(sourceNodeHandle).position;
		const FTransform transform(rotation_in, translation);

		FGraphCopyContext result = FGraphCopyContext::copySubgraph(sourceGraph, targetGraph, oneNode, noEdges, transform);

		// the first node in here is the duplicate of the one we passed in oneNode.
		return result.duplicatedNodes[0];
	}
}

void UMolecularLegoGenerator::_insertElement(FGraphNodeHandle handle, FVector position, float radius, uint32_t segment)
{
	_activeElements.Add(handle);
	elementBVH.insertParticle(handle.index, position, radius);
	segmentForHandle.Add(handle, segment);
}

FBox UMolecularLegoGenerator::_insetBounds()
{
	FBox box = flexSimulation->instanceManagerBounds();

	float particleRadius = flexSimulation->flexParams.solidRestDistance;

	return box.ExpandBy(FVector(-particleRadius, -particleRadius, -particleRadius));
}

void UMolecularLegoGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex)
{
	FVector localPoint = flexSimulationComponent->GetOwner()->GetTransform().InverseTransformPosition(point);

	_initPath();

	_brushPath.conditionallyAddPoint(localPoint, radius, surfaceIndex);
}

void UMolecularLegoGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	FVector localPoint = flexSimulationComponent->GetOwner()->GetTransform().InverseTransformPosition(point);

	if (!_brushPath.conditionallyAddPoint(localPoint, radius, surfaceIndex))
		return;

	auto& points = _brushPath.brushPoints();

	auto selection = exampleSelection();

	if (selection.Num() == 0) return;

	FGraphNodeHandle anyHandle = selection[0];

	// add an element if we're inside the bounds
	if (points.size() == 2 && _insetBounds().IsInside(localPoint) )
	{
		FGraph * exampleGraph = exampleSimulationManager()->graph();
		FGraph * targetGraph = _simulationManager->graph();

		FQuat segmentRotation = _brushPath.segmentRotation(0);

		FGraphNodeHandle newHandle = _copyElement(anyHandle, *exampleGraph, *targetGraph, localPoint, segmentRotation);

		_insertElement(newHandle, localPoint, elementRadius, 0);
	}
}

void UMolecularLegoGenerator::endBrushPath()
{

}

void UMolecularLegoGenerator::_initPath()
{
	_brushPath.segmentLength = segmentLength;
	_brushPath.clear();

	segmentForHandle.Empty();

	_activeElements.Empty();
}
