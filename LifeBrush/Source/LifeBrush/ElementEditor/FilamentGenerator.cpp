//
//  Created by Timothy Davison on 2019-10-20.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "ElementEditor/SwarmGenerator.h"

#include "ShipEditorSimulation/MeshSimulation.h"
#include "Simulation/MeshFilamentSimulation.h"
#include "Simulation/Brownian.h"
#include "MolecularLego/MolecularLego_Relaxation.h"

#include "FilamentGenerator.h"

void UFilamentGenerator::init(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::init(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;

}

void UFilamentGenerator::attach(SynthesisContext * context, UGraphSimulationManager * simulationManager)
{
	Super::attach(context, simulationManager);

	_context = context;
	_simulationManager = simulationManager;

	_simulationManager->detachAllSimulations();

	_simulationManager->registerSimulation<UMeshSimulation>();
	_simulationManager->attachSimulation<UMeshSimulation>();

	_simulationManager->registerSimulation<UMeshFilamentSimulation>();
	_simulationManager->attachSimulation<UMeshFilamentSimulation>();

	// for statics
	_simulationManager->registerSimulation<UStaticPositionSimulation>();
	_simulationManager->attachSimulation<UStaticPositionSimulation>();

	_simulationManager->registerSimulation<UMLElementSimulation>();
	_simulationManager->attachSimulation<UMLElementSimulation>();

	_simulationManager->registerSimulation<USwarmSimulation>();
	_simulationManager->attachSimulation<USwarmSimulation>();

	//_simulationManager->registerSimulation<USwarmSimulation>();
	//_simulationManager->attachSimulation<USwarmSimulation>();

	//_simulationManager->registerSimulation<USingleParticleBrownianSimulation>();
	//_simulationManager->attachSimulation<USingleParticleBrownianSimulation>();

	//flexSimulation->play();

	rand.GenerateNewSeed();

	_readSegments();
}

void UFilamentGenerator::detach()
{
	Super::detach();

	flexSimulation->pause();
}

void UFilamentGenerator::tick(float deltaT)
{
	if (_brushPath.brushPoints().size() == 0) return;

	FGraph * exampleGraph = _exampleGraph();
	FGraph * targetGraph = _simulationManager->graph();

	if (!exampleGraph || !targetGraph) return;

	FGraphNodeHandle nextHandle = _nextPrototypeHandle();

	targetGraph->beginTransaction();
	
	if (nextHandle)
	{
		float radius = _prototype.radius;

		float nextLength = _sequenceLength(_filament.Num() + 1);

		// we can add an element
		if (_brushPath.length() > nextLength)
		{
			FVector startPoint = _brushPath.brushPoints()[0].position;

			LinearPath::PositionRotation frame;

			if (_brushPath.frameAlongPath(startPoint, 0, nextLength, frame))
			{
				FGraphNode& exampleNode = exampleGraph->node(nextHandle);

				FQuat exampleRotation = exampleNode.orientation;
				
				//if (exampleNode.hasComponent<FSGStarRule>())
				//{
				//	FVector forward = frame.rotation.RotateVector(FVector::ForwardVector);

				//	FQuat randomRotation = FQuat(forward, rand.GetFraction() * M_PI * 2.0f);

				//	exampleRotation = randomRotation * exampleRotation;
				//}

				FQuat rotation = frame.rotation * exampleRotation;


				FGraphNodeHandle rootNode = _filament.Num() > 0 ? _filament[0] : FGraphNodeHandle::null;

				FGraphNodeHandle handle = _copyElement(nextHandle, *exampleGraph, *targetGraph, frame.position, rotation, rootNode);

				_filament.Add(handle);
				
				_filamentOriginalRotations.Add(exampleRotation);

				int filamentIndex = _filament.Num() - 1;

				// create connections to previous
				if (filamentIndex > 0)
				{
					FGraphNodeHandle last = _filament[filamentIndex - 1];

					// add an FFilamentConnection
					{
						auto edgeHandle = targetGraph->connectNodes<FFilamentConnection>(last, handle);
						FFilamentConnection& filament = targetGraph->edgeObject<FFilamentConnection>(edgeHandle);

						filament.group = _filamentGroup;
						filament.segmentID = filamentIndex = filamentIndex - 1;
						filament.radius = radius;
					}

					// add a FFlexConnection
					{
						auto edgeHandle = targetGraph->connectNodes<FFlexConnection>(last, handle);
						FFlexConnection& flexConnection = targetGraph->edgeObject<FFlexConnection>(edgeHandle);

						flexConnection.coefficient = 0.5f;
						flexConnection.length = _prototype.segmentLength;
					}

					// doubly linked connections
					if( filamentIndex - 2 >= 0 )
					{ 
						FGraphNodeHandle lastLast = _filament[filamentIndex - 2];

						auto edgeHandle = targetGraph->connectNodes<FFlexConnection>(lastLast, handle);
						FFlexConnection& flexConnection = targetGraph->edgeObject<FFlexConnection>(edgeHandle);

						flexConnection.coefficient = 0.5f;
						flexConnection.length = _prototype.segmentLength * 2.0f;
					}
				}
			}
		}
	}


	_updateFilamentPositions();
	_updateConnections();

	targetGraph->endTransaction();
}

void UFilamentGenerator::tickPaused(float deltaT)
{
	tick(deltaT);
}

float UFilamentGenerator::_correctedSegmentLength(float length)
{
	float radius = length * 0.5f;

	return std::max(flexSimulation->flexParams.solidRestDistance, radius) * 2.0f;
}

FGraphNodeHandle UFilamentGenerator::_nextPrototypeHandle()
{
	int n = _filament.Num();

	if (n < _headSequence.Num())
	{
		return _headSequence[n];
	}
	else if (n < _headSequence.Num() + _prototype.numSegments)
	{
		auto handle = _bodySequence[(n - _headSequence.Num()) % _bodySequence.Num()];

		return handle;
	}
	else if(n < _headSequence.Num() + _prototype.numSegments + _tailSequence.Num())
	{
		int m = (n - (_headSequence.Num() + _prototype.numSegments));

		return _tailSequence[m];
	}

	return FGraphNodeHandle::null;
}


FBox UFilamentGenerator::_insetBounds()
{
	FBox box = flexSimulation->instanceManagerBounds();

	float particleRadius = flexSimulation->flexParams.solidRestDistance;

	return box.ExpandBy(FVector(-particleRadius, -particleRadius, -particleRadius));
}



void UFilamentGenerator::_updateFilamentPositions()
{
	if (_brushPath.brushPoints().size() == 0)
		return;

	FGraph * graph = _simulationManager->graph();

	float t = 0.0f;

	FVector lastPosition = _brushPath.brushPoints()[0].position;
	uint32_t lastSegment = 0;

	bool isFirst = true;

	// idea:
	// walk along the path by segment radii so that segments are separated by the radii of the current
	// and previous elements
	int i = 0;
	for (FGraphNodeHandle handle : _filament)
	{
	
		FGraphNode& node = graph->node(handle);

		FQuat originalRotation = _filamentOriginalRotations[i];

		LinearPath::PositionRotation frame;

		float distance = 0.0f;

		if (isFirst)
		{
			distance = 0.0f;
			isFirst = false;
		}
		else
			distance = _prototype.segmentLength;

		if (_brushPath.frameAlongPath(lastPosition, lastSegment, distance, frame))
		{
			node.position = frame.position;
			node.orientation = frame.rotation * originalRotation;

			lastPosition = frame.position;
			lastSegment = frame.segment;
		}
		else // put the rest of the nodes at the end of the path
		{
			node.position = _brushPath.brushPoints().back().position;
			node.orientation = _brushPath.segmentRotation(_brushPath.brushPoints().size() - 2);
		};

		++i;
	}
}

void UFilamentGenerator::_updateConnections()
{
	FGraph * graph = _simulationManager->graph();

	if (_filament.Num() < 2) return;

	// update filaments
	//for (int i = 0; i < _filament.Num(); ++i)
	//{
	//	FGraphNodeHandle curHandle = _filament[i];
	//	FGraphNodeHandle previousHandle = i > 0 ? _filament[i - 1] : FGraphNodeHandle::null;

	//	FGraphNode& curNode = graph->node(curHandle);

	//	curNode.each<FFilamentConnection>(*graph, [&](FGraphNodeHandle other, FFilamentConnection& connection) {
	//		if (other != previousHandle)
	//		{
	//			connection.segmentID = i;
	//			connection.group = _filamentGroup;
	//		}
	//	});
	//}

	// update flex connections
	for (FGraphNodeHandle handle : _filament)
	{
		FGraphNode& node = graph->node(handle);

		graph->node(handle).each<FFlexConnection>(*graph, [&](FGraphNodeHandle other, FFlexConnection& connection)
		{
			FGraphNode& otherNode = graph->node(other);

			connection.length = FVector::Dist(node.position, otherNode.position);
		});
	}
}



void UFilamentGenerator::beginBrushPath(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	// start the brush path
	FVector localPoint = flexSimulationComponent->GetOwner()->GetTransform().InverseTransformPosition(point);
	
	_initPath();

	// only begin the path if we are inside the bouds
	if (_insetBounds().IsInside(localPoint))
		_beginBrushPath(localPoint, radius, surfaceIndex);
	else
		_firstOutside = true;
}


void UFilamentGenerator::_beginBrushPath(FVector localPoint, float radius, FSurfaceIndex surfaceIndex)
{
	_brushPath.conditionallyAddPoint(localPoint, radius, surfaceIndex);
}

void UFilamentGenerator::addBrushPoint(FVector point, float radius, FSurfaceIndex surfaceIndex /*= FSurfaceIndex::OffSurface*/)
{
	FVector localPoint = flexSimulationComponent->GetOwner()->GetTransform().InverseTransformPosition(point);

	if (_forceEndPath)
		return;

	bool isInside = _insetBounds().IsInside(localPoint);

	// do we still need to start the path?
	if (_brushPath.brushPoints().empty())
	{
		if( isInside )
			_beginBrushPath(localPoint, radius, surfaceIndex);
	}
	else
	{
		if (isInside)
			_brushPath.conditionallyAddPoint(localPoint, radius, surfaceIndex);
		// end the path
		else
		{
			// add a static element to fix the position if we have more to add
			if (_filament.Num() > 0 && _nextPrototypeHandle())
			{
				auto last = _filament.Last();

				FGraph * targetGraph = _simulationManager->graph();

				FGraphNode& lastNode = targetGraph->node(last);

				FStaticPositionObject& staticPosition = lastNode.addComponent<FStaticPositionObject>(*targetGraph, lastNode.position, lastNode.orientation);
			}

			_initPath();
		}
	}

	// add a static element to fix the position if we have more to add
	if (_filament.Num() > 0 && _firstOutside )
	{
		auto first = _filament[0];

		FGraph * targetGraph = _simulationManager->graph();

		FGraphNode& firstNode = targetGraph->node(first);

		FStaticPositionObject& staticPosition = firstNode.addComponent<FStaticPositionObject>(*targetGraph, firstNode.position, firstNode.orientation);
	}
}

void UFilamentGenerator::endBrushPath()
{
	clearFilamentPrototype();

	//_filament.Empty();
	//_filamentOriginalRotations.Empty();
}

TArray<FGraphNodeHandle> UFilamentGenerator::_parseSequence(FString sequence)
{
	TArray<FString> headNames;
	sequence.ParseIntoArray(headNames, TEXT(","));

	TArray<FGraphNodeHandle> handles;

	for (FString& stringName : headNames)
	{
		FName name(*stringName);

		FGraphNodeHandle * headPrototype = _segmentTypeToNode.Find(name);

		if (!headPrototype) continue;

		handles.Add(*headPrototype);
	}

	return handles;
}

void UFilamentGenerator::_initPath()
{
	_filament.Empty();
	_filamentOriginalRotations.Empty();
	
	_forceEndPath = false;
	_firstOutside = false;

	// each filament needs a group
	UMeshFilamentSimulation * filamentSimulation = flexSimulation->simulationManager.simulation<UMeshFilamentSimulation>();
	if (filamentSimulation)
	{
		_filamentGroup = filamentSimulation->nextGroup(_prototype.bodyMaterial);
	}

	_brushPath.segmentLength = brushSegmentLength;
	_brushPath.clear();
}

float UFilamentGenerator::_sequenceLength(int32 numInSequence)
{
	if (numInSequence <= 1)
		return 0.0f;

	return (numInSequence - 1) * _prototype.segmentLength;
}


FGraphNodeHandle UFilamentGenerator::_copyElement(
	FGraphNodeHandle sourceNodeHandle, FGraph& sourceGraph,
	FGraph& targetGraph,
	const FVector position,
	const FQuat rotation_in,
	FGraphNodeHandle rootNode)
{
	USwarmSimulation * swarmSimulation = _simulationManager->simulation<USwarmSimulation>();

	return swarmSimulation->copyElement(sourceNodeHandle, sourceGraph, targetGraph, position, rotation_in, rootNode);
}




void UFilamentGenerator::setFilamentPrototype(FFilamentPrototype filamentPrototype)
{
	clearFilamentPrototype();

	_prototype = filamentPrototype;
	_hasPrototype = true;

	// Correct the length based on the larger of the segment length or the Flex particle radius*2, to avoid physics issues.
	_prototype.segmentLength = _correctedSegmentLength(_prototype.segmentLength);

	_readSegments();
}

void UFilamentGenerator::clearFilamentPrototype()
{
	_hasPrototype = false;
	_firstOutside = false;

	_filament.Empty();
	_filamentOriginalRotations.Empty();
	_brushPath.clear();
}

void UFilamentGenerator::togglePlay()
{
	if (flexSimulation->isPlaying())
		flexSimulation->pause();
	else
	{
		flexSimulation->updateFlexState();

		flexSimulation->begin();

		flexSimulation->play();
	}
}

void UFilamentGenerator::_readSegments()
{
	FGraph * exampleGraph = _exampleGraph();

	if (!exampleGraph) return;

	auto& segmentPrototypes = exampleGraph->componentStorage<FSGRuleTypeName>();

	_segmentTypeToNode.Empty();

	for (FSGRuleTypeName& segment : segmentPrototypes)
	{
		_segmentTypeToNode.Add(segment.typeName, segment.nodeHandle());
	}

	// cache the sequences into node handles
	_headSequence = _parseSequence(_prototype.headSequence);
	_tailSequence = _parseSequence(_prototype.terminusSequence);
	_bodySequence = _parseSequence(_prototype.bodySequence);
}

FGraph * UFilamentGenerator::_exampleGraph()
{
	return &_simulationManager->simulation<USwarmSimulation>()->ruleGraph;
}
