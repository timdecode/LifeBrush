// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Visualization/Timeline.h"

#include "EventVisualizationTool.h"



void UEventVisualizationTool::gainFocus()
{
	Super::gainFocus();

	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	timeline->setGlyphVisibility(true);
}

void UEventVisualizationTool::loseFocus()
{
	Super::loseFocus();

	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	timeline->setGlyphVisibility(false);
}

void UEventVisualizationTool::oneHandStart(UPrimitiveComponent * hand)
{
	Super::oneHandStart(hand);

	_selection.Empty();
}

void UEventVisualizationTool::oneHandEnd(UPrimitiveComponent * hand)
{
	Super::oneHandEnd(hand);

	if( _selection.Num() )
		_traceSelection(_selection);
	else
	{
		UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

		int32 currentFrame = _flexSimulation->graphSimulation.tickCount;

		// show the last 10 seconds of events
		timeline->showAllEvents();

		UVisualization_AgentPathLines * pathlines = _flexSimulation->simulationManager.simulation<UVisualization_AgentPathLines>();

		pathlines->hideTotalHistory();
	}
}

void UEventVisualizationTool::tick(float dt)
{

}

void UEventVisualizationTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	Super::tickOneHand(dt, hand, lastTransform);

	_selectEvent(dt, hand, lastTransform);
}

void UEventVisualizationTool::faceDown_released()
{
	if (_flexSimulation->isPlaying())
		_flexSimulation->pause();
	else
		_flexSimulation->play();
}

void UEventVisualizationTool::faceUp_released(USceneComponent * interactionPoint /* = nullptr */)
{
	physicalInteraction = !physicalInteraction;
}

void UEventVisualizationTool::_selectEvent(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	FTransform toLocal = _flexSimulation->owner->GetRootComponent()->GetComponentTransform();

	FVector localHandPosition = toLocal.InverseTransformPosition(hand->GetComponentLocation());

	auto overlappingEvents = timeline->eventsOverlappingPosition(localHandPosition, _brushRadius());

	for (USEGraphEvent * graphEvent : overlappingEvents)
	{
		_selection.Add(graphEvent);
	}
}

void UEventVisualizationTool::_traceSelection(TSet<USEGraphEvent*>& selection)
{
	std::vector<USEGraphEvent*> events;

	for (USEGraphEvent* e : selection)
		events.push_back(e);

	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	timeline->traceEvents(events);
}












void UAgentPathlineTool::gainFocus()
{
	Super::gainFocus();
}

void UAgentPathlineTool::loseFocus()
{
	Super::loseFocus();
}

void UAgentPathlineTool::oneHandStart(UPrimitiveComponent * hand)
{
	Super::oneHandStart(hand);

	_selection.Empty();
}

void UAgentPathlineTool::oneHandEnd(UPrimitiveComponent * hand)
{
	Super::oneHandEnd(hand);

	if (_selection.Num())
		_pathlinesForSelection(_selection);
	else
	{
		UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

		// show the last 10 seconds of events
		timeline->showAllEvents();

		UVisualization_AgentPathLines * pathlines = _flexSimulation->simulationManager.simulation<UVisualization_AgentPathLines>();

		pathlines->hideTotalHistory();
	}
}

void UAgentPathlineTool::tick(float dt)
{

}

void UAgentPathlineTool::tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	Super::tickOneHand(dt, hand, lastTransform);

	_selectAgent(dt, hand, lastTransform);
}

void UAgentPathlineTool::faceDown_released()
{
	if (_flexSimulation->isPlaying())
		_flexSimulation->pause();
	else
		_flexSimulation->play();
}

void UAgentPathlineTool::faceUp_released(USceneComponent * interactionPoint /* = nullptr */)
{
	physicalInteraction = !physicalInteraction;
}

void UAgentPathlineTool::_selectAgent(float dt, UPrimitiveComponent * hand, FTransform lastTransform)
{
	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	FTransform toLocal = _flexSimulation->owner->GetRootComponent()->GetComponentTransform();

	FVector localHandPosition = toLocal.InverseTransformPosition(hand->GetComponentLocation());

	float brushSqrd = _brushRadius() * _brushRadius();


	for (FGraphNode& node : _flexSimulation->graphSimulation.allNodes)
	{
		if (!node.isValid())
			continue;

		float distSqrd = FVector::DistSquared(localHandPosition, node.position);

		if (distSqrd < brushSqrd)
			_selection.Add(FGraphNodeHandle(node));
	}
}

void UAgentPathlineTool::_pathlinesForSelection(TSet<FGraphNodeHandle>& selection)
{
	std::vector<FGraphNodeHandle> agents;

	for( FGraphNodeHandle& handle : selection)
		agents.push_back(handle);

	UVisualization_AgentPathLines * pathlines = _flexSimulation->simulationManager.simulation<UVisualization_AgentPathLines>();
	UTimelineSimulation * timeline = _flexSimulation->simulationManager.simulation<UTimelineSimulation>();

	int32 curFrame = _flexSimulation->graphSimulation.tickCount;

	pathlines->showTotalHistoryForAgents(agents, 0, curFrame);
}

void UPhysicalInteractionTool::loseFocus()
{
	Super::loseFocus();

	_grabbedNode = FGraphNodeHandle::null;
}

void UPhysicalInteractionTool::faceDown_released()
{
	if (_flexSimulation->isPlaying())
		_flexSimulation->pause();
	else
		_flexSimulation->play();
}

void UPhysicalInteractionTool::oneHandStart(UPrimitiveComponent * handComponent)
{
	Super::oneHandStart(handComponent);

	FVector hand = _flexSimulation->owner->GetTransform().InverseTransformPosition(handComponent->GetComponentLocation());

	if (interactionMode == EPhysicalInteractionType::Grab)
	{
		_cachedCalculatedVelocity.Empty();

		_grabbedNode = FGraphNodeHandle::null;

		FGraph& graph = _flexSimulation->graphSimulation;

		const float interactionDistanceSqrd = std::pow(2.0f, 2.0f);

		float nearestDistanceSqrd = std::numeric_limits<float>::max();
		FGraphNodeHandle nearestHandle;

		for (auto& node : graph.allNodes)
		{
			if (!node.isValid())
				continue;

			float distSqrd = (node.position - hand).SizeSquared();

			if (distSqrd < nearestDistanceSqrd)
			{
				nearestDistanceSqrd = distSqrd;
				nearestHandle = node.handle();
			}
		}

		if (nearestHandle && nearestDistanceSqrd < interactionDistanceSqrd )
		{
			_grabbedNode = nearestHandle;

			FGraphNode& theNode = _grabbedNode(graph);

			_grabbedIdentityMatrix = FRotationMatrix::Make(theNode.orientation.Inverse()) * FTranslationMatrix::Make(-theNode.position);
			_grabbedTransform = FTransform(theNode.orientation, theNode.position, FVector(theNode.scale));

			const FTransform toLocal = _flexSimulation->owner->GetTransform().Inverse();

			_startHand = handComponent->GetComponentTransform() * toLocal;
		}
	}
}

void UPhysicalInteractionTool::oneHandEnd(UPrimitiveComponent * handComponent)
{
	Super::oneHandEnd(handComponent);

	if (!_flexSimulation)
		return;
	
	// put it somewhere far away
	if (interactionMode == EPhysicalInteractionType::Punch)
		_flexSimulation->updateSphereWorldSpace(FVector(10000.0f, 0.0f, 0.0f), 0.1f);
	else if (interactionMode == EPhysicalInteractionType::Grab)
	{
		FGraph& graph = _flexSimulation->graphSimulation;

		auto cachedVelocities = _cachedCalculatedVelocity;
		_flexSimulation->addTickWork([&graph, cachedVelocities] {
			for (auto& pair : cachedVelocities)
			{
				FGraphNode& node = graph.node(pair.Key);

				if (node.hasComponent<FVelocityGraphObject>())
					node.component<FVelocityGraphObject>(graph).linearVelocity = pair.Value;
			}
		});

		_cachedCalculatedVelocity.Empty();
	}
}

void UPhysicalInteractionTool::tickOneHand(float dt, UPrimitiveComponent * handComponent, FTransform lastTransform)
{
	Super::tickOneHand(dt, handComponent, lastTransform);

	const FTransform toLocal = _flexSimulation->owner->GetTransform().Inverse();

	const FTransform lastLocal = lastTransform * toLocal;
	const FTransform local = handComponent->GetComponentTransform() * toLocal;

	const FVector hand = local.GetLocation();
	const FVector lastHand = lastLocal.GetLocation();

	const float invDt = 1.0f / dt;

	if (!_flexSimulation)
		return;

	if (interactionMode == EPhysicalInteractionType::Punch)
		_flexSimulation->updateSphereWorldSpace(handComponent->GetComponentLocation(), _brushRadius());
	else if (interactionMode == EPhysicalInteractionType::Grab)
	{
		if (!_grabbedNode)
			return;

		FGraphNodeHandle grabbedNode = _grabbedNode;

		FGraph& graph = _flexSimulation->graphSimulation;

		FVector velocity = (hand - lastHand) / (dt > 0.0f ? dt : 1.0f);

		FVector offset = (local.GetRotation() * _startHand.GetRotation().Inverse()).RotateVector(_grabbedTransform.GetLocation() - _startHand.GetLocation());

		FVector pr = grabbedNode(graph).position;
		FVector pr_ = hand + offset;

		FQuat rotationalDiff = local.GetRotation() * lastLocal.GetRotation().Inverse();

		auto& cachedVelocity = _cachedCalculatedVelocity;

		_flexSimulation->addTickWork([&graph, &cachedVelocity, grabbedNode, velocity, invDt, pr, pr_, rotationalDiff] {
			FGraphNode& node = graph.node(grabbedNode);

			if (!node.isValid())
				return;

			auto translate = [&](FVector p) -> FVector {
				return rotationalDiff.RotateVector(p - pr) + pr_;
			};

			auto rotate = [&](FQuat rotation) -> FQuat {
				return rotationalDiff * rotation;
			};
			
			if (FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(graph, node.handle()))
			{
				rigidHandle(graph).each<FFlexRigidBodyConnection>(graph, [&](FGraphNodeHandle subRigidHandle, FFlexRigidBodyConnection& rigidConnection) {
					FGraphNode& subNode = subRigidHandle(graph);

					const FVector newPosition = translate(subNode.position);

					// for stability while holding, set linear velocity to the velocity of one node
					subNode.component<FVelocityGraphObject>(graph).linearVelocity = velocity;

					// calculate velocity from the change in position and cache it, we'll set it on oneHandEnd
					cachedVelocity.Add(subRigidHandle, (newPosition - subNode.position) * invDt + velocity);
					subNode.position = newPosition;
					subNode.orientation = rotate(subNode.orientation);
				});
			}
			else
			{
				if (node.hasComponent<FVelocityGraphObject>())
					node.component<FVelocityGraphObject>(graph).linearVelocity = velocity;

				node.position = translate(node.position);
				node.orientation = rotate(node.orientation);
			}
		});

	}
}

bool UPhysicalInteractionTool::shouldShowBrush()
{
	return interactionMode == EPhysicalInteractionType::Punch;
}

