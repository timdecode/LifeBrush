// Copyright (c) 2018 Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Simulation/FlexElements.h"

#include "Visualization/EdgeFactory.h"
#include "Visualization/Timeline.h"

#include "IXRTrackingSystem.h"

void UTimelineSimulation::attach()
{
	_initGlyphIndicesForEventClass();
	_initISMCs();


	_loadAllFrameEvents();

	_currentFrame = NewObject<USEGraphFrame>(this);
}

void UTimelineSimulation::tick(float deltaT)
{
	std::vector<USEGraphEvent*> toRemove;

	// first remove instances
	int32 minFrame = graph->tickCount - maxEventsHistory;
	
	for (USEGraphEvent * event : _loadedEvents)
	{
		if (event->frame < minFrame)
		{
			UInstancedStaticMeshComponent * mesh = event->_transientISMC;

			mesh->RemoveInstance(mesh->GetInstanceCount() - 1);

			toRemove.push_back(event);
		}
	}

	for (USEGraphEvent * event : toRemove)
		_loadedEvents.erase(event);

	// update instances
	std::unordered_map<UInstancedStaticMeshComponent*, int> indexTable;

	for (USEGraphEvent * event : _loadedEvents)
	{
		UInstancedStaticMeshComponent * mesh = event->_transientISMC;

		int& index = indexTable[mesh];

		FTransform transform(event->position);

		transform.SetScale3D(FVector(event->_transientScale));

		mesh->UpdateInstanceTransform(index, transform);

		index++;
	}
}

void UTimelineSimulation::snapshotToActor(AActor * actor)
{
	for (auto& pair : _instancedStaticeMeshes)
	{
		UInstancedStaticMeshComponent * instance = pair.second;

		UInstancedStaticMeshComponent * newInstance = NewObject<UInstancedStaticMeshComponent>(actor);

		newInstance->SetStaticMesh(instance->GetStaticMesh());

		for (int32 mi = 0; mi < instance->GetNumMaterials(); mi++)
			newInstance->SetMaterial(mi, instance->GetMaterial(mi));

		newInstance->SetCollisionEnabled(ECollisionEnabled::NoCollision);
		newInstance->SetCollisionProfileName(TEXT("NoCollision"));
		newInstance->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		actor->AddInstanceComponent(newInstance);

		newInstance->RegisterComponent();

		// copy the instances
		int32 n = instance->GetInstanceCount();
		for (int32 i = 0; i < n; i++)
		{
			FTransform transform;

			instance->GetInstanceTransform(i, transform, true);

			newInstance->AddInstanceWorldSpace(transform);
		}
	}
}

void UTimelineSimulation::preTick(float deltaT)
{
	beginFrame();
}

void UTimelineSimulation::postTick(float deltaT)
{
	endFrame();
}

void UTimelineSimulation::beginFrame()
{
	_didGenerateEventsInFrame = false;

	_currentFrame->number = graph->tickCount;
}

void UTimelineSimulation::endFrame()
{
	if (!_didGenerateEventsInFrame)
		return;

	_currentFrame->snapshot.snapshot(*graph);

	timeline->sparseFrames.Add(_currentFrame->number, _currentFrame);

	_loadFrameEvents(*_currentFrame);

	UVisualization_AgentPathLines * pathlines = simulationManager->simulation<UVisualization_AgentPathLines>();

	pathlines->captureFrame();

	// create a new frame
	_currentFrame = NewObject<USEGraphFrame>(this);
}

void UTimelineSimulation::_loadEvent(USEGraphEvent& graphEvent)
{
	UInstancedStaticMeshComponent& mesh = meshForEvent(graphEvent);

	FTransform transform(graphEvent.position);

	FEventGlyph& glyph = glyphForEventClass(graphEvent.GetClass());

	transform.SetScale3D(FVector(glyph.glyphScale));

	mesh.AddInstance(transform);

	graphEvent._transientISMC = &mesh;
	graphEvent._transientScale = glyph.glyphScale;

	_loadedEvents.insert(&graphEvent);
}

void UTimelineSimulation::_loadFrameEvents(USEGraphFrame& frame)
{
	for (USEGraphEvent * graphEvent : frame.events)
	{
		_loadEvent(*graphEvent);
	}
}

void UTimelineSimulation::_clearVisualizations()
{
	for (auto& pair : _instancedStaticeMeshes)
	{
		pair.second->ClearInstances();
	}
}

void UTimelineSimulation::_loadAllFrameEvents(int32 minFrame)
{
	_clearVisualizations();

	for (auto& pair : timeline->sparseFrames)
	{
		if (pair.Key < minFrame)
			continue;

		USEGraphFrame * frame = pair.Value;

		_loadFrameEvents(*frame);
	}
}

void UTimelineSimulation::_initGlyphIndicesForEventClass()
{
	_glyphIndicesForEventClass.clear();

	TArray<int32> sortedGlyphIndices;

	for (int i = 0; i < glyphPrototypes.Num(); ++i)
		sortedGlyphIndices.Add(i);

	// sort the glyphs, so the deepest ones are first
	sortedGlyphIndices.Sort([this](const int32 a, const int32 b) {
		UClass * aClass = glyphPrototypes[a].eventClass.Get();
		UClass * bClass = glyphPrototypes[b].eventClass.Get();

		return aClass->IsChildOf(bClass);
	});

	for (TObjectIterator<UClass> It; It; ++It)
	{
		if (It->IsChildOf(USEGraphEvent::StaticClass()) && !It->HasAnyClassFlags(CLASS_Abstract))
		{
			UClass * eventClass = *It;

			_glyphIndicesForEventClass[eventClass] = -1;

			// now find the nearest glyph (of which we are a subclass)
			for (int gi : sortedGlyphIndices)
			{
				FEventGlyph& glyph = glyphPrototypes[gi];

				if (eventClass->IsChildOf(glyph.eventClass))
				{
					_glyphIndicesForEventClass[eventClass] = gi;
					break;
				}
			}
		}
	}
}

void UTimelineSimulation::_initISMCs()
{
	// build our mesh keys
	std::vector<FGraphMeshKey> meshKeys;

	for (FEventGlyph& glyph : glyphPrototypes)
	{
		meshKeys.emplace_back();

		FGraphMeshKey& key = meshKeys.back();

		key.staticMesh = glyph.glyphMesh;
		key.material = glyph.glyphMaterial;
	}

	// add a default key too
	{
		meshKeys.emplace_back();

		FGraphMeshKey& key = meshKeys.back();

		key.staticMesh = defaultGlyph.glyphMesh;
		key.material = defaultGlyph.glyphMaterial;
	}

	// create one ISMC for each key (duplicate keys will map to the same ISMC)
	std::unordered_map<FGraphMeshKey, UInstancedStaticMeshComponent*> ismcs;

	for (FGraphMeshKey& key : meshKeys)
	{
		auto found = ismcs.find(key);

		if (found == ismcs.end())
		{
			UInstancedStaticMeshComponent * mesh = NewObject<UInstancedStaticMeshComponent>(actor);
			mesh->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

			mesh->SetStaticMesh(key.staticMesh);
			mesh->SetMaterial(0, key.material);

			mesh->RegisterComponent();

			ismcs[key] = mesh;
		}
	}

	// finally, set a ISMC for each USEGraphEvent class (we can have duplicates)
	for (auto& pair : _glyphIndicesForEventClass)
	{
		UClass * eventClass = pair.first;

		FEventGlyph& glyph = pair.second >= 0 ? glyphPrototypes[pair.second] : defaultGlyph;

		FGraphMeshKey key;
		key.staticMesh = glyph.glyphMesh;
		key.material = glyph.glyphMaterial;

		UInstancedStaticMeshComponent * mesh = ismcs[key];

		_instancedStaticeMeshes[eventClass] = mesh;
	}
}

FEventGlyph& UTimelineSimulation::glyphForEventClass(TSubclassOf<USEGraphEvent> eventClass)
{
	int32 index = _glyphIndicesForEventClass[eventClass.Get()];

	if (index < 0)
		return defaultGlyph;
	else
		return glyphPrototypes[index];
}

UInstancedStaticMeshComponent& UTimelineSimulation::meshForEvent(USEGraphEvent& graphEvent)
{
	return *_instancedStaticeMeshes[graphEvent.GetClass()];
}

std::vector<USEGraphEvent*> UTimelineSimulation::eventsOverlappingPosition(const FVector position, float radius)
{
	std::vector<USEGraphEvent*> overlaps;

	const float radiusSqrd = radius * radius;

	/*for (auto& pair : timeline->sparseFrames)
	{
		USEGraphFrame * graphFrame = pair.Value;

		for (USEGraphEvent * graphEvent : graphFrame->events)
		{
			if (FVector::DistSquared(graphEvent->position, position) < radiusSqrd)
			{
				overlaps.push_back(graphEvent);
			}
		}
	}*/

	for (USEGraphEvent * event : _loadedEvents)
	{
		if (FVector::DistSquared(event->position, position) < radiusSqrd)
		{
			overlaps.push_back(event);
		}
	}

	return overlaps;
}

void UTimelineSimulation::traceEvents(std::vector<USEGraphEvent*> events)
{

	std::vector<int32> sortedFrameNumbers;

	for (auto& pair : timeline->sparseFrames)
	{
		sortedFrameNumbers.push_back(pair.Key);
	}

	std::sort(sortedFrameNumbers.begin(), sortedFrameNumbers.end());

	// trace the events
	TSet<FGraphNodeHandle> allTracedAgents;
	TSet<USEGraphEvent*> allTracedEvents;

	for (USEGraphEvent * graphEvent : events)
	{
		TSet<FGraphNodeHandle> tracedAgents;
		TSet<USEGraphEvent*> tracedEvents;

		_expandFrames2(graphEvent, sortedFrameNumbers, tracedAgents, tracedEvents, 1);

		allTracedAgents.Append(tracedAgents);
		allTracedEvents.Append(tracedEvents);
	}

	// show them (but hide the untraced events)
	_clearVisualizations();

	for (USEGraphEvent * graphEvent : allTracedEvents)
	{
		_loadEvent(*graphEvent);
	}

	// show pathlines
	UVisualization_AgentPathLines * pathlines = simulationManager->simulation<UVisualization_AgentPathLines>();

	_cached_showAgentPathLines = pathlines->showAgentPathLines;

	
	std::vector<FGraphNodeHandle> agentsAsVector;
	for (auto& a : allTracedAgents)
		agentsAsVector.push_back(a);

	int32 minFrame = 0;
	int32 maxFrame = 0;

	for (USEGraphEvent * event : events)
	{
		if (event->frame > maxFrame)
			maxFrame = event->frame;

		if (event->frame < minFrame)
			minFrame = event->frame;
	}

	pathlines->showTotalHistoryForAgents(agentsAsVector, minFrame - maxPathHistory, maxFrame + maxPathHistory);
}

void UTimelineSimulation::_showAgentPathLines(TSet<FGraphNodeHandle>& agentSet, TSet<USEGraphEvent*>& eventSet)
{
	// we need to sort the events by frame

}

void UTimelineSimulation::_expandFrame(USEGraphFrame * frame, TSet<FGraphNodeHandle> &agents_out, TSet<USEGraphEvent *> &includedEvents_out)
{
	// find all agents in each frame that touched this event
	for (USEGraphEvent * frameEvent : frame->events)
	{
		bool shouldAdd = false;

		for (FGraphNodeHandle handle : frameEvent->otherAgents)
		{
			if (agents_out.Contains(handle))
			{
				shouldAdd = true;
				break;
			}
		}

		if (!shouldAdd && agents_out.Contains(frameEvent->triggeringAgent))
		{
			shouldAdd = true;
		}

		if (shouldAdd)
		{
			agents_out.Append(frameEvent->otherAgents);
			agents_out.Add(frameEvent->triggeringAgent);

			includedEvents_out.Add(frameEvent);
		}
	}
}

bool UTimelineSimulation::_parseFrame(USEGraphFrame * frame, TSet<FGraphNodeHandle>& agents_out, TSet<USEGraphEvent*>& events_out)
{
	bool didAdd = false;

	for (USEGraphEvent * event : frame->events)
	{
		bool overlap = false;

		if (agents_out.Contains(event->triggeringAgent))
			overlap = true;

		for (FGraphNodeHandle& agent : event->otherAgents)
		{
			if (agents_out.Contains(agent))
			{
				overlap = true;
				break;
			}
		}

		if (overlap)
		{
			didAdd = true;

			agents_out.Add(event->triggeringAgent);
			agents_out.Append(event->otherAgents);
		}
	}

	return didAdd;
}

void UTimelineSimulation::_expandFrames2(
	USEGraphEvent * graphEvent, 
	const std::vector<int32>& sortedFrames, 
	TSet<FGraphNodeHandle> &agents_out, 
	TSet<USEGraphEvent *> &includedEvents_out,
	int32 maxEventHops
)
{
	includedEvents_out.Add(graphEvent);

	// find our start frame
	int32 start = 0;
	for (start = 0; start < sortedFrames.size(); ++start)
	{
		if (sortedFrames[start] == graphEvent->frame)
			break;
	}

	if (start >= sortedFrames.size())
		return;

	// backwards loop (including start)
	TSet<FGraphNodeHandle> backwardsAgents;

	backwardsAgents.Add(graphEvent->triggeringAgent);
	backwardsAgents.Append(graphEvent->otherAgents);

	size_t didAddCount = 0;

	for (int i = start; i >= 0; --i)
	{
		USEGraphFrame * frame = timeline->sparseFrames[sortedFrames[i]];

		if (_parseFrame(frame, backwardsAgents, includedEvents_out) )
		{
			didAddCount++;

			if (didAddCount >= maxEventHops)
				break;
		}
	}

	// forward loop (past start)
	TSet<FGraphNodeHandle> forwardsAgents;

	forwardsAgents.Add(graphEvent->triggeringAgent);
	forwardsAgents.Append(graphEvent->otherAgents);

	for (int i = start + 1; i < sortedFrames.size(); ++i)
	{
		USEGraphFrame * frame = timeline->sparseFrames[sortedFrames[i]];

		if (_parseFrame(frame, backwardsAgents, includedEvents_out))
		{
			didAddCount++;

			if (didAddCount >= maxEventHops)
				break;
		}
	}

	agents_out.Append(backwardsAgents);
	agents_out.Append(forwardsAgents);
}

void UTimelineSimulation::_expandFrames(
	USEGraphEvent * graphEvent,
	const std::vector<int32>& sortedFrames,
	TSet<FGraphNodeHandle> &agents_out,
	TSet<USEGraphEvent *> &includedEvents_out,
	int32 maxFrameDistance
)
{
	includedEvents_out.Add(graphEvent);

	// find our start frame
	int32 start = 0;
	for (start = 0; start < sortedFrames.size(); ++start)
	{
		if (sortedFrames[start] == graphEvent->frame)
			break;
	}

	if (start >= sortedFrames.size())
		return;

	// backwards loop (including start)
	TSet<FGraphNodeHandle> backwardsAgents;

	backwardsAgents.Add(graphEvent->triggeringAgent);
	backwardsAgents.Append(graphEvent->otherAgents);

	for (int i = start; i >= 0; --i)
	{
		USEGraphFrame * frame = timeline->sparseFrames[sortedFrames[i]];

		if (frame->number < graphEvent->frame - maxFrameDistance)
			break;

		_expandFrame(frame, backwardsAgents, includedEvents_out);
	}

	// forward loop (past start)
	TSet<FGraphNodeHandle> forwardsAgents;

	backwardsAgents.Add(graphEvent->triggeringAgent);
	backwardsAgents.Append(graphEvent->otherAgents);

	for (int i = start + 1; i < sortedFrames.size(); ++i)
	{
		USEGraphFrame * frame = timeline->sparseFrames[sortedFrames[i]];

		if (frame->number > graphEvent->frame + maxFrameDistance)
			break;

		_expandFrame(frame, forwardsAgents, includedEvents_out);
	}

	agents_out.Append(backwardsAgents);
	agents_out.Append(forwardsAgents);
}

void UTimelineSimulation::clear()
{
	_clearVisualizations();

	timeline->sparseFrames.Empty();
}

void UTimelineSimulation::showAllEvents()
{	
	int32 minFrame = graph->tickCount - maxEventsHistory;

	_loadAllFrameEvents(minFrame);
}




void UVisualization_AgentPathLines::attach()
{
	_initPipeFactory();

	_timer = secondsBetweenSegments;
}

void UVisualization_AgentPathLines::clear()
{
	_clearDynamicEdgeFactory();

	totalHistory.Empty();

	lastPositionCache.Empty();
}

void UVisualization_AgentPathLines::_initPipeFactory()
{
	if (!_dynamicEdgeFactory && actor)
	{
		_dynamicEdgeFactory = NewObject<UBlockyEdgeFactory>(actor, TEXT("pipe_edge_factory"));

		if (defaultLineMaterial)
			_dynamicEdgeFactory->material = defaultLineMaterial;
	
	}

	if (!_agentPathFactory && actor)
	{
		_agentPathFactory = NewObject<UColoredLineFactory>(actor, TEXT("line_factory"));

		if( defaultLineMaterial )
			_agentPathFactory->material = defaultLineMaterial;

	}

	_agentPathFactory->RegisterComponent();
	_dynamicEdgeFactory->RegisterComponent();
}

void UVisualization_AgentPathLines::_initCache()
{
	if (lastPositionCache.Num() <= graph->allNodes.Num())
	{
		lastPositionCache.SetNum(graph->allNodes.Num());
	}

	int32 i = 0;

	for (FGraphNode& node : graph->allNodes)
	{
		if (!node.isValid())
			continue;

		lastPositionCache[node.id] = node.position;
	}
}


void UVisualization_AgentPathLines::tick(float deltaT)
{
	_tickVisibility(deltaT);

	if (_timer >= secondsBetweenSegments)
	{
		captureFrame();

		_cachePositions();

		_timer = 0.0f;
	}
	else
	{
		_timer += deltaT;
	}

	if (showAgentPathLines && _dynamicEdgeFactory)
		_updateDynamicPaths();
}

void UVisualization_AgentPathLines::tick_paused(float deltaT)
{
	_tickVisibility(deltaT);
}

void UVisualization_AgentPathLines::snapshotToActor(AActor * actor)
{
	URuntimeMeshComponent * rmc = _agentPathFactory->_runtimeMeshComponent;

	if (rmc->GetNumSections() == 0)
		return;

	URuntimeMeshComponent * newRMC = NewObject<URuntimeMeshComponent>(actor);

	newRMC->SetShouldSerializeMeshData(true);

	int32 newRMCSection = 0;

	auto sectionIds = rmc->GetRuntimeMesh()->GetSectionIds();

	for (auto& pair : _materialToSection)
	{
		int32 section = pair.second;
		UMaterialInterface * material = pair.first;

		if (!rmc->DoesSectionExist(section))
			continue;

		auto readonlySection = rmc->GetSectionReadonly(section);

		if (!readonlySection || readonlySection->NumIndices() == 0)
			continue;

		auto builder = MakeRuntimeMeshBuilder(*readonlySection);

		readonlySection->CopyTo(builder, true);

		newRMC->CreateMeshSection(newRMCSection, builder, false, EUpdateFrequency::Infrequent);
		newRMC->SetSectionMaterial(newRMCSection, material);
		newRMC->SetMaterial(newRMCSection, material);

		newRMCSection++;
	}


	newRMC->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

	actor->AddInstanceComponent(newRMC);

	newRMC->RegisterComponent();


}

void UVisualization_AgentPathLines::_updateDynamicPaths()
{
	TArray<UEdgeFactory::EdgeParameter> edges;

	edges.SetNumUninitialized(graph->numNodes());

	int32 i = 0;

	for (FGraphNode& node : graph->allNodes)
	{
		if (!node.isValid())
			continue;

		UEdgeFactory::EdgeParameter& edge = edges[i++];

		edge.a = lastPositionCache[node.id];
		edge.b = node.position;
		edge.halfExtents = pathRadius;

	}

	_dynamicEdgeFactory->processEdges(edges, FVector::UpVector, historyIndex, defaultLineMaterial);

	historyIndex = (historyIndex + 1) % historySize;
}

void UVisualization_AgentPathLines::_tickVisibility(float deltaT)
{
	if (!camera)
		return;

	//FQuat cameraOrientation;
	//FVector cameraPosition;
	//GEngine->XRSystem->GetCurrentPose(0, cameraOrientation, cameraPosition);

	//APawn * playerPawn = UGameplayStatics::GetPlayerPawn(GetWorld(), 0);

	//cameraPosition = playerPawn->GetTransform().TransformPosition(cameraPosition);

	FVector cameraPosition = camera->GetComponentLocation();

	cameraPosition = actor->GetTransform().InverseTransformPosition(cameraPosition);

	auto& meshes = graph->componentStorage<FGraphMesh>();

	for (FGraphMesh& mesh : meshes)
	{
		FGraphNodeHandle handle(mesh.nodeIndex);

		auto found = _cachedVisibility.find(handle);
		if (found == _cachedVisibility.end())
			_cachedVisibility[handle] = mesh.visible;

		bool& cachcedVisibility = found->second;

		// don't hide the agent, if it was one of the total history agents
		if (_totalHistoryAgents.find(handle) != _totalHistoryAgents.end())
		{
			mesh.visible = cachcedVisibility;
			continue;
		}

		FGraphNode& node = graph->node(mesh.nodeIndex);

		FVector dir = (node.position - cameraPosition);

		bvh::Ray ray;
		ray.origin = node.position;
		ray.direction = dir.GetSafeNormal();

		bvh::IntersectionInfo intersection;

		bool hit = pathBVH.getIntersection(ray, intersection, true, true);

		mesh.visible = !hit;
	}
}


void UVisualization_AgentPathLines::buildBVH(ColoredQuadFactory& quadFactory)
{
	std::vector<bvh::Triangle> triangles;

	for (int32 i = 0; i < quadFactory.indices.Num(); i += 3)
	{
		triangles.emplace_back();

		bvh::Triangle& triangle = triangles.back();

		triangle.vertices[0] = quadFactory.vertices[i + 0];
		triangle.vertices[1] = quadFactory.vertices[i + 1];
		triangle.vertices[2] = quadFactory.vertices[i + 2];
	}

	pathBVH.build(triangles);
}

int32 UVisualization_AgentPathLines::_sectionForMaterial(UMaterialInterface* material)
{
	auto found = _materialToSection.find(material);

	if (found == _materialToSection.end())
	{
		int32 section = _materialToSection.size();

		_materialToSection[material] = section;

		return section;
	}
	else
		return found->second;
}

FLinearColor UVisualization_AgentPathLines::_colorForNode(FGraphNodeHandle node)
{
	return FColor::Red;
}

std::unordered_map<UMaterialInterface*, FColoredLineBuilder>
UVisualization_AgentPathLines::coloredLinesForAgents(const std::vector<FGraphNodeHandle>& agents, float radius, int32 minFrame, int32 maxFrame)
{
	// map from frames with agent positions to an array of positions for each agent
	std::unordered_map<FGraphNodeHandle, FColoredLineBuilder> lineBuilders;

	FPositionSnapshot * last = nullptr;

	for (auto& current : totalHistory)
	{
		if (current.frameNumber > maxFrame)
			break;

		if (last && current.frameNumber > minFrame)
		{
			for (FGraphNodeHandle agent : agents)
			{
				if (agent.index >= current.agentPositions.Num() || agent.index >= last->agentPositions.Num())
					continue;

				auto& builder = lineBuilders[agent];

				FVector p = current.agentPositions[agent.index];

				if (builder.lineSegments.Num() == 0)
					builder.begin(p, radius);
				else
					builder.addPoint(p, radius);
			}
		}

		last = &current;
	}

	// find the longest
	int32 longest = 0;
	for (auto& pair : lineBuilders)
	{
		FGraphNodeHandle node = pair.first;
		auto& builder = pair.second;

		if (builder.lineSegments.Num() > longest)
			longest = builder.lineSegments.Num();
	}

	// cap them with an end
	for (auto& pair : lineBuilders)
	{
		FGraphNodeHandle node = pair.first;
		auto& builder = pair.second;

		// one one element, we need to nuke it
		auto n = builder.lineSegments.Num();
		if (n == 1)
			builder.clear();
		else if (n > 1)
		{
			builder.lineSegments.Last().type = FColoredLineBuilder::LineElement::SegmentType::End;
		}

		// color the line segments
		const FLinearColor baseColor = _colorForNode(node);

		int i = 0;
		for (auto& element : builder.lineSegments)
		{
			float scale = float(i) / float(longest);


			element.color = (baseColor * scale).ToFColor(false);

			i++;
		}
	}
	
	// Merge builders into one builder per material
	std::unordered_map<UMaterialInterface*, FColoredLineBuilder> builderByMaterial;
	{
		// allocate
		std::unordered_map<UMaterialInterface*, int32> nums;

		for (auto& pair : lineBuilders)
		{
			FGraphNodeHandle node = pair.first;
			auto& builder = pair.second;

			if (!node(graph).hasComponent<FGraphMesh>())
				continue;

			UMaterialInterface * material = node(graph).component<FGraphMesh>(graph).material;

			auto& n = nums[material];

			n += builder.lineSegments.Num();
		}

		for (auto& pair : nums)
		{
			FColoredLineBuilder& mergedBuilder = builderByMaterial[pair.first];

			mergedBuilder.lineSegments.Reserve(pair.second);
		}

		// append
		for (auto& pair : lineBuilders)
		{
			FGraphNodeHandle node = pair.first;
			auto& builder = pair.second;

			if (!node(graph).hasComponent<FGraphMesh>())
				continue;

			UMaterialInterface * material = node(graph).component<FGraphMesh>(graph).material;

			auto& mergedBuilder = builderByMaterial[material];

			mergedBuilder.lineSegments.Append(builder.lineSegments);
		}
	}

	return builderByMaterial;
}


void UVisualization_AgentPathLines::showTotalHistoryForAgents(const std::vector<FGraphNodeHandle>& agents, int32 minFrame, int32 maxFrame)
{
	showAgentPathLines = false;

	_clearDynamicEdgeFactory();

	historyIndex = 0;


	auto coloredBuilders = coloredLinesForAgents(agents, pathRadius, minFrame, maxFrame);

	for (auto& pair : coloredBuilders)
	{
		UMaterialInterface * material = pair.first;
		auto& builder = pair.second;

		int32 section = _sectionForMaterial(material);

		_agentPathFactory->commitSection(builder, section, material);
	}

	// refactor coloredBuilders into fatEdges
	ColoredQuadFactory fatEdgeFactory;

	for (auto& pair : coloredBuilders)
	{
		UMaterialInterface * material = pair.first;
		auto& builder = pair.second;

		// grow the radius
		for (auto& e : builder.lineSegments)
			e.radius += 2.0f;

		builder.appendToQuadFactory(fatEdgeFactory, _agentPathFactory->uvBottomY, _agentPathFactory->uvTopY, _agentPathFactory->uvXScale);
	}

	buildBVH(fatEdgeFactory);

	_totalHistoryAgents.insert(agents.begin(), agents.end());
}

void UVisualization_AgentPathLines::_clearDynamicEdgeFactory()
{
	for (int i = 0; i < _dynamicEdgeFactory->numSections(); ++i)
	{
		_dynamicEdgeFactory->clearSection(i);
	}
}

void UVisualization_AgentPathLines::_clearAgentPathFactory()
{
	for (int i = 0; i < _agentPathFactory->numSections(); ++i)
	{
		_agentPathFactory->clearSection(i);
	}
}

void UVisualization_AgentPathLines::hideTotalHistory()
{
	_clearAgentPathFactory();

	_totalHistoryAgents.clear();

	historyIndex = 0;

	pathBVH.clear();

	auto& meshes = graph->componentStorage<FGraphMesh>();

	for (FGraphMesh& mesh : meshes)
	{
		mesh.visible = _cachedVisibility[FGraphNodeHandle(mesh.nodeIndex)];
	}
}

void UVisualization_AgentPathLines::_cachePositions()
{
	if (lastPositionCache.Num() < graph->allNodes.Num())
		_initCache();

	for (FGraphNode& node : graph->allNodes)
	{
		if (!node.isValid())
			continue;

		lastPositionCache[node.id] = node.position;
	}
}

void UVisualization_AgentPathLines::captureFrame()
{
	_cachePositions();

	totalHistory.Emplace();
	FPositionSnapshot& snapshot = totalHistory.Last();

	snapshot.agentPositions.SetNum(graph->allNodes.Num());
	for (FGraphNode& node : graph->allNodes)
	{
		if (!node.isValid())
			continue;

		snapshot.agentPositions[node.id] = node.position;
	}

	snapshot.frameNumber = graph->tickCount;
}

void UVisualization_AgentPathLines::nodeAdded(FGraphNode& node)
{
	if (lastPositionCache.Num() <= node.id)
	{
		lastPositionCache.SetNumZeroed(node.id + 1);
	}

	lastPositionCache[node.id] = node.position;

	if (node.hasComponent<FGraphMesh>())
	{
		_cachedVisibility[FGraphNodeHandle(node)] = node.component<FGraphMesh>(graph).visible;
	}
}

void UVisualization_AgentPathLines::nodeRemoved(FGraphNode& oldNode)
{
	// don't leak memory if we have less nodes now
	if (graph->numNodes() < lastPositionCache.Num())
	{
		lastPositionCache.SetNum(graph->numNodes());
	}


	if (_cachedVisibility.find(FGraphNodeHandle(oldNode)) != _cachedVisibility.end())
		_cachedVisibility.erase(FGraphNodeHandle(oldNode));
}
