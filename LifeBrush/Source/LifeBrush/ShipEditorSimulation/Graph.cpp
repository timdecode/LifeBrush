// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"
#include "Graph.h"

void FGraph::init()
{
	if (_didInit)
		return;

	// allocate our component storage for all loaded component classes
	_initComponentStorage();

	_initEdgeStorage();

	_assignComponents();

	_didInit = true;
}

void FGraph::_assignComponents()
{
	// update the graph node component arrays
	for (auto& graphNode : allNodes)
		graphNode.components.Empty();

	for (auto& storage : _componentStorage)
	{
		ComponentType type = FGraphObject::componentType(storage.componentClass);

		for (int i = 0; i < storage.size(); ++i)
		{
			FGraphObject * object = storage.at(i, storage.componentClass);

			if( !object->isValid() )
				continue;

			FGraphNode& graphNode = node(object->nodeIndex);

			graphNode.components.AddUnique(type);
		}
	}
}

void FGraph::_edgeObjectAdded(FGraphEdgeHandle handle, const EdgeObjectType type)
{
	GraphTransaction * backTransaction = _backTransaction();

	// transactional dispatch
	if (backTransaction)
		backTransaction->addedEdgeObjects[type].emplace_back(handle);

	if (!hasActiveTransaction())
	{
		TransactionContext& transactionContext = currentTransactionContext();

		auto& edgeObjectListeners = transactionContext._edgeObjectListeners;

		auto found = edgeObjectListeners.find(type);

		if (found == edgeObjectListeners.end())
			return;

		for (auto listener : found->second)
			listener->edgeObjectAdded(handle, type);
	}
}

void FGraph::_edgeObjectRemoved(FGraphEdgeHandle handle, EdgeObjectType type)
{
	GraphTransaction * backTransaction = _backTransaction();

	// transactional dispatch
	if (backTransaction)
		backTransaction->removedEdgeObjects[type].emplace_back(handle);

	if (!hasActiveTransaction())
	{
		TransactionContext& transactionContext = currentTransactionContext();

		auto& edgeObjectListeners = transactionContext._edgeObjectListeners;

		auto found = edgeObjectListeners.find(type);

		if (found == edgeObjectListeners.end())
			return;

		for (auto listener : found->second)
			listener->edgeObjectRemoved(handle, type);
	}
}

void FGraph::clear()
{
	beginTransaction();

	for(FGraphNode& node : allNodes)
	{
		if( node.isValid() )
			removeNode( node.id );
	}

	endTransaction();

	allNodes.Empty();
	recycledNodes.Empty();
}

void FGraph::_initComponentStorage()
{
	// get rid of null component storages
	_componentStorage.RemoveAll([](const FComponentStorage& storage) {
		return storage.componentClass == nullptr;
	});

	TArray<UScriptStruct*> derivedComponents;
	ShipEditorUtility::derivedStructs(FGraphObject::StaticStruct(), derivedComponents);
	
	// add any missing storage
	for (UScriptStruct * scriptStruct : derivedComponents)
	{
		bool shouldAdd = true;

		for (FComponentStorage& storage : _componentStorage)
		{
			if (storage.componentClass == scriptStruct)
			{
				shouldAdd = false;
				break;
			}
		}

		if (shouldAdd)
		{
			_componentStorage.Emplace();
			FComponentStorage& storage = _componentStorage.Last();

			storage.componentClass = scriptStruct;
		}
	}

	// reorder _componentStorage so that it is indexed by ComponentType
	_componentStorage.Sort([](const FComponentStorage& a, const FComponentStorage& b) {
		return FGraphObject::componentType(a.componentClass) < FGraphObject::componentType(b.componentClass);
	});
}

void FGraph::_initEdgeStorage()
{
	// get rid of null component storages
	_edgeStorage.RemoveAll([](const FEdgeStorage& storage) {
		return storage._scriptStruct == nullptr;
	});

	TArray<UScriptStruct*> derivedEdgeScriptStructs;
	ShipEditorUtility::derivedStructs(FGraphEdgeObject::StaticStruct(), derivedEdgeScriptStructs);

	// add any missing storage
	for (UScriptStruct * scriptStruct : derivedEdgeScriptStructs)
	{
		bool shouldAdd = true;

		for (FEdgeStorage& storage : _edgeStorage)
		{
			if (storage._scriptStruct == scriptStruct)
			{
				shouldAdd = false;
				break;
			}
		}

		if (shouldAdd)
		{
			_edgeStorage.Emplace();
			FEdgeStorage& storage = _edgeStorage.Last();

			storage._scriptStruct = scriptStruct;
		}
	}

	// reorder _edgeStorage so that it is indexed by EdgeType
	_edgeStorage.Sort([](const FEdgeStorage& a, const FEdgeStorage& b) {
		return FGraphEdgeObject::edgeType(a._scriptStruct) < FGraphEdgeObject::edgeType(b._scriptStruct);
	});

	for (auto& storage : _edgeStorage)
	{
		storage._structSize = storage._scriptStruct->GetStructureSize();
	}
}

FComponentStorage& FGraph::componentStorage( ComponentType type )
{
	checkfSlow(_didInit, TEXT("The graph was accessed without calling init()"));

	FComponentStorage& storage = _componentStorage[type];

	return storage;
}

FEdgeStorage& FGraph::rawEdgeStorage(EdgeObjectType type)
{
	return _edgeStorage[type];
}

FGraphNodeHandle FGraph::anyNode()
{
	for(FGraphNode& node : allNodes)
	{
		if(!node.isValid())
			continue;

		return FGraphNodeHandle( node.id );
	}

	return FGraphNodeHandle();
}

int32 FGraph::addNode(FVector position, FQuat orientation, float scale )
{
	int32 nodeIndex = RecyclingArray::emplace( allNodes, recycledNodes, position, orientation, scale );

	FGraphNode& node = allNodes[nodeIndex];
	node.id = nodeIndex;

	TransactionContext& transactionContext = currentTransactionContext();

	GraphTransaction * backTransaction = _backTransaction();

	if (backTransaction)
	{
		backTransaction->addedNodes.push_back(nodeIndex);
	}

	if( !hasActiveTransaction() )
	{
		TransactionContext& transactionContext = currentTransactionContext();

		for (auto& listener : transactionContext._nodeListeners)
			listener->nodeAdded(node);
	}

	return nodeIndex;
}

void FGraph::removeNode( int32 index )
{
	if(!allNodes[index].isValid())
		return;

	FGraphNode& node = allNodes[index];

	allNodes[index]._isValid = false;

	GraphTransaction * backTransaction = _backTransaction();

	TransactionContext& transactionContext = currentTransactionContext();

	// remove connections first
	for (auto edgeHandle : node.edges)
	{
		removeConnection(FGraphEdgeHandle(edgeHandle));
	}

	if (backTransaction)
	{
		// invalidate the components
		for (auto componentType : node.components)
		{
			auto object = component(node.handle(), componentType);

			if (!object->isValid())
				continue;

			object->invalidate();

			backTransaction->removedComponents[componentType].push_back(node.handle());
		}

		backTransaction->removedNodes.push_back(node.id);
	}

	if( !hasActiveTransaction())
	{
		TransactionContext& transactionContext = currentTransactionContext();

		auto& componentListeners = transactionContext._componentListeners;
		auto& nodeListeners = transactionContext._nodeListeners;

		// dispatch removed components
		for (ComponentType type : node.components)
		{
			FGraphObject * graphObject = component(node.handle(), type);

			graphObject->invalidate();

			auto& listeners = componentListeners[type];

			for (auto& listener : listeners)
			{
				listener->componentRemoved(node.handle(), type);
			}
		}

		// dispatch the removed node
		for (auto& listener : nodeListeners)
		{
			listener->nodeRemoved(node);
		}
	}

}

void FGraph::markNodeDirty( int32 nodeIndex )
{

}

FGraphEdgeHandle FGraph::connectNodes(FGraphNodeHandle a, FGraphNodeHandle b)
{
	int32 index = RecyclingArray::emplace( privateEdges, recycledEdges, a.index, b.index );

	FGraphNode& nodeA = node(a);
	FGraphNode& nodeB = node(b);

	nodeA.edges.Add( index );
	nodeB.edges.Add( index );

	GraphTransaction * backTransaction = _backTransaction();

	if (backTransaction)
	{
		backTransaction->addedConnections.push_back(index);
	}

	if( !hasActiveTransaction() )
		_connectionAdded( index );

	return FGraphEdgeHandle(index);
}

void FGraph::removeConnection(FGraphEdgeHandle handle)
{
    if(!privateEdges[handle.index].isValid())
        return;

	FGraphEdge& edge = privateEdges[handle.index];
	edge.invalidate();

	{
		EdgeObjectType type = 0;
		for (auto& storage : _edgeStorage)
		{
			if (!storage.contains(handle))
				continue;

			removeEdgeObject(handle, type);

			type++;
		}
	}

	GraphTransaction * backTransaction = _backTransaction();

	if (backTransaction)
	{
		backTransaction->removedConnections.push_back(handle.index);
	}

	if (!hasActiveTransaction())
		_connectionRemoved(handle.index);
}

TArray<FGraphEdgeHandle> FGraph::edgesBetweenNodes(TArray<FGraphNodeHandle>& nodes)
{
	auto set = TSet<FGraphNodeHandle>(nodes);

	TSet<FGraphEdgeHandle> edgeSet;

	for (auto node : nodes)
	{
		for (auto ei : node(this).edges)
		{
			FGraphEdgeHandle edgeHandle = FGraphEdgeHandle(ei);

			FGraphEdge& theEdge = edge(edgeHandle);

			if (set.Contains(FGraphNodeHandle(theEdge.a)) && set.Contains(FGraphNodeHandle(theEdge.b)))
			{
				edgeSet.Add(edgeHandle);
			}
		}
	}

	return edgeSet.Array();
}

FGraph& FGraph::operator=( const FGraph& other )
{
	allNodes = other.allNodes;
	recycledNodes = other.recycledNodes;
	privateEdges = other.privateEdges;
	recycledEdges = other.recycledEdges;

	_componentStorage = other._componentStorage;

	return *this;
}


FGraphNode& FGraphObject::node(FGraph& graph)
{
	return graph.node(nodeIndex);
}

// ------------------------------------------------------------
// FGraphObject
// ------------------------------------------------------------

UScriptStruct* FGraphObject::componentStruct( ComponentType type )
{
	return StructTypeManager<FGraphObject>::structForType( type );
}

ComponentType FGraphObject::componentType( UScriptStruct* typeStruct )
{
	return StructTypeManager<FGraphObject>::typeForStruct( typeStruct );
}


UScriptStruct* FGraphEdgeObject::edgeStruct(EdgeObjectType type)
{
	return StructTypeManager<FGraphEdgeObject>::structForType(type);
}

EdgeObjectType FGraphEdgeObject::edgeType(UScriptStruct* typeStruct)
{
	return StructTypeManager<FGraphEdgeObject>::typeForStruct(typeStruct);
}

FGraphNodeHandle::FGraphNodeHandle( FGraphNode& node )
{
	index = node.id;
}

FGraphNode& FGraphNodeHandle::node(FGraph& graph)
{
	return graph.node(*this);
}

FGraphNode& FGraphNodeHandle::operator()( FGraph& graph )
{
	return graph.allNodes[index];
}

FGraphNode& FGraphNodeHandle::operator()( FGraph* graph )
{
	return graph->allNodes[index];
}

const FGraphNodeHandle FGraphNodeHandle::null(-1);

FGraphNodeHandle FGraphEdge::other( FGraphNodeHandle node )
{
	bool isA = node.index == a;
	bool isB = node.index == b;

	assert( isA || isB );

	return FGraphNodeHandle( isA ? b : a );
}

FGraphEdge& FGraphEdgeHandle::operator()(FGraph& graph)
{
	return graph.privateEdges[index];
}

FGraphEdge& FGraphEdgeHandle::operator()(FGraph* graph)
{
	return graph->privateEdges[index];
}

const FGraphEdgeHandle FGraphEdgeHandle::null(-1);


bool FComponentStorage::Serialize(FArchive& Ar)
{
	Ar.UsingCustomVersion(FGraphVersion::GUID);

	Ar << componentClass;
	Ar << _size;

	const int64 location_endOfStructs = Ar.Tell();
	int64 endOfStructs = 0;

	if (Ar.IsLoading() || Ar.IsSaving())
	{
		Ar << endOfStructs;
	}

	// The componentClass is gone, we have to seek past the meaningless data.
	if (componentClass == nullptr && Ar.IsLoading())
	{
		Ar.Seek(endOfStructs);

		pool.resize(_size, 0);

		_size = 0;
		return true;
	}

	if(componentClass)
		pool.resize(_size, componentClass->GetStructureSize());

	// Serialize each component
	for (size_t i = 0; i < _size; ++i)
	{
		FGraphObject * object = at(i, componentClass);

		// very important to call the C++ constructor
		if (Ar.IsLoading())
			componentClass->GetCppStructOps()->Construct(object);

		componentClass->SerializeItem(Ar, object, nullptr);

		_objectIndexToElement.emplace(object->nodeIndex, i);
	}

	// If we are saving, we have written our component data. Now, we know where our data ends. We need
	// this information for the future possibility where we don't have the component class anymore, so
	// that we can skip over the data.
	if(Ar.IsSaving())
	{
		// write where we finished writing, so we can skip to here on loading if the component class is gone
		endOfStructs = Ar.Tell();
		Ar.Seek(location_endOfStructs);
		Ar << endOfStructs;
		Ar.Seek(endOfStructs);
	}

	return true;
}

void FComponentStorage::PostSerialize(const FArchive& Ar)
{
	_objectIndexToElement.clear();

	for (size_t i = 0; i < _size; ++i)
	{
		FGraphObject * object = at(i, componentClass);

		_objectIndexToElement.emplace(object->nodeIndex, i);
	}
}

bool FEdgeStorage::Serialize(FArchive& Ar)
{
	Ar.UsingCustomVersion(FGraphVersion::GUID);

	Ar << _scriptStruct;
	Ar << _size;

	_structSize = _scriptStruct ? _scriptStruct->GetStructureSize() : 0;

	// endObjectObjects stores the end of this edge-storage block in the archive.
	// We'll use it if we don't have the _scripStruct anymore, on load.
	const int64 location_endOfObjects = Ar.Tell();
	int64 endOfObjects = 0;
	Ar << endOfObjects;

	if (_scriptStruct == nullptr && Ar.IsLoading())
	{
		// If we lost the componentClass, skip over whatever might have been in the archive
		Ar.Seek(endOfObjects);

		_size = 0;
		_objectIndexToHandle.Empty();
		_objectIndexIsValid.Empty();
		_handleToObject.Empty();
		pool.resize(0, 0);

		return true;
	}

	Ar << _objectIndexToHandle;

	if (_scriptStruct)
		pool.resize(_size, _scriptStruct->GetStructureSize());

	// Serialize each edge object
	for (size_t i = 0; i < _size; ++i)
	{
		FGraphEdgeObject * object = at(i);

		// very important to call the C++ constructor
		if (Ar.IsLoading())
			_scriptStruct->GetCppStructOps()->Construct(object);

		_scriptStruct->SerializeItem(Ar, object, nullptr);

		FGraphEdgeHandle handle = _objectIndexToHandle[i];

		_handleToObject.Add(handle, i);
	}

	// If we are saving, we now know where the end of the edge storage block is. Write
	// it back to the header, at location_endOfObjects.
	if (_scriptStruct && Ar.IsSaving())
	{
		endOfObjects = Ar.Tell();
		Ar.Seek(location_endOfObjects);
		Ar << endOfObjects;
		Ar.Seek(endOfObjects);
	}

	return true;
}

void FEdgeStorage::PostSerialize(const FArchive& Ar)
{
	// I'm not sure this is necessary
	_handleToObject.Empty();

	for (size_t i = 0; i < _size; ++i)
	{
		FGraphEdgeHandle handle = _objectIndexToHandle[i];

		_handleToObject.Add(handle, i);
	}

	// hack for now
	_objectIndexIsValid.Init(true, _size);
}
