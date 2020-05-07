// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include <vector>
#include <stdint.h>
#include <limits>
#include <memory>
#include <cstdint>

#include <unordered_map>
#include <unordered_set>

#include "ObjectRecord.h"
#include "TypeId.h"
#include "Pool.h"
#include "GraphUtility.h"
#include "HalfEdgeMesh.h"
#include "GraphVersion.h"

#include "BitArray.h"

#include "Graph.generated.h"

class UGraphComponent;
class UObjectSimulation;
struct FGraph;

struct FGraphObject;
struct FGraphNode;

struct FGraphEdge;
struct FGraphEdgeObject;


// A convenience class for keeping simulation state objects associated with FGraphObject indices.
//
// The order of the elements in the array will change as elements are added or removed, however
// a secondary map is used to uniquely map each element to an unique index on emplace.
// The elements must have an nodeIndex property for this to work.
// This enables very fast iteration of elements for simulation, while also maintaining fast
// insert/delete and a mapping to FGraphObject indices.
template<typename ElementType>
struct MappedArray
{
private:
	TArray<ElementType> * _elements;
	std::unordered_map<int32, int32> _objectIndexToElement;

public:
	void init( TArray<ElementType>& elements )
	{
		_elements = &elements; 

		int32 n = elements.Num();

		_objectIndexToElement.clear();
		for(int32 elementIndex = 0; elementIndex < n; ++elementIndex)
		{
			ElementType& element = elements[elementIndex];
			_objectIndexToElement.emplace( element.nodeIndex, elementIndex );
		}
	}

	template<class... Args>
	ElementType& emplace( int32 nodeIndex, Args&&... args )
	{
		int32 elementIndex = _elements->Emplace( std::forward<Args>( args )... );
		ElementType& element = (*_elements)[elementIndex];
		element.nodeIndex = nodeIndex;

		_objectIndexToElement.emplace( nodeIndex, elementIndex );

		return (*_elements)[elementIndex];
	}

	void erase( int32 nodeIndex )
	{
		auto found = _objectIndexToElement.find( nodeIndex );
		if(found == _objectIndexToElement.end())
			return;

		int32 n = _elements->Num();

		if(n == 1)
		{
			_elements->RemoveAt( 0 );
			_objectIndexToElement.erase( found );
		}
		else
		{
			int32 elementIndex = found->second;

			// last element, we can just kill it
			if(elementIndex == n - 1)
				_elements->RemoveAt( elementIndex );
			else
			{
				// middle element, so we efficiently swap in the last element of the array
				// and remove the element while update the _objectIndexToElement map with
				// the swapped element index
				_elements->RemoveAtSwap( elementIndex );

				ElementType& element = (*_elements)[elementIndex];

				_objectIndexToElement[element.nodeIndex] = elementIndex;
			}

			_objectIndexToElement.erase( found );
		}
	}

	int32 elementIndex( int32 nodeIndex )
	{
		return _objectIndexToElement[nodeIndex];
	}

	ElementType& at( int32 nodeIndex )
	{
		return (*_elements)[_objectIndexToElement[nodeIndex]];
	}

	bool contains( int32 nodeIndex )
	{
		return _objectIndexToElement.find( nodeIndex ) != _objectIndexToElement.end();
	}

	void clear()
	{
		_objectIndexToElement.clear();
		_elements->Empty();
	}
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraphNodeHandle
{
	GENERATED_BODY()

public:
	FGraphNodeHandle() : index(-1) {}
	explicit FGraphNodeHandle( int32 i ) : index( i ) {}

	explicit FGraphNodeHandle( FGraphNode& node );

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 index = -1;

	inline explicit operator bool() const
	{
		return index >= 0;
	}

	FGraphNode& node(FGraph& graph);

	FGraphNode& operator() ( FGraph& graph );
	FGraphNode& operator() ( FGraph* graph );

	friend bool operator<( const FGraphNodeHandle& a, const FGraphNodeHandle& b ) 
	{
		return a.index < b.index;
	}

	bool operator==(const FGraphNodeHandle& other)
	{
		return index == other.index;
	}

	bool operator!=(const FGraphNodeHandle& other)
	{
		return index != other.index;
	}

	friend inline uint32 GetTypeHash(const FGraphNodeHandle& Key)
	{
		return GetTypeHash(Key.index);
	}

	static const FGraphNodeHandle null;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FGraphEdgeHandle
{
	GENERATED_BODY()

public:
	FGraphEdgeHandle() : index(-1) {}

	explicit FGraphEdgeHandle(int32 i) : index(i) {}


	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 index = -1;

	inline explicit operator bool() const
	{
		return index >= 0;
	}

	FGraphEdge& operator() (FGraph& graph);
	FGraphEdge& operator() (FGraph* graph);

	friend bool operator<(const FGraphEdgeHandle& a, const FGraphEdgeHandle& b) 
	{
		return a.index < b.index;
	}

	bool operator==(const FGraphEdgeHandle& other) 
	{
		return index == other.index;
	}

	bool operator!=(const FGraphEdgeHandle& other) 
	{
		return index != other.index;
	}

	friend inline uint32 GetTypeHash(const FGraphEdgeHandle& Key) 
	{
		return GetTypeHash(Key.index);
	}

	static const FGraphEdgeHandle null;
};

FORCEINLINE FArchive& operator<<(FArchive &Ar, FGraphEdgeHandle& handle)
{
	Ar << handle.index;

	return Ar;
}

namespace std
{
	template<> struct hash<FGraphEdgeHandle>
	{
		std::size_t operator()(const FGraphEdgeHandle& handle) const
		{
			return std::hash<decltype(handle.index)>()(handle.index);
		}
	};
}

inline bool operator==(const FGraphEdgeHandle& A, const FGraphEdgeHandle& B)
{
	return A.index == B.index;
}

typedef StructTypeManager<FGraphObject>::TypeId ComponentType;

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraphObject 
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 nodeIndex = 0;
	
	FGraphNodeHandle nodeHandle() { return FGraphNodeHandle(nodeIndex); }
	FGraphNode& node(FGraph& graph);

	bool _isValid = true;
	bool isValid() { return _isValid; }
	void invalidate() { _isValid = false; }

	static UScriptStruct* componentStruct( ComponentType type );

	static ComponentType componentType( UScriptStruct* typeStruct );
};

template<class TComponent>
inline static ComponentType componentType()
{
	static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must be derived from TGraphComponent.");

	return StructTypeManager<FGraphObject>::typeId<TComponent>(); 
}

template<class TComponent>
inline static UScriptStruct* componentStruct()
{
	static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must be derived from TGraphComponent.");

	ComponentType type = componentType<TComponent>();
	return StructTypeManager<FGraphObject>::structForType( type );
}



namespace std
{
	template<> struct hash<FGraphNodeHandle>
	{
		std::size_t operator()( const FGraphNodeHandle& handle ) const
		{
			return std::hash<decltype(handle.index)>()(handle.index);
		}
	};
}

inline bool operator==(const FGraphNodeHandle& A, const FGraphNodeHandle& B)
{
	return A.index == B.index;
}

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraphNode
{
	GENERATED_USTRUCT_BODY()

public:
	FGraphNode() : position( FVector::ZeroVector )
	{
		edges.Empty();
	};

	FGraphNode( FVector& position_in ) : position( position_in )
	{
		edges.Empty();
	}

	FGraphNode( FVector& position_in, FQuat& orientation_in, float scale_in ) : position( position_in ), orientation( orientation_in ), scale( scale_in )
	{
		edges.Empty();
		
	}

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 id = 0;

	FGraphNodeHandle handle() { return FGraphNodeHandle(id); }

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Meta = (MakeEditWidget = true), Category = "ShipEditor" )
	FVector position = FVector::ZeroVector;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	FQuat orientation = FQuat::Identity;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float scale = 1.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<int32> edges;

	template<class TComponent>
	TComponent& component(FGraph* graph);

	template<class TComponent>
	TComponent& component( FGraph& graph );

	FGraphObject* component( FGraph& graph, ComponentType type );

	template<class TComponent>
	bool hasComponent();

	bool hasComponent( ComponentType type );

	FGraphObject* addComponent( FGraph& graph, ComponentType type );

	template<class TComponent, class... TArgs>
	TComponent& addComponent( FGraph& graph, TArgs... args );

	template<class TComponent>
	void removeComponent( FGraph& graph );

	void removeComponent( FGraph& graph, ComponentType type );

	void removeComponents( FGraph& graph )
	{
        TArray<int32> copy(components);

		for(ComponentType type : copy)
			removeComponent( graph, type );
	}

	// Calls the func on each edge that has the desired TEdgeObject.
	// ```func(FGraphNodeHandle otherNode, TEdgeObject& edgeObject)```
	template<class TEdgeObject, typename Func>
	void each(FGraph& graph, Func func);

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	bool _isValid = true;

	bool isValid() const { return _isValid; };
	void invalidate() { _isValid = false; edges.Empty(); };

	// Component types attached to this node, this will be populated during post serialize.
	// For simplicity of dealing with added/removed component classes, during serialization, and
	// the changing component type values (they are only valid at runtime and can change from run-to-run),
	// this property is transient.
	UPROPERTY( EditAnywhere, BlueprintReadOnly, Transient, Category = "ShipEditor" )
	TArray<int32> components; 
	operator FGraphNodeHandle() { return FGraphNodeHandle( id ); }
};


typedef int32 NodeIndex;

USTRUCT( BlueprintType )
struct FComponentStorage
{
	GENERATED_USTRUCT_BODY()

protected:


	std::unordered_map<int32, int32> _objectIndexToElement;

public:
	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	FPoolBase pool;

	// The number of objects
	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	int32 _size = 0;

	UPROPERTY( EditAnywhere, BlueprintReadOnly, Category = "ShipEditor" )
	UScriptStruct * componentClass = nullptr;

	bool Serialize(FArchive& Ar);

	void PostSerialize( const FArchive& Ar );

	FComponentStorage& operator=( const FComponentStorage& other )
	{
		_size = other._size;
		componentClass = other.componentClass;

		pool = other.pool;

		_objectIndexToElement.clear();

		for(size_t i = 0; i < _size; ++i)
		{
			FGraphObject * object = at( i, componentClass );

			_objectIndexToElement.emplace( object->nodeIndex, i );
		}

		return *this;
	}

	bool Identical( const FComponentStorage* other, uint32 PortFlags ) const
	{
		if(!other)
			return false;

		if(componentClass != other->componentClass)
			return false;

		if(!componentClass)
			return true;

		if(_size != other->_size)
			return false;

		return pool == other->pool;
	}

	void AddStructReferencedObjects( class FReferenceCollector& collector ) const
	{
		auto unconstThis = const_cast<FComponentStorage*>(this);

		collector.AddReferencedObject( unconstThis->componentClass );

		if(componentClass)
		{
			for(size_t i = 0; i < _size; ++i)
			{
				FGraphObject * object = unconstThis->at( i, componentClass );

				componentClass->SerializeBin( collector.GetVerySlowReferenceCollectorArchive(), object );
			}


			//FSimpleObjectReferenceCollectorArchive objectReferenceCollector( nullptr, collector );

			//for(size_t i = 0; i < _size; ++i)
			//{
			//	FGraphObject * object = unconstThis->at( i, componentClass );

			//	componentClass->SerializeBin( objectReferenceCollector, object );
			//}
		}
	}

	size_t size()
	{
		return _size;
	}

	FGraphObject* emplace( FGraphNode& node, UScriptStruct * typeStruct )
	{
		size_t componentSize = typeStruct->GetStructureSize();
		auto nodeIndex = node.id;

		pool.resize( _size + 1, componentSize );

		FGraphObject * component = static_cast<FGraphObject*>(pool.get( _size, componentSize ));

		typeStruct->InitializeStruct( (uint8*)component ); 

		int32 componentIndex = _size;

		_size++;

		component->nodeIndex = nodeIndex;

		_objectIndexToElement.emplace( nodeIndex, componentIndex );

		return component;
	}


	template<class TComponent, class... Args>
	TComponent& emplace( FGraphNode& node, Args&&... args )
	{
		static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must derive from FGraphObject.");

		auto nodeIndex = node.id;

		Pool<TComponent> * typedPool = static_cast< Pool<TComponent>* >(&pool);

		typedPool->resize( _size + 1 );

		TComponent * component = new(typedPool->get( _size )) TComponent( std::forward<Args>( args )... );

		int32 componentIndex = _size;

		_size++;

		component->nodeIndex = nodeIndex;

		_objectIndexToElement.emplace( nodeIndex, componentIndex );

		return *component;
	}

	void erase( FGraphNodeHandle nodeIndex, ComponentType type )
	{
		auto found = _objectIndexToElement.find( nodeIndex.index );
		if(found == _objectIndexToElement.end())
			return;

		UScriptStruct * scriptStruct = FGraphObject::componentStruct( type );

		size_t componentSize = scriptStruct->GetStructureSize();

		if(_size == 1)
		{
			scriptStruct->DestroyStruct( pool.get( 0, componentSize ) );
			
			_size = 0;

			_objectIndexToElement.erase( found );
		}
		else
		{
			int32 componentIndex = found->second;

			// last element, we can just kill it
			if(componentIndex == _size - 1)
			{
				scriptStruct->DestroyStruct( pool.get( componentIndex, componentSize ) );
				_size--;
			}
			else
			{
				scriptStruct->DestroyStruct( pool.get( componentIndex, componentSize ) );

				// copy the last element to this spot in the pool
				FMemory::Memcpy(
					pool.get( componentIndex, componentSize ),
					pool.get( _size - 1, componentSize ),
					componentSize
				);

				_size--;

				FGraphObject * component = static_cast<FGraphObject*>(pool.get( componentIndex, componentSize ));

				_objectIndexToElement[component->nodeIndex] = componentIndex;
			}

			_objectIndexToElement.erase( found );
		}

		pool.resize( _size, componentSize );
	}

	template<class TComponent>
	void erase(FGraphNodeHandle nodeIndex)
	{
		static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must derive from FGraphObject.");

		se::Pool<TComponent> * typedPool = static_cast< se::Pool<TComponent>* >(&pool);

		auto found = _objectIndexToElement.find( nodeIndex );
		if(found == _objectIndexToElement.end())
			return;

		if(_size == 1)
		{
			typedPool->destroy( 0 );
			_size = 0;

			_objectIndexToElement.erase( found );
		}
		else
		{
			int32 componentIndex = found->second;

			// last element, we can just kill it
			if(componentIndex == _size - 1)
			{
				typedPool->destroy( componentIndex );
				_size--;
			}
			else
			{
				typedPool->destroy( componentIndex );

				// copy the last element to this spot in the pool
				FMemory::Memcpy(
					typedPool->get( componentIndex ),
					typedPool->get( _size - 1),
					sizeof( TComponent )
				);

				_size--;

				TComponent * component = typedPool->get( componentIndex );

				_objectIndexToElement[component->nodeIndex] = componentIndex;
			}

			_objectIndexToElement.erase( found );
		}

		typedPool->resize( _size );
	}

	int32 elementIndex(FGraphNodeHandle nodeHandle )
	{
		return _objectIndexToElement[nodeHandle.index];
	}



	template<class TComponent>
	TComponent& at( int32 componentIndex )
	{
		static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must derive from FGraphObject.");
		Pool<TComponent> * typedPool = static_cast< Pool<TComponent>* >(&pool);

		TComponent * component = typedPool->get( componentIndex );
		return *component;
	}

	FGraphObject* at( int32 componentIndex, UScriptStruct * componentStruct )
	{
		size_t size = componentStruct->GetStructureSize();

		return static_cast<FGraphObject*>(pool.get( componentIndex, size ));
	}

	FGraphObject* componentForNode(FGraphNodeHandle nodeHandle)
	{
		auto found = _objectIndexToElement.find(nodeHandle.index);

		if (found == _objectIndexToElement.end())
			return nullptr;

		return static_cast<FGraphObject*>(pool.get(found->second, componentClass->GetStructureSize()));
	}

	FGraphObject* componentForNode( FGraphNode& node, ComponentType type )
	{
		auto found = _objectIndexToElement.find( node.id );
		
		if(found == _objectIndexToElement.end())
			return nullptr;
		 
		UScriptStruct * typeStruct = FGraphObject::componentStruct( type );

		return static_cast<FGraphObject*>(pool.get( found->second, typeStruct->GetStructureSize() ));
	}

	bool contains( int32 nodeIndex )
	{
		return _objectIndexToElement.find( nodeIndex ) != _objectIndexToElement.end();
	}

	template<class TComponent>
	void clear()
	{
		static_assert(std::is_base_of<FGraphObject, TComponent>::value, "TComponent must derive from FGraphObject.");
		Pool<TComponent>& typedPool = static_cast<Pool<TComponent>>(pool);

		_objectIndexToElement.clear();

		for(size_t i = 0; i < _size; ++i)
		{
			typedPool.destroy( i );
		}

		_size = 0;
		typedPool.resize( 0, 0 );
	}
};

template<>
struct TStructOpsTypeTraits<FComponentStorage> : public TStructOpsTypeTraitsBase2<FComponentStorage>
{
	enum
	{
		WithSerializer = true,
		WithCopy = true,
		WithIdentical = true,
		//WithAddStructReferencedObjects = true,
		WithAddStructReferencedObjects = true,
		WithPostSerialize = true
	};
};

template<class TComponent>
struct TypedComponentStorage : FComponentStorage
{
	TComponent* begin()
	{
		Pool<TComponent> * typedPool = static_cast<Pool<TComponent>*>(&pool);

		return typedPool->begin();
	}

	TComponent* end()
	{
		Pool<TComponent> * typedPool = static_cast<Pool<TComponent>*>(&pool);

		return typedPool->end();
	}

	TComponent& componentForNode( FGraphNode& node )
	{
		ComponentType type = componentType<TComponent>();

		return static_cast<TComponent&>(*FComponentStorage::componentForNode( node, type ));
	}

	TComponent& componentForNode(FGraphNodeHandle node)
	{
		return static_cast<TComponent&>(*FComponentStorage::componentForNode(node));
	}

	// returns null if the node doesn't have the component
	TComponent* componentPtrForNode(FGraphNodeHandle node)
	{
		return static_cast<TComponent*>(FComponentStorage::componentForNode(node));
	}

	TComponent& at( int32 i )
	{
		return FComponentStorage::at<TComponent>( i );
	}

	TComponent& operator()(size_t i)
	{
		at(i);
	}


	TComponent& operator[] (size_t i)
	{
		return at(i);
	}

	void swap(size_t objectIndexA, size_t objectIndexB)
	{
		Pool<TComponent>& typedPool = static_cast<Pool<TComponent>&>(pool);

		auto& handleA = _objectIndexToElement[objectIndexA];
		auto& handleB = _objectIndexToElement[objectIndexB];

		std::swap(handleA, handleB);
		std::swap(*typedPool.get(objectIndexA), *typedPool.get(objectIndexB));
	}

	// predicate(TComponent&, TComponent&): Sorts the connection objects relative to the predicate.
	template<typename TSortPredicate>
	void sort(TSortPredicate predicate)
	{
		// Inspired by: sparse_set::sort in https://github.com/skypjack/entt/blob/master/src/entt/entity/sparse_set.hpp

		Pool<TComponent>& typedPool = static_cast<Pool<TComponent>&>(pool);

		// we'll sort based on the values of the typed pool into a separate sorted index
		std::vector<size_t> indices(_size);
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(), [&typedPool, &predicate](size_t a, size_t b)
		{
			return predicate(*typedPool.get(a), *typedPool.get(b));
		});

		// apply the sorted indices to our containers
		for (size_t i = 0, last = indices.size(); i < last; ++i)
		{
			auto current = i;
			auto next = indices[current];

			while (current != next)
			{
				const auto lhs = indices[current];
				const auto rhs = indices[next];

				swap(lhs, rhs);

				indices[current] = current;
				current = next;
				next = indices[current];
			}
		}
	}
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraphEdge
{
	GENERATED_USTRUCT_BODY()

public:
    FGraphEdge() {};
	FGraphEdge( int32 aIn, int32 bIn ) : a( aIn ), b( bIn ) {};
    FGraphEdge( int32 aIn, int32 bIn, int32 type_in, float halfExtents_in ) : a( aIn ), b( bIn ), type(type_in), halfExtents(halfExtents_in) {}; 

	bool _isValid = true;
	bool isValid() const { return _isValid; };
	void invalidate() { _isValid = false; };

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 a = 0;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 b = 0;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 type = 0;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float halfExtents = 0.0f;

	// given a node in this edge, returns the other node
	// will assert if the node is not in this edge
	FGraphNodeHandle other( FGraphNodeHandle node );
};



template<typename iterator_type, typename value_type>
class IsValidArray_iterator
{
	iterator_type _current;

public:
	IsValidArray_iterator( const iterator_type& begin) : _current( begin )
	{
		// move to the first valid item
		while((_current) && !(*_current).isValid())
			++_current;
	}

	const value_type& operator*() const
	{
		return *_current;
	}

	const value_type* operator->() const
	{
		return _current.operator->();
	}

	IsValidArray_iterator& operator++()
	{
		do
		{
			++_current;
		} while((_current) && !(*_current).isValid());

		return *this;
	}

	bool operator==( const IsValidArray_iterator& rhs ) const
	{
		return _current == rhs._current;
	}

	bool operator!=( const IsValidArray_iterator& rhs ) const
	{
		return !(operator==( rhs ));
	}

	/** conversion to "bool" returning true if the iterator has not reached the last element. */
	FORCEINLINE explicit operator bool() const
	{
		return (bool)_current;
	}
	/** inverse of the "bool" operator */
	FORCEINLINE bool operator !() const
	{
		return !(bool)*this;
	}

	int32 GetIndex() const
	{
		return _current.GetIndex();
	}
};

// A static struct with static methods to emplace/remove elements from an array without resizing the array.
// The indices of recycled elements are kept in an array.
//
// Unreal's property system does not support templated structs. So, this is a compromise on that where one
// passes the storage array, and the recycled index array to emplace and remove elements.
struct RecyclingArray
{
	template<typename ElementType, class... Args>
	static int32 emplace( TArray<ElementType>& elements, TSet<int32>& recycledIndices, Args&&... args )
	{
		int32 index = -1;

		int32 n = recycledIndices.Num();
		if(n)
		{
			auto iterator = recycledIndices.CreateIterator();

			index = *iterator;

			iterator.RemoveCurrent();

			new(elements.GetData() + index) ElementType( std::forward<Args>(args)... );
		}
		else
		{
			index = elements.Num();
			elements.Emplace( std::forward<Args>( args )... );
		}
			
		return index;
	}

	// returns true if we invalidated the object at index, returns false if the object was already invalidated 
	template<typename ElementType>
	static bool remove( TArray<ElementType>& elements, TSet<int32>& recycledIndices, int32 index )
	{
		ElementType& e = elements[index];

		e.invalidate();

		recycledIndices.Add( index );

		return true;
	}
};


/**
 * View
 */

template<typename TComponent>
class View final 
{
	friend struct FGraph;

public:
	template<typename Func>
	void each(Func func)  
	{
		std::for_each(storage->begin(), storage->end(), [&, func = std::move(func)](auto& component) mutable {
			FGraphNode& node = graph->node(component.nodeHandle());
			func(node, component);
		});
	}

	TypedComponentStorage<TComponent> * storage;
	FGraph * graph;
};


typedef StructTypeManager<FGraphEdgeObject>::TypeId EdgeObjectType;

USTRUCT(BlueprintType)
struct FGraphEdgeObject
{
	GENERATED_BODY()
public:
	static UScriptStruct* edgeStruct(EdgeObjectType type);

	static EdgeObjectType edgeType(UScriptStruct* typeStruct);

};

template<class TEdgeObject>
inline static EdgeObjectType edgeType()
{
	static_assert(std::is_base_of<FGraphEdgeObject, TEdgeObject>::value, "TEdgeObject must be derived from FGraphEdgeObject.");

	static auto theType = StructTypeManager<FGraphEdgeObject>::typeId<TEdgeObject>();

	return theType;
}

USTRUCT( BlueprintType )
struct FEdgeStorage
{
	friend FGraph;

	GENERATED_BODY()

protected:
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	UScriptStruct * _scriptStruct = nullptr;
	size_t _structSize;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	int32 _size = 0;

	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	int32 _validSize = 0;

	// The pool is transient, we actually write directly to the archive during serialization, and 
	// reconstruct the pool during serialization.
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Transient, Category = "ShipEditor")
	FPoolBase pool;

	// A parallel array, where each entry is a handle for the corresponding edge-object
	UPROPERTY(EditAnywhere, BlueprintReadOnly, Category = "ShipEditor")
	TArray<FGraphEdgeHandle> _objectIndexToHandle;

	// A parallel array, to store whether a graph edge object is valid. It is not valid, if it is
	// removed, but not yet delete (because all listeners haven't been processed).

	TBitArray<> _objectIndexIsValid;

	// Not a property, we can reconstruct this during serialization
	TMap<FGraphEdgeHandle, int32> _handleToObject;

public:
	// The number of valid objects
	size_t validSize() { return _validSize; }
	// The total allocated number of objects, with invalids included. Use this for enumeration in an index based for loop.
	size_t totalSize() { return _size; }

	UScriptStruct * scriptStruct() { return _scriptStruct; }

public: // Access
	inline FGraphEdgeObject * add(FGraphEdgeHandle handle, FGraphEdgeObject& objectToCopy)
	{
		pool.resize(_size + 1, _structSize);

		FGraphEdgeObject * object = static_cast<FGraphEdgeObject*>(pool.get(_size, _structSize));

		_scriptStruct->CopyScriptStruct(object, &objectToCopy);

		int32 objectIndex = _size;

		_size++;
		_validSize++;

		_handleToObject.Add(handle, objectIndex);
		_objectIndexToHandle.Add(handle);
		_objectIndexIsValid.Add(true);

		return object;
	}

	inline FGraphEdgeObject * emplace(FGraphEdgeHandle handle)
	{
		pool.resize(_size + 1, _structSize);

		FGraphEdgeObject * object = static_cast<FGraphEdgeObject*>(pool.get(_size, _structSize));

		_scriptStruct->InitializeStruct((uint8*)object);

		int32 objectIndex = _size;

		_size++;
		_validSize++;

		_handleToObject.Add(handle, objectIndex);
		_objectIndexToHandle.Add(handle);
		_objectIndexIsValid.Add(true);

		return object;
	}

	// Erases elements at the object-index associated with the handle, through a swap
	inline void erase(FGraphEdgeHandle handle)
	{
		int32 objectIndex = -1;
		bool found = _handleToObject.RemoveAndCopyValue(handle, objectIndex);
		_objectIndexToHandle.RemoveAtSwap(objectIndex);

		auto valid = _objectIndexIsValid[objectIndex];

		_objectIndexIsValid.RemoveAtSwap(objectIndex);

		checkfSlow(found, TEXT("Tried to remove an object that doesn't exist"));

		if (_size == 1)
		{
			_scriptStruct->DestroyStruct(pool.get(0, _structSize));

			_size = 0;
			_validSize = 0;
		}
		else
		{
			// last element, we can just kill it
			if (objectIndex == _size - 1)
			{
				_scriptStruct->DestroyStruct(pool.get(objectIndex, _structSize));
				_size--;
				if (valid) _validSize--;
			}
			else
			{
				_scriptStruct->DestroyStruct(pool.get(objectIndex, _structSize));

				// Swap the last element to this spot in the pool
				FMemory::Memcpy(
					pool.get(objectIndex, _structSize),
					pool.get(_size - 1, _structSize),
					_structSize
				);

				_size--;
				if (valid) _validSize--;

				// update maps
				FGraphEdgeHandle oldHandle = _objectIndexToHandle[objectIndex];

				_handleToObject.Add(oldHandle, objectIndex);
			}
		}

		pool.resize(_size, _structSize);
	}

	inline void clear()
	{
		for (int i = 0; i < _size; ++i)
		{
			FGraphEdgeObject * edgeObject = static_cast<FGraphEdgeObject*>(pool.get(i, _structSize));

			_scriptStruct->DestroyStruct(edgeObject);
		}

		pool.resize(0, _structSize);

		_handleToObject.Reset();
		_objectIndexToHandle.Reset();
		_objectIndexIsValid.Reset();

		_size = 0;
		_validSize = 0;
	}

public:
	inline FGraphEdgeObject * at(int32 objectIndex)
	{
		size_t size = _scriptStruct->GetStructureSize();

		return static_cast<FGraphEdgeObject*>(pool.get(objectIndex, size));
	}

	inline FGraphEdgeObject * at(FGraphEdgeHandle handle)
	{
		size_t size = _scriptStruct->GetStructureSize();

		auto objectIndex = _handleToObject.Find(handle);

		if (!objectIndex)
			return nullptr;

		return static_cast<FGraphEdgeObject*>(pool.get(*objectIndex, size));
	}

	inline bool contains(FGraphEdgeHandle handle)
	{
		return _handleToObject.Contains(handle);
	}

	int32 objectIndex(FGraphEdgeObject * object)
	{
		auto objectI = reinterpret_cast<std::uintptr_t>(object);
		auto beginI = reinterpret_cast<std::uintptr_t>(pool.begin(_structSize));

		size_t index = (objectI - beginI) / _structSize;

		checkfSlow(index < _size, TEXT("The object is not part of the edge storage you nunny."));

		return index;
	}

	bool isValid(FGraphEdgeObject * object)
	{
		auto index = objectIndex(object);

		return _objectIndexIsValid[index];
	}

	bool isValid(FGraphEdgeHandle handle)
	{
		auto objectIndex = _handleToObject.Find(handle);

		if (objectIndex == nullptr) return false;

		return _objectIndexIsValid[*objectIndex];
	}

	bool isValid(int32 objectIndex)
	{
		return _objectIndexIsValid[objectIndex];
	}

	void invalidate(FGraphEdgeObject * object)
	{
		auto index = objectIndex(object);

		auto valid = _objectIndexIsValid[index];

		if (valid) _validSize--;

		_objectIndexIsValid[index] = false;
	}

	void invalidate(FGraphEdgeHandle handle)
	{
		auto index = _handleToObject.Find(handle);

		checkfSlow(index != nullptr, TEXT("The object is not part of the edge storage."));

		auto valid = _objectIndexIsValid[*index];

		if (valid) _validSize--;

		_objectIndexIsValid[*index] = false;
	}


public: // Serialization
	bool Serialize(FArchive& Ar);

	void PostSerialize(const FArchive& Ar);

public: // Other traits
	void AddStructReferencedObjects(class FReferenceCollector& collector) const
	{
		auto unconstThis = const_cast<FEdgeStorage*>(this);

		collector.AddReferencedObject(unconstThis->_scriptStruct);

		if (_scriptStruct)
		{
			for (size_t i = 0; i < _size; ++i)
			{
				FGraphEdgeObject * object = unconstThis->at(i);

				_scriptStruct->SerializeBin(collector.GetVerySlowReferenceCollectorArchive(), object);
			}
		}
	}

	bool Identical(const FEdgeStorage* other, uint32 PortFlags) const
	{
		if (!other)
			return false;

		if (_scriptStruct != other->_scriptStruct)
			return false;

		if (!_scriptStruct)
			return true;

		if (_size != other->_size)
			return false;

		return pool == other->pool;
	}
};

template<>
struct TStructOpsTypeTraits<FEdgeStorage> : public TStructOpsTypeTraitsBase2<FEdgeStorage>
{
	enum
	{
		WithSerializer = true,
		WithCopy = true,
		WithIdentical = true,
		//WithAddStructReferencedObjects = true,
		WithAddStructReferencedObjects = true,
		WithPostSerialize = true
	};
};



//// Implements observer-forwarding for addition/removal events on top of FEdgeStorage.
//USTRUCT(BlueprintType)
//struct FForwardingEdgeStorage : public FEdgeStorage
//{
//	friend FGraph;
//
//	GENERATED_BODY()
//public:
//	FGraph * graph;
//
//public:
//	bool contains(FGraphEdgeHandle handle)
//	{
//		return FEdgeStorageBase::contains(handle);
//	}
//
//	FGraphEdgeObject * add(FGraphEdgeHandle handle, FGraphEdgeObject& objectToCopy)
//	{
//		FGraphEdgeObject * object = FEdgeStorageBase::add(handle, objectToCopy);
//
//		graph->_connectionObjectAdded(handle, _edgeObjectType);
//
//		return object;
//	}
//
//	FGraphEdgeObject * emplace(FGraphEdgeHandle handle)
//	{
//		FGraphEdgeObject * object = FEdgeStorageBase::emplace(handle);
//
//		graph->_connectionObjectAdded(handle, _edgeObjectType);
//
//		return object;
//	}
//
//	// Erases elements at the object-index associated with the handle, through a swap
//	void erase(FGraphEdgeHandle handle)
//	{
//		FEdgeStorageBase::erase(handle);
//
//		graph->_connectionObjectRemoved(handle, _edgeObjectType);
//	}
//
//	void clear()
//	{
//		size_t structSize = _scriptStruct->GetStructureSize();
//
//		for (int i = 0; i < _size; ++i)
//		{
//			FGraphEdgeObject * edgeObject = static_cast<FGraphEdgeObject*>(pool.get(i, structSize));
//
//			_scriptStruct->DestroyStruct(edgeObject);
//		}
//
//		pool.resize(0, structSize);
//
//		_handleToObject.Reset();
//		_objectIndexToHandle.Reset();
//	}
//}



template<class TObject>
struct TypedEdgeStorage : FEdgeStorage
{
	TObject * begin()
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		return typedPool.begin();
	}

	TObject * end()
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		return typedPool.end();
	}

	template<class... Args>
	TObject& emplaceOrReplace(FGraphEdgeHandle handle, bool& isNew_out, Args&&... args)
	{
		auto objectIndex = _handleToObject.Find(handle);

		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		TObject * object = nullptr;

		if (objectIndex)
		{
			_scriptStruct->DestroyStruct(object);

			object = new(typedPool.get(*objectIndex)) TObject(std::forward<Args>(args)...);

			auto isValid = _objectIndexIsValid[*objectIndex];

			if (!isValid) _validSize++;
			isValid = true;

			isNew_out = false;
		}
		else
		{
			typedPool.resize(_size + 1);

			object = new(typedPool.get(_size)) TObject(std::forward<Args>(args)...);

			int32 objectIndex = _size;

			_size++;
			_validSize++;

			_objectIndexToHandle.Add(handle);
			_objectIndexIsValid.Add(true);
			_handleToObject.Add(handle, objectIndex);

			isNew_out = true;
		}

		return *object;
	}

	TObject * objectPtr(FGraphEdgeHandle handle)
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		auto * index = _handleToObject.Find(handle);

		if (!index)
			return nullptr;
		else
			return typedPool.get(*index);
	}

	TObject& operator[](size_t i)
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		return typedPool[i];
	}

	auto at(size_t i) -> std::pair< std::reference_wrapper<TObject>, FGraphEdgeHandle>
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		TObject& object = typedPool[i];

		return std::make_pair(std::reference_wrapper<TObject>(object), _objectIndexToHandle[i]);
	}

	FGraphEdgeHandle edge(size_t i)
	{
		return _objectIndexToHandle[i];
	}

	void swap(size_t objectIndexA, size_t objectIndexB)
	{
		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		auto bitReferenceA = _objectIndexIsValid[objectIndexA];
		auto bitReferenceB = _objectIndexIsValid[objectIndexB];

		bool temp = bitReferenceA;
		bitReferenceA = bool(bitReferenceB);
		bitReferenceB = temp;

		auto& handleA = _objectIndexToHandle[objectIndexA];
		auto& handleB = _objectIndexToHandle[objectIndexB];

		_handleToObject[handleA] = objectIndexB;
		_handleToObject[handleB] = objectIndexA;
				
		std::swap(handleA, handleB);
		std::swap(*typedPool.get(objectIndexA), *typedPool.get(objectIndexB));
	}

	// predicate(TObject&, TObject&): Sorts the connection objects relative to the predicate.
	template<typename TSortPredicate>
	void sort(TSortPredicate predicate)
	{
		// Inspired by: sparse_set::sort in https://github.com/skypjack/entt/blob/master/src/entt/entity/sparse_set.hpp

		Pool<TObject>& typedPool = static_cast<Pool<TObject>&>(pool);

		// we'll sort based on the values of the typed pool into a separate sorted index
		std::vector<size_t> indices(_size);
		std::iota(indices.begin(), indices.end(), 0);

		std::sort(indices.begin(), indices.end(), [&typedPool, &predicate, this](size_t a, size_t b)
		{
			bool aValid = isValid(a);
			bool bValid = isValid(b);

			if (aValid && bValid)
				return predicate(*typedPool.get(a), *typedPool.get(b));
			else if (aValid && !bValid)
				return true;
			else if (!aValid && bValid)
				return false;
			else 
				return a < b;		
		});

		// apply the sorted indices to our containers
		for (size_t i = 0, last = indices.size(); i < last; ++i)
		{
			auto current = i;
			auto next = indices[current];

			while (current != next) 
			{
				const auto lhs = indices[current];
				const auto rhs = indices[next];

				swap(lhs, rhs);

				indices[current] = current;
				current = next;
				next = indices[current];
			}
		}
	}
};

template<class TEdgeObject>
struct TypedImmutableEdgeStorage
{
protected:
	TypedEdgeStorage<TEdgeObject> & storage;

public:
	TypedImmutableEdgeStorage(TypedEdgeStorage<TEdgeObject>& storage) : storage(storage) {}

	TEdgeObject * begin()
	{
		return storage.begin();
	}

	TEdgeObject * end()
	{
		return storage.end();
	}

	TEdgeObject * objectPtr(FGraphEdgeHandle handle)
	{
		return storage.objectPtr(handle);
	}

	bool isValid(FGraphEdgeHandle handle)
	{
		return storage.isValid(handle);
	}

	bool isValid(TEdgeObject& object)
	{
		return storage.isValid(&object);
	}

	bool isValid(int32 index)
	{
		return storage.isValid(index);
	}

	TEdgeObject& operator[](size_t i)
	{
		return storage[i];
	}

	auto at(size_t i) -> std::pair< std::reference_wrapper<TEdgeObject>, FGraphEdgeHandle>
	{
		return storage.at(i);
	}

	size_t validSize()
	{
		return storage.validSize();
	}

	size_t totalSize()
	{
		return storage.totalSize();
	}
};


template<typename TEdgeObject>
class EdgeView final
{
	friend struct FGraph;

public:
	// func(FGraphEdgeObject&, FGraphEdgeHandle)
	template<typename Func>
	void eachHandle(Func func)
	{
		size_t n = storage.totalSize();

		for (int32 i = 0; i < n; ++i)
		{
			auto pair = storage.at(i);
			FGraphEdge& edge = edges[pair.second.index];

			if( edge.isValid() )
				func(pair.first, FGraphEdgeHandle(pair.second.index));
		}
	}

	// func(FGraphEdgeObject&, FGraphEdge&)
	template<typename Func>
	void each(Func func)
	{
		size_t n = storage.totalSize();

		for (int32 i = 0; i < n; ++i)
		{
			if (storage.isValid(i))
			{
				auto pair = storage.at(i);
				FGraphEdge& edge = edges[pair.second.index];

				if (edge.isValid())
					func(pair.first, edge);
			}
		}
	}

	TypedEdgeStorage<TEdgeObject>& storage;
	TArray<FGraphEdge>& edges;

	EdgeView(TypedEdgeStorage<TEdgeObject>& storage_in, TArray<FGraphEdge>& edges_in) : storage(storage_in), edges(edges_in) {}
};

struct NodeListener
{
	virtual void nodeAdded(FGraphNode& node) {};
	virtual void nodeUpdated(FGraphNode& node) {};
	virtual void nodeRemoved(FGraphNode& oldNode) {};
};

struct ComponentListener
{
	virtual void componentAdded(FGraphNodeHandle node, ComponentType type) {};
	virtual void componentRemoved(FGraphNodeHandle node, ComponentType type) {};
	virtual void componentUpdated(FGraphNodeHandle node, ComponentType type) {};

	virtual void connectionAdded(int32 edgeIndex, FGraphEdge& edge, ComponentType type) {};
	virtual void connectionRemoved(int32 edgeIndex, FGraphEdge& oldEdge, ComponentType type) {};
};

struct EdgeListener
{
	virtual void connectionAdded(int32 edgeIndex, FGraphEdge& edge) {};
	virtual void connectionUpdated(int32 edgeIndex, FGraphEdge& edge) {};
	virtual void connectionRemoved(int32 edgeIndex, FGraphEdge& oldEdge) {};
};

struct EdgeObjectListener
{
	virtual void edgeObjectAdded(FGraphEdgeHandle handle, EdgeObjectType type) {};
	virtual void edgeObjectRemoved(FGraphEdgeHandle handle, EdgeObjectType type) {};
};

/**
 * Graph
 */
USTRUCT( BlueprintType )
struct LIFEBRUSH_API FGraph
{
	GENERATED_BODY()

public:
	typedef IsValidArray_iterator<TArray<FGraphNode>::TConstIterator, FGraphNode> nodeIterator;
	typedef IsValidArray_iterator<TArray<FGraphEdge>::TConstIterator, FGraphEdge> edgeIterator;

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<FGraphNode> allNodes;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TSet<int32> recycledNodes;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<FGraphEdge> privateEdges;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TSet<int32> recycledEdges;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<FComponentStorage> _componentStorage;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<FEdgeStorage> _edgeStorage;

	int32 tickCount = 0;

	TSet<FGraphNodeHandle> unrealEditorSelection;

public:
	void init();

	void clear();


	void resetListeners()
	{
		_transactionContexts.erase(_transactionContexts.begin() + 1, _transactionContexts.end());

		_transactionContexts[0].resetListeners();
	}

	// -------------------------------------------------------------------------------
	// Components
	// -------------------------------------------------------------------------------

	FComponentStorage& componentStorage( ComponentType type );

	template<class TComponent>
	TypedComponentStorage<TComponent>& componentStorage()
	{
		static const auto type = componentType<TComponent>();

		return static_cast<TypedComponentStorage<TComponent>&>(componentStorage( type ));
	}



	FGraphNodeHandle anyNode();

	// the returned index is only valid until the next change to the graph
	int32 addNode( FVector position, FQuat orientation = FQuat::Identity, float scale = 1.0f);
	void removeNode(FGraphNodeHandle node) { removeNode(node.index); }
	void removeNode( int32 index );

	FGraphNode& node( int32 index ) { return allNodes[index]; }
	FGraphNode& node(FGraphNodeHandle handle) { return allNodes[handle.index]; }

	template<class TComponent>
	TComponent& component(FGraphNodeHandle nodeHandle);

	template<class TComponent>
	TComponent * componentPtr(FGraphNodeHandle nodeHandle);

	FGraphObject * component(FGraphNodeHandle nodeHandle, ComponentType componentType)
	{
		return componentStorage(componentType).componentForNode(nodeHandle);
	}

    void markNodeDirty( int32 index );


	// -------------------------------------------------------------------------------
	// Edges
	// -------------------------------------------------------------------------------

	// This is a semi-private accessor. Don't use it unless you are updating the rest of the graph state.
	// It's here so that it is easy for other classes to build functionality on top of the graph, such
	// as implementing snapshots.
	FEdgeStorage& rawEdgeStorage(EdgeObjectType type);

	template<class TEdgeObject>
	TypedEdgeStorage<TEdgeObject>& rawEdgeStorage()
	{
		static const auto type = edgeType<TEdgeObject>();

		return static_cast<TypedEdgeStorage<TEdgeObject>&>(_edgeStorage[type]);
	}

	template<class TEdgeObject>
	TypedImmutableEdgeStorage<TEdgeObject> edgeStorage()
	{
		static const auto type = edgeType<TEdgeObject>();

		auto& rawStorage = rawEdgeStorage<TEdgeObject>();

		return TypedImmutableEdgeStorage<TEdgeObject>(rawStorage);
	}

	FGraphEdgeHandle connectNodes(FGraphNodeHandle a, FGraphNodeHandle b);

	template<class TEdgeObject, class... Args>
	FGraphEdgeHandle connectNodes(FGraphNodeHandle a, FGraphNodeHandle b, Args&&... args)
	{
		FGraphEdgeHandle handle = connectNodes(a, b);

		bool isNew = false;
		rawEdgeStorage<TEdgeObject>().emplaceOrReplace(handle, isNew, std::forward<Args>(args)...);

		static const auto type = edgeType<TEdgeObject>();

		_edgeObjectAdded(handle, type);

		return handle;
	}

	void removeConnection(FGraphEdgeHandle handle );

	FGraphEdge& edge(FGraphEdgeHandle handle ) { return privateEdges[handle.index]; };

	template<class TEdgeObject, class... Args>
	TEdgeObject& addOrReplaceEdgeObject(FGraphEdgeHandle handle, Args&&... args)
	{
		static const auto type = edgeType<TEdgeObject>();

		bool isNew = false;
		auto& object = rawEdgeStorage<TEdgeObject>().emplaceOrReplace(handle, isNew, std::forward<Args>(args)...);

		if( isNew )
			_edgeObjectAdded(handle, type);

		return object;
	}

	void removeEdgeObject(FGraphEdgeHandle handle, EdgeObjectType type)
	{
		auto& storage = rawEdgeStorage(type);

		FGraphEdgeObject * object = storage.at(handle);

		if (object && storage.isValid(object))
		{
			storage.invalidate(object);
			_edgeObjectRemoved(handle, type);
		}
	}

	template<class TEdgeObject>
	void removeEdgeObject(FGraphEdgeHandle handle)
	{
		static const auto type = edgeType<TEdgeObject>();

		auto& storage = rawEdgeStorage<TEdgeObject>();

		TEdgeObject * object = storage.objectPtr(handle);

		if (object && storage.isValid(object))
		{
			storage.invalidate(object);
			_edgeObjectRemoved(handle, type);
		}
	}

	template<class TEdgeObject>
	TEdgeObject& edgeObject(FGraphEdgeHandle handle)
	{
		auto ptr = rawEdgeStorage<TEdgeObject>().objectPtr(handle);

		checkfSlow(ptr, TEXT("The FGraphEdgeObject at handle %d does not exist"), handle.index);

		return *ptr;
	}

	// \return nullptr if the FGraphEdgeObject of the requested type is not created for the handle.
	template<class TEdgeObject>
	TEdgeObject* edgeObjectPtr(FGraphEdgeHandle handle)
	{
		return rawEdgeStorage<TEdgeObject>().objectPtr(handle);
	}

	template<class TEdgeObject>
	bool hasEdgeObject(FGraphEdgeHandle handle)
	{
		return rawEdgeStorage<TEdgeObject>().objectPtr(handle) != nullptr;
	}

	TArray<FGraphEdgeHandle> edgesBetweenNodes(FGraphNodeHandle a, FGraphNodeHandle b);

	TArray<FGraphEdgeHandle> edgesBetweenNodes(TArray<FGraphNodeHandle>& nodes);

	// -------------------------------------------------------------------------------
	// Misc.
	// -------------------------------------------------------------------------------



	// returns all of the valid nodes
	nodeIterator nodes() 
	{
		return nodeIterator( allNodes.CreateConstIterator() );
	}

	int32 numNodes()
	{
		return allNodes.Num() - recycledNodes.Num();
	}

	edgeIterator edges()
	{
		return edgeIterator( privateEdges.CreateConstIterator() ); 
	}

	int32 numEdges()
	{
		return privateEdges.Num() - recycledEdges.Num();
	}

	std::vector<FGraphNodeHandle> nodeHandles()
	{
		std::vector<FGraphNodeHandle> handles;
		handles.reserve(numNodes());

		for (FGraphNode& n : allNodes)
		{
			if (n.isValid())
				handles.push_back(n.handle());
		}

		return handles;
	}

	std::vector<FGraphEdgeHandle> edgeHandles()
	{
		std::vector<FGraphEdgeHandle> handles;
		handles.reserve(numEdges());

		for (int32 i = 0; i <privateEdges.Num(); ++i)
		{
			FGraphEdge& edge = privateEdges[i];

			if (edge.isValid())
				handles.push_back(FGraphEdgeHandle(i));
		}

		return handles;
	}

private:
	void _initComponentStorage();
	void _initEdgeStorage();

	void _assignComponents();

	// Notifies the transaction system of an added edge object.
	void _edgeObjectAdded(FGraphEdgeHandle handle, const  EdgeObjectType type);
	// Notifies the transaction system of a removed edge object.
	void _edgeObjectRemoved(FGraphEdgeHandle handle, EdgeObjectType type);

private:
	bool _didInit = false;


	// ---------------------------------------------------------------------------
	// Transaction API
	// ---------------------------------------------------------------------------
public:
	void addComponentListener(ComponentListener * listener, ComponentType type)
	{
		currentTransactionContext().addListener(listener, type);
	}

	template<class TComponent>
	void addComponentListener(ComponentListener * listener)
	{
		addComponentListener(listener, componentType<TComponent>());
	}

	void addNodeListener(NodeListener * listener)
	{
		currentTransactionContext().addListener(listener);
	}

	template<class TComponent>
	void removeComponentListener(ComponentListener * listener)
	{
		auto type = componentType<TComponent>();

		TransactionContext& context = currentTransactionContext();

		auto found = context._componentListeners.find(type);

		if (found != context._componentListeners.end())
		{
			auto& listeners = found->second;

			listeners.erase(std::remove_if(listeners.begin(), listeners.end(), [=](ComponentListener * l) {
				return l == listener;
			}), listeners.end());
		}
	}

	// -------------------------------------------------------------
	// Edge Listeners
	// -------------------------------------------------------------

	void addEdgeListener(EdgeListener * listener)
	{
		currentTransactionContext().addListener(listener);
	}

	void addEdgeObjectListener(EdgeObjectListener * listener, EdgeObjectType type)
	{
		currentTransactionContext().addListener(listener, type);
	}

	template<class TEdgeObject>
	void addEdgeObjectListener(EdgeObjectListener * listener)
	{
		static_assert(std::is_base_of<FGraphEdgeObject, TEdgeObject>::value, "TEdgeObject must be derived from FGraphEdgeObject.");

		addEdgeObjectListener(listener, edgeType<TEdgeObject>());
	}


	template<class TEdgeObject>
	void removeEdgeObjectListener(EdgeObjectListener * listener)
	{
		static_assert(std::is_base_of<FGraphEdgeObject, TEdgeObject>::value, "TEdgeObject must be derived from FGraphEdgeObject.");

		auto type = edgeType<TEdgeObject>();

		TransactionContext& context = currentTransactionContext();

		auto found = context._edgeObjectListeners.find(type);

		if (found != context._edgeObjectListeners.end())
		{
			auto& listeners = found->second;

			listeners.erase(std::remove_if(listeners.begin(), listeners.end(), [=](EdgeObjectListener * l) {
				return l == listener;
			}), listeners.end());
		}
	}

private:
	// Transactions group operations together and defer the applications of the operations.
    struct GraphTransaction
	{
		friend FGraph;
	protected:
		std::vector<int32> addedNodes;
		std::vector<int32> removedNodes;   // this has to be a copy, since the node may have already been overwritten
		std::vector<int32> dirtyNodes;

		std::unordered_map<ComponentType, std::vector<FGraphNodeHandle>> addedComponents;
		std::unordered_map<ComponentType, std::vector<FGraphNodeHandle>> removedComponents;

		std::vector<int32> addedConnections;
		std::vector<int32> removedConnections; // old index and a copy of the old edge

		std::unordered_map<EdgeObjectType, std::vector<FGraphEdgeHandle>> addedEdgeObjects;
		std::unordered_map<EdgeObjectType, std::vector<FGraphEdgeHandle>> removedEdgeObjects;

	public:

		void append(GraphTransaction& other)
		{
			for (auto& pair : other.addedComponents)
			{
				auto& handles = addedComponents[pair.first];

				handles.insert(handles.end(), pair.second.begin(), pair.second.end());
			}

			for (auto& pair : other.removedComponents)
			{
				auto& handles = removedComponents[pair.first];

				handles.insert(handles.end(), pair.second.begin(), pair.second.end());
			}

			for (auto& pair : other.removedEdgeObjects)
			{
				auto& handles = removedEdgeObjects[pair.first];

				handles.insert(handles.end(), pair.second.begin(), pair.second.end());
			}

			for (auto& pair : other.addedEdgeObjects)
			{
				auto& handles = addedEdgeObjects[pair.first];

				handles.insert(handles.end(), pair.second.begin(), pair.second.end());
			}

			addedConnections.insert(addedConnections.end(), other.addedConnections.begin(), other.addedConnections.end());
			removedConnections.insert(removedConnections.end(), other.removedConnections.begin(), other.removedConnections.end());

			addedNodes.insert(addedNodes.end(), other.addedNodes.begin(), other.addedNodes.end());
			removedNodes.insert(removedNodes.end(), other.removedNodes.begin(), other.removedNodes.end());
			dirtyNodes.insert(dirtyNodes.end(), other.dirtyNodes.begin(), other.dirtyNodes.end());
		}
	};

	// A transaction context dispatches transaction operations from the GraphTransactions in its stack.
	// The actual operations are pushed down the transaction stack when the context is popped.
	struct TransactionContext
	{
	public:
		std::vector<NodeListener*> _nodeListeners;
		std::unordered_map<ComponentType, std::vector<ComponentListener*>> _componentListeners;

		std::vector<EdgeListener*> _edgeListeners;
		std::unordered_map<EdgeObjectType, std::vector<EdgeObjectListener*>> _edgeObjectListeners;

		void addListener(NodeListener * listener)
		{
			_nodeListeners.push_back(listener);
		}

		void addListener(ComponentListener * listener, ComponentType type)
		{
			_componentListeners[type].push_back(listener);
		}

		void addListener(EdgeListener * listener)
		{
			_edgeListeners.push_back(listener);
		}

		void addListener(EdgeObjectListener * listener, EdgeObjectType type)
		{
			_edgeObjectListeners[type].push_back(listener);
		}

		void resetListeners()
		{
			_nodeListeners.clear();
			_componentListeners.clear();
			
			_edgeListeners.clear();
			_edgeObjectListeners.clear();
		}

		uint32 _graphTransactions = 0;
	};

	std::vector<GraphTransaction> _graphTransactions;

	// we always have at least one context
	std::vector<TransactionContext> _transactionContexts = { TransactionContext() };

	TransactionContext& currentTransactionContext()
	{
		return _transactionContexts.back();
	}

	bool hasActiveTransaction()
	{
		return currentTransactionContext()._graphTransactions > 0;
	}

	GraphTransaction * _activeTransaction()
	{
		TransactionContext& transactionContext = currentTransactionContext();

		uint32 n = transactionContext._graphTransactions;
		if (n)
			return &_graphTransactions.back();
		else
			return nullptr;
	}

	GraphTransaction * _backTransaction()
	{
		size_t n = _graphTransactions.size();

		if (n)
			return &_graphTransactions.back();
		else
			return nullptr;
	}

	GraphTransaction * _nextTransaction()
	{
		size_t n = _graphTransactions.size();

		if (n < 2)
			return nullptr;
		else
			return &_graphTransactions[n - 2];
	}

	void _beginTransaction()
	{
		_graphTransactions.emplace_back();

	}

	void _endTransaction()
	{
		TransactionContext& transactionContext = currentTransactionContext();

		GraphTransaction * back = _backTransaction();

		if (!back)
			return;

		auto& componentListeners = transactionContext._componentListeners;
		auto& edgeObjectListeners = transactionContext._edgeObjectListeners;
		auto& nodeListeners = transactionContext._nodeListeners;
		auto& edgeListeners = transactionContext._edgeListeners;

		// dispatch added components
		for (auto& pair : back->addedComponents)
		{
			auto type = pair.first;
			auto& handles = pair.second;

			auto found = componentListeners.find(type);

			if (found == componentListeners.end())
				continue;

			for (auto listener : found->second)
			{
				for (auto& handle : handles)
					listener->componentAdded(handle, type);
			}
		}

		// dispatch removed components
		for (auto& pair : back->removedComponents)
		{
			auto type = pair.first;
			auto& handles = pair.second;

			auto found = componentListeners.find(type);

			if (found == componentListeners.end())
				continue;

			for (auto listener : found->second)
			{
				for (auto& handle : handles)
					listener->componentRemoved(handle, type);
			}
		}

		// dispatch added graph-edge objects
		for (auto& pair : back->addedEdgeObjects)
		{
			auto type = pair.first;
			auto& edgeHandles = pair.second;

			auto found = edgeObjectListeners.find(type);

			if (found == edgeObjectListeners.end())
				continue;

			for (auto listener : found->second)
			{
				for (auto& handle : edgeHandles)
					listener->edgeObjectAdded(handle, type);
			}
		}

		// dispatch removed graph-edge objects
		for (auto& pair : back->removedEdgeObjects)
		{
			auto type = pair.first;
			auto& edgeHandles = pair.second;

			auto found = edgeObjectListeners.find(type);

			if (found == edgeObjectListeners.end())
				continue;

			for (auto listener : found->second)
			{
				for (auto& handle : edgeHandles)
					listener->edgeObjectRemoved(handle, type);
			}
		}

		// dispatch added connections
		for (auto edgeIndex : back->addedConnections)
			_connectionAdded(edgeIndex);

		// dispatch removed connections
		for (auto edgeIndex : back->removedConnections)
			_connectionRemoved(edgeIndex);

		// dispatch added nodes
		for (auto nodeIndex : back->addedNodes)
		{
			for (auto listener : nodeListeners)
				listener->nodeAdded(allNodes[nodeIndex]);
		}

		// dispatch removed nodes
		for (auto& nodeIndex : back->removedNodes)
		{
			for (auto listener : nodeListeners)
				listener->nodeRemoved(allNodes[nodeIndex]);
		}

		// dispatch dirty nodes
		for (auto nodeIndex : back->dirtyNodes)
		{
			for (auto listener : nodeListeners)
				listener->nodeUpdated(allNodes[nodeIndex]);
		}

		GraphTransaction * nextTransaction = _nextTransaction();

		if (nextTransaction)
		{
			// trickle down transactions to the next transaction in the stack
			nextTransaction->append(*back);
		}
		else
		{
			// if this is the last transaction, then its time to purge
			_purgeTransaction(*back);
		}

		_graphTransactions.pop_back();
	}

public:
	// A transaction context dispatches componentAdded/nodeAdded/connectionAdded events to all listeners registered
	// between a pushTransactionContext and popTransactionContext.
	// All of the events between, are also saved, on popTransactionContext, they are passed to the next transaction
	// context on the stack.
	void pushTransactionContext()
	{
		// every context has an invisible transaction to collect events, when the context is
		// popped, we dispatch those events, if there are no beginTransaction s for the current
		// context, then immediate dispatch is performed, for the listeners of the context
		_beginTransaction();

		_transactionContexts.emplace_back();
	}

	void popTransactionContext()
	{
		if (_transactionContexts.size() == 1)
			std::logic_error("too many popTransactionContext calls");

		if (_transactionContexts.back()._graphTransactions)
			std::logic_error("the transaction context still has GraphTransactions. Missing endTransaction.");

		_transactionContexts.pop_back();

		// We must do this after we pop the _transactionContexts.
		// The idea: we want the collected GraphTransaction events to pop to the next transaction in the stack, not
		// to the one we just popped.
		_endTransaction();
	}

	void beginTransaction()
	{
		TransactionContext& transactionContext = currentTransactionContext();

		transactionContext._graphTransactions++;

		_beginTransaction();
	}

	void endTransaction()
	{
		_endTransaction();

		TransactionContext& transactionContext = currentTransactionContext();

		if (transactionContext._graphTransactions <= 0)
			std::logic_error("too many endTransactions called, for the current TransactionContext. Did you call endTransaction after a popTransactionContext?");

		transactionContext._graphTransactions--;
	}

	void didUpdateNode( NodeIndex index )
	{
		//GraphTransaction& back = _graphTransactions.back();
		//
		//back.dirtyNodes.push_back( index );

		UE_LOG( LogTemp, Warning, TEXT( "didUpdateNode not implemented in FGraph" ) );
	}

protected:
	// removes all of the nodes, edges and components referenced in the transaction
	// This is destructive to the GraphTransaction object, it's state will change, so
	// only call this once and don't rely on the validity of the GraphTransaction afterwards.
	void _purgeTransaction(GraphTransaction& transaction);

public:

	void componentAdded( FGraphNode& node, ComponentType type )
	{
		GraphTransaction * backTransaction = _backTransaction();

		// transactional dispatch
		if (backTransaction)
			backTransaction->addedComponents[type].emplace_back(node.id);
		
		if (!hasActiveTransaction())
		{
			TransactionContext& transactionContext = currentTransactionContext();

			auto& componentListeners = transactionContext._componentListeners;

			auto found = componentListeners.find( type );

			if(found == componentListeners.end())
				return;

			for(auto listener : found->second)
				listener->componentAdded( node, type );
		}
	}

	void componentRemoved( FGraphNode& node, ComponentType type )
	{
		GraphTransaction * backTransaction = _backTransaction();

		// transactional dispatch
		if(backTransaction)
			backTransaction->removedComponents[type].emplace_back(node.id);
		
		if(!hasActiveTransaction())
		{
			TransactionContext& transactionContext = currentTransactionContext();

			auto& componentListeners = transactionContext._componentListeners;

			auto found = componentListeners.find(type);

			if(found == componentListeners.end())
			{
				for(auto listener : found->second)
					listener->componentRemoved( node, type );
			}

			componentStorage( type ).erase( node, type );
		}
	}

	void _connectionAdded( int32 edgeIndex )
	{
		TransactionContext& context = currentTransactionContext();

		FGraphEdge& edge = privateEdges[edgeIndex];

		FGraphNode& nodeA = allNodes[edge.a];
		FGraphNode& nodeB = allNodes[edge.b];

		std::unordered_set<ComponentType> dispatchedTypes;

		auto& componentListeners = context._componentListeners;
		auto& connectionListeners = context._edgeListeners;

        // dispatch to typed listeners
		for(ComponentType type : nodeA.components)
		{
			auto found = context._componentListeners.find(type);

			if(found == componentListeners.end())
				continue;

			for( auto listener : found->second )
			{
				listener->connectionAdded( edgeIndex, edge, type);
				dispatchedTypes.insert( type );
			}
		}

		for(ComponentType type : nodeB.components)
		{
			auto found = componentListeners.find( type );

			if( dispatchedTypes.find(type) != dispatchedTypes.end() || found == componentListeners.end())
				continue;

			for(auto listener : found->second)
				listener->connectionAdded( edgeIndex, edge, type );
		}

        // dispatch to regular listeners
        for(auto listener : connectionListeners)
            listener->connectionAdded( edgeIndex, edge );
    }

	void _connectionRemoved( int32 edgeIndex )
	{
		TransactionContext& context = currentTransactionContext();

		auto& componentListeners = context._componentListeners;
		auto& connectionListeners = context._edgeListeners;

        FGraphEdge& edge = privateEdges[edgeIndex];

        FGraphNode& nodeA = allNodes[edge.a];
        FGraphNode& nodeB = allNodes[edge.b];

        std::unordered_set<ComponentType> dispatchedTypes;

        // dispatch to typed listeners
        for(ComponentType type : nodeA.components)
        {
            auto found = componentListeners.find( type );

            if(found == componentListeners.end())
                continue;

            for(auto listener : found->second)
            {
                listener->connectionAdded( edgeIndex, edge, type );
                dispatchedTypes.insert( type );
            }
        }

        for(ComponentType type : nodeB.components)
        {
            auto found = componentListeners.find( type );

            if(dispatchedTypes.find( type ) != dispatchedTypes.end() || found == componentListeners.end())
                continue;

            for(auto listener : found->second)
                listener->connectionAdded( edgeIndex, edge, type );
        }

        // dispatch to regular listeners
        for(auto listener : connectionListeners)
            listener->connectionAdded( edgeIndex, edge );
	}



	//---------------------------------------------------------------------------
	// TStructOpsTypeTraitsBase
	//---------------------------------------------------------------------------

	FGraph& operator=( const FGraph& Other );

	// --------------------------------------------------------------
	// Serialization
	// --------------------------------------------------------------
	void PostSerialize( const FArchive& Ar ) 
	{
		// the order of component storage containers may have changed due to componentType<TComponent>()
		TArray<UScriptStruct*> derived;
		ShipEditorUtility::derivedStructs( FGraphObject::StaticStruct(), derived );

		// before trying to get the component type for these structs, we have to register them with typeForStruct
		for(UScriptStruct * scriptStruct : derived)
			FGraphObject::componentType( scriptStruct );

		// first prune any component storage that doesn't have a componentClass (the class was renamed or removed)
		_componentStorage.RemoveAll([](FComponentStorage& storage) {
			return storage.componentClass == nullptr;
		});

		// determine the new classes to add
		{
			std::vector<UScriptStruct*> toAdd;

			for(UScriptStruct * componentClass : derived)
			{
				FComponentStorage* found = _componentStorage.FindByPredicate( [&]( FComponentStorage& item )->bool
				{
					return item.componentClass == componentClass;
				} );

				if(found == nullptr)
					toAdd.push_back( componentClass );
			}

			// add them
			for(UScriptStruct * componentClass : toAdd)
			{
				int32 newIndex = _componentStorage.Emplace();

				_componentStorage[newIndex].componentClass = componentClass;
			}
		}

		// sort the component array by componentTypes
		// if we allocated some new component storage, we'll deal with it after, but we'll sort it to the end of the array
		_componentStorage.Sort( [&]( const FComponentStorage& a, const FComponentStorage& b ) {
			ComponentType typeA = FGraphObject::componentType( a.componentClass );
			ComponentType typeB = FGraphObject::componentType( b.componentClass );

			return typeA < typeB;
		} );

		//// for some reason, Unreal will have called PostSerialize, when only the component storage has been set on the graph, but the nodes
		//// are missing (presumably they are going to be copied, rather than serialized), in which case we are in trouble!
		//// abort, if there are no nodes
		//if(allNodes.Num() == 0)
		//	return;

		// update the graph node component arrays
		for(auto& graphNode : allNodes)
			graphNode.components.Empty();
		
		for(auto& storage : _componentStorage)
		{
			ComponentType type = FGraphObject::componentType( storage.componentClass );

			for(int i = 0; i < storage.size(); ++i)
			{
				FGraphObject * object = storage.at( i, storage.componentClass );

				FGraphNode& graphNode = node( object->nodeIndex );

				graphNode.components.AddUnique( type );
			}
		}
	}

public:
	//template<TComponent,...>
	//View<FGraphNodeHandle, TComponent...> view()
	//{
	//	return { &componentStorage<TComponent>()... };
	//}
	template<typename TComponent>
	View<TComponent> view()
	{
		View<TComponent> v = { &(componentStorage<TComponent>()), this };

		return v;
	}

	// Enumerates each valid FGraphNodeObject and corresponding FGraphNode
	template<typename TNodeObject, typename Func>
	void each_node_object(Func func)
	{
		auto& storage = componentStorage<TNodeObject>();

		std::for_each(storage.begin(), storage.end(), [&, func = std::move(func)](TNodeObject& object) mutable
		{
			if (!object.isValid())
				return;

			FGraphNode& node = object.node(*this);

			if (node.isValid())
				func(node, object);
		});
	}

	template<typename TEdgeObject>
	EdgeView<TEdgeObject> edgeView()
	{
		return EdgeView<TEdgeObject>(rawEdgeStorage<TEdgeObject>(), privateEdges);
	}
};

template<class TComponent>
TComponent& FGraph::component(FGraphNodeHandle nodeHandle)
{
	return componentStorage<TComponent>().componentForNode(nodeHandle);
}

template<class TComponent>
TComponent* FGraph::componentPtr(FGraphNodeHandle nodeHandle)
{
	return componentStorage<TComponent>().componentPtrForNode(nodeHandle);
}


template<>
struct TStructOpsTypeTraits<FGraph> : public TStructOpsTypeTraitsBase2<FGraph>
{
	enum
	{
		WithPostSerialize = false,
	};
};

// ------------------------------------------------------------
// FGraphObject implementation
// ------------------------------------------------------------

template<class TComponent, class... TArgs >
inline TComponent& FGraphNode::addComponent( FGraph& graph, TArgs... args )
{
	TypedComponentStorage<TComponent>& storage = graph.componentStorage<TComponent>();


	ComponentType type = componentType<TComponent>();

	TComponent * comp = storage.componentPtrForNode(handle());

	// reinitialize the component
	if (comp)
	{
		*comp = TComponent(std::forward<TArgs>(args)...);

		comp->nodeIndex = id;
	}
	else
	{
		comp = &storage.emplace<TComponent>(*this, std::forward<TArgs>(args)...);
	}

	if (!components.Contains(type)) components.Add(type);

	graph.componentAdded(*this, type);

	return *comp;
}

inline FGraphObject* FGraphNode::addComponent( FGraph& graph, ComponentType type )
{
	FComponentStorage& storage = graph.componentStorage( type );

	UScriptStruct * typeStruct = FGraphObject::componentStruct( type );

	FGraphObject * comp = storage.componentForNode(handle());

	// reinitialize the component
	if (comp)
	{
		typeStruct->GetCppStructOps()->Destruct(comp);
		typeStruct->GetCppStructOps()->Construct(comp);

		// restore the node index
		comp->nodeIndex = id;

	}
	else
	{
		comp = storage.emplace(*this, typeStruct);
	}

	if (!components.Contains(type)) components.Add(type);

	graph.componentAdded( *this, type );

	return comp;
}

inline FGraphObject* FGraphNode::component( FGraph& graph, ComponentType type )
{
	return graph.componentStorage(type).componentForNode( *this, type );
}

template<class TComponent>
inline TComponent& FGraphNode::component(FGraph* graph)
{
	return graph->componentStorage<TComponent>().componentForNode(*this);
}

template<class TComponent>
inline TComponent& FGraphNode::component( FGraph& graph )
{
	return graph.componentStorage<TComponent>().componentForNode( *this );  
}

template<class TComponent>
inline void FGraphNode::removeComponent( FGraph& graph )
{
	ComponentType type = componentType<TComponent>();

    removeComponent( graph, type );
}

inline void FGraphNode::removeComponent( FGraph& graph, ComponentType type )
{
    FGraphObject * component = this->component( graph, type );

	component->invalidate();

	graph.componentRemoved( *this, type );

	components.Remove( type );
}

template<class TComponent>
inline bool FGraphNode::hasComponent()
{
	ComponentType type = componentType<TComponent>();

	return components.Contains( type );
}

inline bool FGraphNode::hasComponent( ComponentType type )
{
	return components.Contains( type );
}

template<class TEdgeObject, typename Func>
inline void FGraphNode::each(FGraph& graph, Func func)
{
	auto edgeObjects = graph.edgeStorage<TEdgeObject>();

	for (auto ei : edges)
	{
		FGraphEdgeHandle edgeHandle(ei);

		if (auto edgeObjectPtr = edgeObjects.objectPtr(edgeHandle))
		{
			FGraphEdge& edge = graph.edge(edgeHandle);

			auto other = edge.other(handle());

			if(edgeObjects.isValid(edgeHandle))
				func(other, *edgeObjectPtr);
		}
	}
}