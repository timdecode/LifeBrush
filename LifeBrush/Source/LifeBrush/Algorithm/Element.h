//
//  Element.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-06-25.
//  Copyright (c)  Timothy Davison. All rights reserved.
//

#pragma once

#include "Eigen/Dense"
#include "TimStructBox.h"
#include "NeighbourhoodParameters.h"
#include "SurfaceIndex.h"
#include <stdint.h>

#include "ShipEditorSimulation/Graph.h"
#include "Utility.h"

#include "Element.generated.h"

USTRUCT(BlueprintType)
struct FElementObject : public FGraphObject
{
	GENERATED_BODY()
public:


	UPROPERTY(EditAnywhere) float radius = 70.0f;
	UPROPERTY(EditAnywhere) int16 entityType = 0;		// entity type, cached from the parent entity for fast lookups
	UPROPERTY(EditAnywhere) int16 type = 0;			// element type
	UPROPERTY(EditAnywhere) bool generative = true;
	UPROPERTY(EditAnywhere) int32 entityIndex = 0;	// index of the parent entity

	UPROPERTY(EditAnywhere) FSurfaceIndex surfaceIndex; // if surface walking, the surface index that the element is on


	UPROPERTY(EditAnywhere) FNeighbourhoodParameters generationParameters;
	UPROPERTY(EditAnywhere) FNeighbourhoodParameters optimizationParameters;

	UPROPERTY(EditAnywhere) float minAssignmentDistance = 1.0f;
	UPROPERTY(EditAnywhere) float freespaceRadius = 1.0f;

	UPROPERTY(EditAnywhere) float generationInnerRadius = 1.0f;
};

struct ElementTuple
{
	ElementTuple(FGraphNodeHandle handle, FGraph& graph) :
		element(graph.component<FElementObject>(handle)),
		node(graph.node(handle))
	{}

	ElementTuple(FGraphNode& node, FGraph& graph) : 
		element(graph.component<FElementObject>(node.handle())),
		node(node)
	{}

	ElementTuple(FGraphNode& node, FElementObject& element) : 
		element(element),
		node(node)
	{}

	FElementObject& element;
	FGraphNode& node;

	Eigen::Vector3f position() { return eigen(node.position); }
	Eigen::Quaternionf rotation() { return eigen(node.orientation); }

	FGraphNodeHandle handle() { return node.handle(); }
};

namespace Element
{
	template<typename Container, typename Func>
	static void eachTuple(const Container& handles, FGraph& graph, Func func)
	{
		auto& storage = graph.componentStorage<FElementObject>();

		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &graph, &storage](FGraphNodeHandle h) mutable
		{
			FElementObject * element = storage.componentPtrForNode(h);
			FGraphNode& node = h(graph);

			if (element)
				func(ElementTuple(node, *element));
		});
	}

	template<typename Container, typename Func>
	static void eachTuple(Container& handles, FGraph& graph, Func func)
	{
		auto& storage = graph.componentStorage<FElementObject>();

		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &graph, &storage](FGraphNodeHandle h) mutable
		{
			FElementObject * element = storage.componentPtrForNode(h);
			FGraphNode& node = h(graph);

			if (element)
				func(ElementTuple(node, *element));
		});
	}

	template<typename Container, typename Func>
	static void eachNode(const Container& handles, FGraph& graph, Func func)
	{
		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &graph](FGraphNodeHandle h) mutable
		{
			func(graph.node(h));
		});
	}

	template<typename Container, typename Func>
	static void eachNode(Container& handles, FGraph& graph, Func func)
	{
		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &graph](FGraphNodeHandle h) mutable
		{
			func(graph.node(h));
		});
	}

	template<typename Container, typename Func>
	static void eachElement(const Container& handles, FGraph& graph, Func func)
	{
		auto& storage = graph.componentStorage<FElementObject>();

		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &storage](FGraphNodeHandle h) mutable
		{
			FElementObject * element = storage.componentPtrForNode(h);

			if(element)
				func(*element);
		});
	}

	template<typename Container, typename Func>
	static void eachElement(Container& handles, FGraph& graph, Func func)
	{
		auto& storage = graph.componentStorage<FElementObject>();

		std::for_each(handles.begin(), handles.end(), [func = std::move(func), &storage](FGraphNodeHandle h) mutable
		{
			FElementObject * element = storage.componentPtrForNode(h);

			if (element)
				func(*element);
		});
	}
};

USTRUCT(BlueprintType)
struct FElement
{
	GENERATED_BODY()
public:

    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    Eigen::Quaternionf orientation = Eigen::Quaternionf::Identity();
    
	UPROPERTY( EditAnywhere ) float radius = 70.0f;
	UPROPERTY( EditAnywhere ) int16 entityType = 0;		// entity type, cached from the parent entity for fast lookups
	UPROPERTY( EditAnywhere ) int16 type = 0;			// element type
	UPROPERTY( EditAnywhere ) bool generative = true;
	UPROPERTY( EditAnywhere ) int32 entityIndex = 0;	// index of the parent entity

	UPROPERTY(EditAnywhere) FSurfaceIndex surfaceIndex; // if surface walking, the surface index that the element is on


	UPROPERTY( EditAnywhere ) FNeighbourhoodParameters generationParameters;
	UPROPERTY( EditAnywhere ) FNeighbourhoodParameters optimizationParameters;

	UPROPERTY( EditAnywhere ) float minAssignmentDistance = 1.0f;
	UPROPERTY( EditAnywhere ) float freespaceRadius = 1.0f;

	UPROPERTY( EditAnywhere ) float generationInnerRadius = 1.0f;

	UPROPERTY( EditAnywhere ) float hackScale = 1.0f;

	UPROPERTY( EditAnywhere ) TArray<FTimStructBox> graphObjects;
    
    //bool operator==(const Element& other) const
    //{
    //    return position == other.position &&
    //           rotation == other.rotation &&
    //           radius == other.radius &&
    //           entityType == other.entityType &&
    //           type == other.type &&
    //           generative == other.generative &&
    //           entityIndex == other.entityIndex &&
    //           faceIndex == other.faceIndex &&
    //           gradient == other.gradient &&
    //           dotProduct == other.dotProduct;
    //}
};
