// Copyright (c) 2018 Timothy Davison. All rights reserved.

#pragma once

#include "Domain.h"
#include "tcodsMeshInterface.h"
#include "ctpl_stl.h"
#include "ShipEditorSimulation/Graph.h"

struct SynthesisContext
{
public:
	SynthesisContext(FGraph& graph_in) : domain(graph_in)
	{
		unsigned int nThreads = std::thread::hardware_concurrency();

		threadPool = std::make_shared<ctpl::thread_pool>();
		threadPool->resize(nThreads);
	}

	Domain domain;

	std::shared_ptr<tcodsMeshInterfaceBase> meshInterface;

	std::shared_ptr<ctpl::thread_pool> threadPool;

	FBox limits;

	inline FGraph& graph() { return domain.graph; }
};