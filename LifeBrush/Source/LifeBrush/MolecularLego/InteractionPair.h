// Copyright (c) 2019 Timothy Davison. All rights reserved.
#pragma once

#include "ShipEditorSimulation/Graph.h"

struct InteractionPair
{
public:
	InteractionPair(FGraphNodeHandle a, FGraphNodeHandle b) : a(a), b(b) {}

	FGraphNodeHandle a;
	FGraphNodeHandle b;

	bool operator== (const InteractionPair& Other) const
	{
		return (a == Other.a && b == Other.b) || (a == Other.b && b == Other.a);
	}
};

FORCEINLINE uint32 GetTypeHash (const InteractionPair& s)
{
	const auto a_ = std::min(s.a.index, s.b.index);
	const auto b_ = std::max(s.a.index, s.b.index);

	auto ha = GetTypeHash(a_);
	auto hb = GetTypeHash(b_);

	return HashCombine(ha, hb);
}