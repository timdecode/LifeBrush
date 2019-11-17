// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"
#include "ShipEditorSimulation/ObjectSimulation.h"

#include "Actin.generated.h"

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLActin : public FGraphObject
{
	GENERATED_BODY()

public:
	static constexpr std::size_t _bondsSize = 4;
	std::array<FGraphNodeHandle, _bondsSize> _bonds;
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UMLActinSimulation : public UObjectSimulation
{
	GENERATED_BODY()

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void tick(float deltaT) override;

	void _tickRules();
	void _tickRules2();

protected:
	void _cacheBonds();


protected:
	bool _bondsDirty = true;
};