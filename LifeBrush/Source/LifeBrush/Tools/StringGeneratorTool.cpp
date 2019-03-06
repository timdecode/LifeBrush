//
//  Created by Timothy Davison on 2018-12-28.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementEditor/SwarmGenerator.h"

#include "StringGeneratorTool.h"

void UStringGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	UStringGenerator * stringGenerator = elementEditor->generator<UStringGenerator>();

	if (!stringGenerator)
		stringGenerator = NewObject<UStringGenerator>(this, TEXT("stringGenerator"));

	_generator = stringGenerator;
}

void UCollagenGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	UCollagenGenerator * generator = elementEditor->generator<UCollagenGenerator>();

	if (!generator)
		generator = NewObject<UCollagenGenerator>(this, TEXT("collagenGenerator"));

	_generator = generator;
}

void USwarmGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	USwarmGenerator * generator = elementEditor->generator<USwarmGenerator>();

	if (!generator)
		generator = NewObject<USwarmGenerator>(this, TEXT("swarmGenerator"));

	_generator = generator;
}

USwarmGenerator * USwarmGeneratorTool::generator()
{
	return Cast<USwarmGenerator>(_generator);
}

void USwarmGeneratorTool::faceUp_released(USceneComponent * interactionPoint /* = nullptr */)
{
	auto& brushType = generator()->brushType;

	if (brushType == ESwarmGenerator_BrushType::Anchor)
		brushType = ESwarmGenerator_BrushType::Star;
	else
		brushType = ESwarmGenerator_BrushType::Anchor;
}
