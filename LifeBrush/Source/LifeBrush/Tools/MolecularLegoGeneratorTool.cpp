//
//  Created by Timothy Davison on 2019-10-14
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"

#include "ElementEditor/MolecularLegoGenerator.h"

#include "MolecularLegoGeneratorTool.h"

void UMolecularLegoGeneratorTool::init(FRGC_UToolInitProperties& initProperties)
{
	Super::init(initProperties);

	if (!elementEditor)
		return;

	molecularGenerator = elementEditor->generator<UMolecularLegoGenerator>();

	if (!molecularGenerator)
		molecularGenerator = NewObject<UMolecularLegoGenerator>(this, TEXT("MolecularLegoGenerator"));

	_generator = molecularGenerator;
}

void UMolecularLegoGeneratorTool::faceUp_released(USceneComponent * interactionPoint /*= nullptr*/)
{
	if (molecularGenerator->mode == EMolecularLegotGeneratorMode::AlongPath)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::Grow;
	else if (molecularGenerator->mode == EMolecularLegotGeneratorMode::Grow)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::PlaceOne;
	else if (molecularGenerator->mode == EMolecularLegotGeneratorMode::PlaceOne)
		molecularGenerator->mode = EMolecularLegotGeneratorMode::AlongPath;
}
