//
//  Created by Timothy Davison on 2019-10-14.
//  Copyright (c) 2019 Timothy Davison. All rights reserved.
//

#pragma once

#include "RegionGrowingGeneratorTool.h"

#include "ElementEditor/StringGenerator.h"
#include "ElementEditor/MolecularLegoGenerator.h"

#include "MolecularLegoGeneratorTool.generated.h"

class UMolecularLegoGenerator;

UCLASS(Blueprintable)
class UMolecularLegoGeneratorTool : public UGeneratorTool
{
	GENERATED_BODY()

protected:
	UMolecularLegoGenerator * molecularGenerator = nullptr;

public:
	virtual ~UMolecularLegoGeneratorTool() {}

	void init(FRGC_UToolInitProperties& initProperties);

	virtual void faceUp_released(USceneComponent * interactionPoint = nullptr) override;

};
