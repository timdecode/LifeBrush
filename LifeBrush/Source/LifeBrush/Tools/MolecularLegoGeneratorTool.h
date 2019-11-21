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

	virtual void oneHandStart(UPrimitiveComponent * hand) override;
	virtual void oneHandEnd(UPrimitiveComponent * hand) override;

	virtual void tickOneHand(float dt, UPrimitiveComponent * hand, FTransform lastTransform) override;

	virtual void loseFocus() override;

protected:
	TArray<AElementActor*> _overlappingPrototypes(UPrimitiveComponent * hand, float radius);

	void _tickSelection(float dt, UPrimitiveComponent * hand, FTransform lastTransform);

	void _hideSelection();
	void _showSelection();

protected:
	enum class Mode {
		Selecting,
		Painting
	};

	Mode _mode = Mode::Selecting;

	AElementActor * _selection = nullptr;
};
