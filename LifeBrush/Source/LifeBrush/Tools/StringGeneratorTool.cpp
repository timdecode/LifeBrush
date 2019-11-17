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

