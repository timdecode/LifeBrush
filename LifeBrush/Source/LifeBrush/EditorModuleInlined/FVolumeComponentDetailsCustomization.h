// Copyright 2017 Code Monkey Castle, all rights reserved.

#pragma once

#if WITH_EDITOR

#include "Editor/PropertyEditor/Public/IDetailCustomization.h"


#include "IPropertyTypeCustomization.h"
#include "IDetailCustomization.h"

class FVolumeComponentDetailsCustomization : public IDetailCustomization
{
public:
	static TSharedRef<IDetailCustomization> MakeInstance();

	/** IDetailCustomization interface */
	virtual void CustomizeDetails( IDetailLayoutBuilder& DetailBuilder ) override;
};

class FChunkedVolumeComponentDetailsCustomization : public IDetailCustomization
{
public:
	static TSharedRef<IDetailCustomization> MakeInstance();

	/** IDetailCustomization interface */
	virtual void CustomizeDetails(IDetailLayoutBuilder& DetailBuilder) override;
};

#endif // WITH_EDITOR