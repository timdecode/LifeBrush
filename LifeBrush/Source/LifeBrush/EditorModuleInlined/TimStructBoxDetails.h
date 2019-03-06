// Copyright 2018, Timothy Davison. All rights reserved.

#pragma once

#if WITH_EDITOR

#include "IPropertyTypeCustomization.h"
#include "IDetailCustomization.h"

class IPropertyHandle;

class FTimStructBoxCustomization : public IPropertyTypeCustomization
{
public:

	static TSharedRef<IPropertyTypeCustomization> MakeInstance();

	virtual void CustomizeHeader( TSharedRef<class IPropertyHandle> StructPropertyHandle, class FDetailWidgetRow& HeaderRow, IPropertyTypeCustomizationUtils& StructCustomizationUtils ) override;
	virtual void CustomizeChildren( TSharedRef<class IPropertyHandle> StructPropertyHandle, class IDetailChildrenBuilder& StructBuilder, IPropertyTypeCustomizationUtils& StructCustomizationUtils ) override;

protected:
	void OnChildValueChanged();

	TSharedPtr<IPropertyHandle> _structPropertyHandle;
};

#endif // WITH_EDITOR