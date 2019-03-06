// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#if WITH_EDITOR

#include "Editor/PropertyEditor/Public/IDetailCustomization.h"

#include "SlateBasics.h"
#include "SlateExtras.h"
#include "LevelEditor.h"

#include "Editor/PropertyEditor/Public/PropertyEditorModule.h"
#include "Editor/PropertyEditor/Public/PropertyEditing.h"

/**
 * 
 */

class  URegionGrowingComponentEditor : public IDetailCustomization
{
public:
    void LayoutDetails( IDetailLayoutBuilder& );
    
    static TSharedRef<IDetailCustomization> MakeInstance();
    
    virtual void CustomizeDetails(IDetailLayoutBuilder& DetailLayout) override;

};

#endif // WITH_EDITOR