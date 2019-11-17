// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "Engine.h"

#if WITH_EDITOR
#include "SlateBasics.h"
#include "SlateExtras.h"
#include "LevelEditor.h"

#include "Editor/PropertyEditor/Public/PropertyEditorModule.h"
#include "Editor/PropertyEditor/Public/PropertyEditing.h"
#endif

class FLifeBrush : public IModuleInterface
{
public:
    virtual bool IsGameModule() const override
    {
        return true;
    }
    
    virtual void StartupModule() override;
    virtual void ShutdownModule() override;
    
private:
#if WITH_EDITOR
    void registerCustomClassLayout(FName className, FOnGetDetailCustomizationInstance detailLayoutDelegate);


#endif
    
private:
    TSet<FName> registeredClassNames;
 };


