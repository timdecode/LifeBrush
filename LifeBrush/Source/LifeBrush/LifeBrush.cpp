// Fill out your copyright notice in the Description page of Project Settings.

#include "LifeBrush.h"


#ifdef __APPLE__
#include <dlfcn.h>
#include <cstring>
#include <string>
#endif

#include "ShipEditorSimulation/GraphVersion.h"

#if WITH_EDITOR
#include "EditorModuleInlined/RegionGrowingComponentEditor.h" 
#include "EditorModuleInlined/TCODSEditorMode.h"
#include "EditorModuleInlined/FVolumeComponentDetailsCustomization.h"
#include "EditorModuleInlined/TimStructBoxDetails.h"
#include "EditorModuleInlined/GraphPropertyCustomization.h"
#include "EditorModuleInlined/GraphSimulationEdMode.h"
#endif

IMPLEMENT_PRIMARY_GAME_MODULE( FLifeBrush, LifeBrush, "LifeBrush" );

FCustomVersionRegistration GRegisterGraphCustomVersion(FGraphVersion::GUID, FGraphVersion::LatestVersion, TEXT("FGraphVersion latest"));


static void anAdressToFindTheDll() {}

void FLifeBrush::StartupModule()
{
#ifdef __APPLE__
    Dl_info info;
    if( dladdr((void*)&anAdressToFindTheDll, &info) )
    {
        
        const size_t cSize = strlen(info.dli_fname);
        std::wstring ws(cSize, L'#');
        mbstowcs(&ws[0], info.dli_fname, cSize);
        FString dylibName(cSize, &ws[0]);
        
        UE_LOG(LogTemp, Warning, TEXT("FRegionGrowing::StartModule() loaded from %s"), *dylibName);
    }
#endif

#if WITH_EDITOR
    registerCustomClassLayout( "RegionGrowingComponent", FOnGetDetailCustomizationInstance::CreateStatic(&URegionGrowingComponentEditor::MakeInstance));
	registerCustomClassLayout( "VolumeComponent", FOnGetDetailCustomizationInstance::CreateStatic( &FVolumeComponentDetailsCustomization::MakeInstance ) );
	registerCustomClassLayout( "ChunkedVolumeComponent", FOnGetDetailCustomizationInstance::CreateStatic(&FChunkedVolumeComponentDetailsCustomization::MakeInstance));


    // Unreal doesn't want to Hot Reload classes correctly. We often seem to end up with two different versions of a library loaded in memory.
    FPropertyEditorModule& propertyModule = FModuleManager::GetModuleChecked<FPropertyEditorModule>("PropertyEditor");

	FName structBoxName = FTimStructBox::StaticStruct()->GetFName();
	FName graphName = FGraph::StaticStruct()->GetFName();
	FName graphNodeName = FGraphNode::StaticStruct()->GetFName();

	propertyModule.RegisterCustomPropertyTypeLayout( structBoxName, FOnGetPropertyTypeCustomizationInstance::CreateStatic( &FTimStructBoxCustomization::MakeInstance ) );
	propertyModule.RegisterCustomPropertyTypeLayout( graphName, FOnGetPropertyTypeCustomizationInstance::CreateStatic(&FGraphPropteryCustomization::MakeInstance));
	propertyModule.RegisterCustomPropertyTypeLayout( graphNodeName, FOnGetPropertyTypeCustomizationInstance::CreateStatic(&FGraphNodePropteryCustomization::MakeInstance));

	propertyModule.NotifyCustomizationModuleChanged();

	FEditorModeRegistry::Get().RegisterMode<FTCODSEditorMode>( FTCODSEditorMode::EM_Geometry,
		NSLOCTEXT( "EditorModes", "TrivialConnectionsMode", "Trivial Connections Mesh Editing" ),
		FSlateIcon( FEditorStyle::GetStyleSetName(), "LevelEditor.BspMode", "LevelEditor.BspMode.Small" ),
		true, 500);

	FEditorModeRegistry::Get().RegisterMode<FGraphSimulationEdMode>(FGraphSimulationEdMode::EM_GraphSimulationEdModeId,
		NSLOCTEXT("EditorModes", "UGraphComponent editor", "UGraphComponent editor"),
		FSlateIcon(FEditorStyle::GetStyleSetName(), "LevelEditor.BspMode", "LevelEditor.BspMode.Small"),
		true, 500);

#endif
}

void FLifeBrush::ShutdownModule()
{
#if WITH_EDITOR
    FPropertyEditorModule& propertyModule = FModuleManager::GetModuleChecked<FPropertyEditorModule>("PropertyEditor");
    
    for( auto name : registeredClassNames )
        propertyModule.UnregisterCustomClassLayout(name);
    
	FName structBoxName = FTimStructBox::StaticStruct()->GetFName();
	FName graphName = FGraph::StaticStruct()->GetFName();
	FName graphNodeName = FGraphNode::StaticStruct()->GetFName();

	propertyModule.UnregisterCustomPropertyTypeLayout( structBoxName );
	propertyModule.UnregisterCustomPropertyTypeLayout( graphName );
	propertyModule.UnregisterCustomPropertyTypeLayout( graphNodeName );

    propertyModule.NotifyCustomizationModuleChanged();
    
    FEditorModeRegistry::Get().UnregisterMode(FTCODSEditorMode::EM_Geometry);
	FEditorModeRegistry::Get().UnregisterMode(FGraphSimulationEdMode::EM_GraphSimulationEdModeId);
#endif
}

#if WITH_EDITOR
void FLifeBrush::registerCustomClassLayout(FName className, FOnGetDetailCustomizationInstance detailLayoutDelegate)
{
    registeredClassNames.Add( className );
    
    static FName PropertyEditor("PropertyEditor");
    FPropertyEditorModule& propertyModule = FModuleManager::GetModuleChecked<FPropertyEditorModule>(PropertyEditor);
    propertyModule.RegisterCustomClassLayout( className, detailLayoutDelegate );
}
#endif
