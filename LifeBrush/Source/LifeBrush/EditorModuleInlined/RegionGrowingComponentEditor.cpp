// Fill out your copyright notice in the Description page of Project Settings.

#include "LifeBrush.h"

#if WITH_EDITOR

#include "RegionGrowingComponentEditor.h"
#include "LifeBrush/RegionGrowingComponent.h"

#include "Editor/PropertyEditor/Public/IDetailCustomization.h"

#define LOCTEXT_NAMESPACE "URegionGrowingComponentEditor"

void URegionGrowingComponentEditor::LayoutDetails( IDetailLayoutBuilder& )
{
    
}

TSharedRef<IDetailCustomization> URegionGrowingComponentEditor::MakeInstance()
{
    UE_LOG(LogTemp, Log, TEXT("MakeInstance") );
    
    
    return MakeShareable(new URegionGrowingComponentEditor);
}


void URegionGrowingComponentEditor::CustomizeDetails(IDetailLayoutBuilder& detailLayout)
{
    IDetailCategoryBuilder& category = detailLayout.EditCategory("Functions", LOCTEXT("Functions", "Functions"), ECategoryPriority::Important);
    
    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
    .ValueContent()
    [
        SNew(SButton)
        .Text(LOCTEXT("Load Exemplar", "Load Exemplar"))
        .OnClicked_Lambda([&detailLayout]()->FReply{
            TArray<TWeakObjectPtr<UObject> > outObjects;
            detailLayout.GetObjectsBeingCustomized(outObjects);
        
            URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>(outObjects.Last().Get());
            regionGrowingComponent->LoadExemplar();
        
            return FReply::Handled();
        })
    ];

    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
    .ValueContent()
    [
        SNew(SButton)
        .Text(LOCTEXT("Generate", "Generate"))
        .OnClicked_Lambda([&detailLayout]()->FReply{
            TArray<TWeakObjectPtr<UObject> > outObjects;
            detailLayout.GetObjectsBeingCustomized(outObjects);
        
            URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>(outObjects.Last().Get());

			regionGrowingComponent->Generate();
        
            return FReply::Handled();
        })
    ];
    
    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
    .ValueContent()
    [
        SNew(SButton)
        .Text(LOCTEXT("Global Optimization", "Global Optimization"))
        .OnClicked_Lambda([&detailLayout]()->FReply{
            TArray<TWeakObjectPtr<UObject> > outObjects;
            detailLayout.GetObjectsBeingCustomized(outObjects);
        
            URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>(outObjects.Last().Get());
            regionGrowingComponent->GlobalOptimization();
        
            return FReply::Handled();
        })
    ];

    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
    .ValueContent()
    [
        SNew(SButton)
        .Text(LOCTEXT("Load Parameters", "Load Parameters"))
        .OnClicked_Lambda([&detailLayout]()->FReply{
            TArray<TWeakObjectPtr<UObject> > outObjects;
            detailLayout.GetObjectsBeingCustomized(outObjects);
        
            URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>(outObjects.Last().Get());
            regionGrowingComponent->loadParameters();
        
            return FReply::Handled();
        })
    ];

    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
    .ValueContent()
    [
        SNew(SButton)
        .Text(LOCTEXT("Clear Elements Keep Paths", "Clear Elements Keep Paths"))
        .OnClicked_Lambda([&detailLayout]()->FReply{
            TArray<TWeakObjectPtr<UObject> > outObjects;
            detailLayout.GetObjectsBeingCustomized(outObjects);
        
            URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>(outObjects.Last().Get());
            regionGrowingComponent->ClearElementsKeepPaths();
        
            return FReply::Handled();
        })
    ];
    
    category.AddCustomRow(LOCTEXT("A Function", "A header name"))
    .NameContent()
    [
        SNew(STextBlock)
        .Text(LOCTEXT("FunctionHeader", ""))
        .Font(IDetailLayoutBuilder::GetDetailFont())
    ]
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Export ISMC Actor", "Expoort ISMC Actor" ) )
		.OnClicked_Lambda( [&detailLayout]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailLayout.GetObjectsBeingCustomized( outObjects );

			UWorld * world = GEditor->EditorWorld;

			URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>( outObjects.Last().Get() );
			regionGrowingComponent->createGraphicalSnapshotActor(world);

			return FReply::Handled();
		} )
	];

	category.AddCustomRow( LOCTEXT( "A Function", "A header name" ) )
	.NameContent()
	[
		SNew( STextBlock )
		.Text( LOCTEXT( "FunctionHeader", "" ) )
		.Font( IDetailLayoutBuilder::GetDetailFont() )
	]
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Output Stats", "Output Stats" ) )
		.OnClicked_Lambda( [&detailLayout]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailLayout.GetObjectsBeingCustomized( outObjects );

			URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>( outObjects.Last().Get() );
			regionGrowingComponent->OutputStats();

			return FReply::Handled();
		} )
	];

	category.AddCustomRow( LOCTEXT( "A Function", "A header name" ) )
	.NameContent()
	[
		SNew( STextBlock )
		.Text( LOCTEXT( "FunctionHeader", "" ) )
		.Font( IDetailLayoutBuilder::GetDetailFont() )
	]
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Tab Separated Value Stats", "Tab Separated Value Stats" ) )
		.OnClicked_Lambda( [&detailLayout]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailLayout.GetObjectsBeingCustomized( outObjects );

			URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>( outObjects.Last().Get() );
			regionGrowingComponent->OutputSpaceSeparatedStats();

			return FReply::Handled();
		} )
	];
    
	category.AddCustomRow( LOCTEXT( "A Function", "A header name" ) )
	.NameContent()
	[
		SNew( STextBlock )
		.Text( LOCTEXT( "FunctionHeader", "" ) )
		.Font( IDetailLayoutBuilder::GetDetailFont() )
	]
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Log Parameters", "Log Parameters" ) )
		.OnClicked_Lambda( [&detailLayout]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailLayout.GetObjectsBeingCustomized( outObjects );

			URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>( outObjects.Last().Get() );
			regionGrowingComponent->OutputParameters();

			return FReply::Handled();
		} )
	];

	category.AddCustomRow( LOCTEXT( "A Function", "A header name" ) )
	.NameContent()
	[
		SNew( STextBlock )
		.Text( LOCTEXT( "FunctionHeader", "" ) )
		.Font( IDetailLayoutBuilder::GetDetailFont() )
	]
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Trim to trimBounds", "Trim to trimBounds" ) )
		.OnClicked_Lambda( [&detailLayout]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailLayout.GetObjectsBeingCustomized( outObjects );

			URegionGrowingComponent * regionGrowingComponent = Cast<URegionGrowingComponent>( outObjects.Last().Get() );
			regionGrowingComponent->Trim();

			return FReply::Handled();
		} )
	];
}

#undef LOCTEXT_NAMESPACE

#endif // WITH_EDITOR
