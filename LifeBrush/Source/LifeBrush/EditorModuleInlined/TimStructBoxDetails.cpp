// Copyright 2018, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#if WITH_EDITOR

#include "TimStructBoxDetails.h"

#include "TimStructBox.h"
#include "ShipEditorSimulation/MeshSimulation.h"

#define LOCTEXT_NAMESPACE "FTimStructBoxCustomization"

TSharedRef<IPropertyTypeCustomization> FTimStructBoxCustomization::MakeInstance()
{
	return MakeShareable( new FTimStructBoxCustomization );
}




void FTimStructBoxCustomization::CustomizeHeader( TSharedRef<class IPropertyHandle> StructPropertyHandle, FDetailWidgetRow & HeaderRow, IPropertyTypeCustomizationUtils & StructCustomizationUtils )
{

}

void FTimStructBoxCustomization::OnChildValueChanged()
{
	TArray<void*> rawData;

	_structPropertyHandle->AccessRawData( rawData );

	FTimStructBox& baseStruct = *reinterpret_cast<FTimStructBox*>(rawData.Top());

	if(baseStruct.IsValid())
	{
		for(void * data : rawData)
		{
			if(data == &baseStruct)
				continue;

			FTimStructBox& otherStruct = *reinterpret_cast<FTimStructBox*>(data);

			otherStruct = baseStruct;
		}
	}
}

void FTimStructBoxCustomization::CustomizeChildren(
	TSharedRef<class IPropertyHandle> structPropertyHandle,
	class IDetailChildrenBuilder& structBuilder,
	IPropertyTypeCustomizationUtils& structCustomizationUtils 
)
{
	_structPropertyHandle = structPropertyHandle;

	TArray<void*> rawData;

	structPropertyHandle->AccessRawData( rawData );

	FTimStructBox& templateStruct = *reinterpret_cast<FTimStructBox*>(rawData.Top());

	if(templateStruct.IsValid())
	{
		UScriptStruct * boxedScriptStruct = templateStruct.scriptStruct;
		TSharedRef<FStructOnScope> structData( new FStructOnScope( boxedScriptStruct, templateStruct.structMemory ) );

		IDetailPropertyRow * row = structBuilder.AddExternalStructure( structData );

		auto externalPropertyHandle = row->GetPropertyHandle();
		
		auto valueChangedDelegate = FSimpleDelegate::CreateSP( this, &FTimStructBoxCustomization::OnChildValueChanged );
		externalPropertyHandle->SetOnChildPropertyValueChanged( valueChangedDelegate );
	}
	else
	{
		IDetailCategoryBuilder * categoryBuilder = &structBuilder.GetParentCategory();

		// populate our menu delegate
		TSharedPtr< const FUICommandList > commandsEmpty;
		FMenuBuilder componentTypesBuilder(true, commandsEmpty);
		{
			TArray<UScriptStruct*> componentStructs;
			for(TObjectIterator<UScriptStruct> it; it; ++it)
			{
				if(it->IsChildOf( FGraphObject::StaticStruct() ))
					componentStructs.Add( *it );
			}

			for(UScriptStruct * componentStruct : componentStructs)
			{
				// ----------------------
				// Create a component
				// ----------------------
				FUIAction action( FExecuteAction::CreateLambda( [rawData, componentStruct, categoryBuilder]() {
					for(void * data : rawData)
					{
						FTimStructBox& structBox = *reinterpret_cast<FTimStructBox*>(data);

						structBox.scriptStruct = componentStruct;

						structBox.initStruct();
					}

					categoryBuilder->GetParentLayout().ForceRefreshDetails();
				} ));

				componentTypesBuilder.AddMenuEntry(
					FText::FromName( componentStruct->GetFName() ),
					FText::FromName( componentStruct->GetFName() ),
					FSlateIcon(),
					action,
					NAME_None,
					EUserInterfaceActionType::Button 
				);
			}
		}

		// build the menu widget
		structBuilder.AddCustomRow( LOCTEXT( "Class", "Class" ) )
		.NameContent()[
			SNew( STextBlock )
			.Text( LOCTEXT( "Init FGraphNode class", "Init FGraphNode clas" ) )
		]
		.ValueContent()[
			SNew( SComboButton )
			.ButtonContent()
			[
				SNew( STextBlock )
				.Text( LOCTEXT("Class", "Class") )
			]
			.MenuContent()
			[
				componentTypesBuilder.MakeWidget()
			]
		];
	}
}

#endif // WITH_EDITOR

#undef LOCTEXT_NAMESPACE
