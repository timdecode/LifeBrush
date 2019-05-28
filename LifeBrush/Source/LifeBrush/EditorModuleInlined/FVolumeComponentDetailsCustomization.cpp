// Copyright 2017 Code Monkey Castle, all rights reserved.

#include "LifeBrush.h"

#if WITH_EDITOR

#include "RuntimeMeshComponent.h"
#include "FVolumeComponentDetailsCustomization.h"
#include "PropertyEditing.h"


#include "VolumeComponent.h"

#define LOCTEXT_NAMESPACE "FVolumeComponentDetailsCustomization"


TSharedRef<IDetailCustomization> FVolumeComponentDetailsCustomization::MakeInstance()
{
	return MakeShareable( new FVolumeComponentDetailsCustomization );
}

void FVolumeComponentDetailsCustomization::CustomizeDetails( IDetailLayoutBuilder& detailBuilder )
{
	IDetailCategoryBuilder& category = detailBuilder.EditCategory( "Functions", LOCTEXT( "Functions", "Functions" ), ECategoryPriority::Important );

	category.AddCustomRow( LOCTEXT( "Build Geometry", "Build Geometry" ) )
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Build Geometry", "Build Geometry" ) )
		.OnClicked_Lambda( [&detailBuilder]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailBuilder.GetObjectsBeingCustomized( outObjects );

			UVolumeComponent * blocksComponent = Cast<UVolumeComponent>( outObjects.Last().Get() );
			blocksComponent->build();

			return FReply::Handled();
		} )
	];
}

TSharedRef<IDetailCustomization> FChunkedVolumeComponentDetailsCustomization::MakeInstance()
{
	return MakeShareable(new FChunkedVolumeComponentDetailsCustomization);
}

void FChunkedVolumeComponentDetailsCustomization::CustomizeDetails(IDetailLayoutBuilder& detailBuilder)
{
	IDetailCategoryBuilder& category = detailBuilder.EditCategory("Functions", LOCTEXT("Functions", "Functions"), ECategoryPriority::Important);

	category.AddCustomRow( LOCTEXT( "Build Geometry", "Build Geometry" ) )
	.ValueContent()
	[
		SNew( SButton )
		.Text( LOCTEXT( "Initialize test volume", "Initialize test volume" ) )
		.OnClicked_Lambda( [&detailBuilder]()->FReply {
			TArray<TWeakObjectPtr<UObject> > outObjects;
			detailBuilder.GetObjectsBeingCustomized( outObjects );

			UChunkedVolumeComponent * volumeComponent = Cast<UChunkedVolumeComponent>( outObjects.Last().Get() );

			volumeComponent->initializeTestVolume();

			volumeComponent->writeAccessGrid([&](ChunkGrid<float>& grid) 
			{
				FIntVector min = grid.gridIndexFromChunkIndex(grid.minChunkndex());
				FIntVector max = grid.gridIndexFromChunkIndex(grid.maxChunkIndex()) + grid.chunkDimensions();

				volumeComponent->markDirtyIndex(min, max);
			});

			return FReply::Handled();
		} )
	];
}

#undef LOCTEXT_NAMESPACE

#endif // WITH_EDITOR