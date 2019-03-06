// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"
#include "GraphUtility.h"

// Copied (with modification) from LevelStreaming.cpp
ULevelStreaming * ShipEditorUtility::createInstance( ULevelStreaming * InLevelStreaming, FString InstanceUniqueName )
{
	ULevelStreaming* StreamingLevelInstance = nullptr;

	UWorld* InWorld = InLevelStreaming->GetWorld();
	if(InWorld)
	{
		// Create instance long package name 
		FString InstanceShortPackageName = InWorld->StreamingLevelsPrefix + FPackageName::GetShortName( InstanceUniqueName );
		FString InstancePackagePath = FPackageName::GetLongPackagePath( InLevelStreaming->GetWorldAssetPackageName() ) + TEXT( "/" );
		FName	InstanceUniquePackageName = FName( *(InstancePackagePath + InstanceShortPackageName) );

		// check if instance name is unique among existing streaming level objects
		const bool bUniqueName = (InWorld->StreamingLevels.IndexOfByPredicate( ULevelStreaming::FPackageNameMatcher( InstanceUniquePackageName ) ) == INDEX_NONE);

		if(bUniqueName)
		{
			StreamingLevelInstance = NewObject<ULevelStreaming>( InWorld, InLevelStreaming->GetClass(), NAME_None, RF_Transient, NULL );
			
			// new level streaming instance will load the same map package as this object
			StreamingLevelInstance->PackageNameToLoad = (InLevelStreaming->PackageNameToLoad == NAME_None ? InLevelStreaming->GetWorldAssetPackageFName() : InLevelStreaming->PackageNameToLoad);
			
			// under a provided unique name
			StreamingLevelInstance->SetWorldAssetByPackageName( InstanceUniquePackageName );
			StreamingLevelInstance->bShouldBeLoaded = false;
			StreamingLevelInstance->bShouldBeVisible = false;
			StreamingLevelInstance->LevelTransform = InLevelStreaming->LevelTransform;

			// add a new instance to streaming level list
			InWorld->StreamingLevels.Add( StreamingLevelInstance );
		}
		else
		{
			UE_LOG( LogStreaming, Warning, TEXT( "Provided streaming level instance name is not unique: %s" ), *InstanceUniquePackageName.ToString() );
		}
	}

	return StreamingLevelInstance;
}
