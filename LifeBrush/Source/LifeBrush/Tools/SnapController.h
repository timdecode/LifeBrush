// Copyright 2016 Timothy Davison, all rights reserved.
// This was written for the game Ship Editor.

#pragma once

#include "SnapController.generated.h"

/**
 * 
 */
USTRUCT( BlueprintType )
struct LIFEBRUSH_API FSnapController
{
	GENERATED_BODY()

public:
	FSnapController() {}

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	bool snapPosition = false;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float positionSpacing = 5.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	bool snapRotation = false;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float rotationSpacing = 11.25f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float shipScale = 1.0f;

	FQuat snappedRotation( FQuat quat, bool forceSnap = false )
	{
		if(!snapRotation && !forceSnap)
			return quat;
		else
		{
			FRotator gridSettings( rotationSpacing, rotationSpacing, rotationSpacing );

			return quat.Rotator().GridSnap( gridSettings ).Quaternion();
		}
	}

	FVector snappedPosition( FVector position, bool forceSnap = false )
	{
		if(!snapPosition && !forceSnap)
			return position;
		else
			return position.GridSnap( positionSpacing / shipScale );
	}
};