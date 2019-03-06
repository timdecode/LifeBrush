//
//  SurfaceIndex.h
//  RegionGrowing
//
//  Created by Timothy Davison on 2018-10-18.
//  Copyright (c) 2018 Timothy Davison. All rights reserved.
//

#pragma once


#include "SurfaceIndex.generated.h"


USTRUCT(BlueprintType)
struct FSurfaceIndex
{
	GENERATED_BODY()
public:
	FSurfaceIndex() {}
	FSurfaceIndex(int32 sectionIndex, int32 faceIndex) : sectionIndex(sectionIndex), faceIndex(faceIndex) {}

	UPROPERTY(EditAnywhere) int32 faceIndex = -1;	  // if surface walking, the index of the face that the element is on
	UPROPERTY(EditAnywhere) int32 sectionIndex = -1;  // if surface walking, the mesh section index of the face that the element is on

	bool isOnSurface() const { return faceIndex > 0 && sectionIndex > 0; }

	void setOffSurface() { faceIndex = -1; sectionIndex = -1; }

	static const FSurfaceIndex OffSurface;
};