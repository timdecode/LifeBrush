// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include <Eigen/Dense>

#include "GraphUtility.h"

#include <vector>
#include <algorithm>
#include <numeric>

// From: http://www.highprogrammer.com/alan/windev/visualstudio.html
// Statements like:
//		#pragma message(Reminder "Fix this problem!")
// Which will cause messages like:
//		C:\Source\Project\main.cpp(47): Reminder: Fix this problem!
// to show up during compiles.  Note that you can NOT use the
// words "error" or "warning" in your reminders, since it will
// make the IDE think it should abort execution.  You can double
// click on these messages and jump to the line in question.
#define Stringize( L )			#L
#define MakeString( M, L )		M(L)
#define $Line					\
	MakeString( Stringize, __LINE__ )
#define Reminder				\
	__FILE__ "(" $Line ") : Reminder: "

/**
 * 
 */
class LIFEBRUSH_API ShipEditorUtility
{
public:
	template <typename T>
	static FORCEINLINE T* LoadObjectFromPath( const FName& Path )
	{
		if(Path == NAME_None)
			return nullptr;

		return Cast<T>( StaticLoadObject( T::StaticClass(), NULL, *Path.ToString() ) );
	}

	static FORCEINLINE UScriptStruct* structNamed( FString name )
	{
		UObject * package = ANY_PACKAGE;

		if(UScriptStruct * scriptStruct = FindObject<UScriptStruct>( package, *name ))
			return scriptStruct;

		if(UObjectRedirector * redirector = FindObject<UObjectRedirector>( package, *name ))
			return (UScriptStruct*)CastChecked<UScriptStruct>( redirector->DestinationObject );

		return nullptr;
	}

	static FORCEINLINE void derivedStructs( UScriptStruct * base, TArray<UScriptStruct*> & derived )
	{
		for(TObjectIterator<UScriptStruct> it; it; ++it)
		{
			if(it->IsChildOf( base ))
				derived.Add( *it );
		}
	}

	static FORCEINLINE FVector scale( const FIntVector& intVector, const FVector& floatVector )
	{
		FVector result;
		for(int c = 0; c < 3; ++c)
			result[c] = intVector[c] * floatVector[c];

		return result;
	}

	// For whatever stupid reason, Epic doesn't export ULevelStreaming::CreateInstance (this class is using the MinimalAPI specifier).
	// So, this is a copy/pasted static version of that.
	static ULevelStreaming * createInstance( ULevelStreaming * levelStreaming, FString uniqueName );
};








// -----------------------------------------------------------------------------
// - Face Windings
// -----------------------------------------------------------------------------

static std::vector<size_t> sortPointsClockwise( std::vector<FVector> points, FPlane plane, FVector centre, bool counterClockWise = false )
{
	std::vector<FVector> projected;

	FVector normal = FVector( plane.X, plane.Y, plane.Z );

	for(FVector& point : points)
	{
		FVector p = FVector::PointPlaneProject( point, plane );
		projected.push_back( p );
	}

	centre = FVector::PointPlaneProject( centre, plane );

	// sort the faces
	// we don't sort _face, we sort an auxiliary array
	std::vector<size_t> indices( points.size() );
	std::iota( indices.begin(), indices.end(), 0 );

	std::sort( indices.begin(), indices.end(),
		[&]( size_t ai, size_t bi ) -> bool
	{
		FVector& a = projected[ai];
		FVector& b = projected[bi];

		// from http://stackoverflow.com/a/14371081
		// > 0 is counter-clockwise, < 0 is clockwise
		float dir = FVector::DotProduct(normal, FVector::CrossProduct( a - centre, b - centre ) );

		return counterClockWise ? dir > 0.0f : dir <= 0.0f;
	} );

	return indices;
}

inline static float triangleAreaSqrd( FVector& v0, FVector& v1, FVector& v2 )
{
	FVector e0 = v1 - v0;
	FVector e1 = v2 - v0;
	FVector n = FVector::CrossProduct( e0, e1 );

	return n.SizeSquared() * 0.25f;
}

