// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "MirrorPoints.generated.h"

UENUM( Blueprintable, Meta = (Bitflags) )
enum class EMirrorAxis : uint8
{
	X,
	Y,
	Z,
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FMirrorAxis
{
	GENERATED_USTRUCT_BODY()

public:


	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor", meta = (Bitmask, BitmaskEnum = "EMirrorAxis") )
	int32 flags;

	void set( EMirrorAxis mirroredAxis ) { flags |= (1 << static_cast<uint32>(mirroredAxis)); }
	void clear( EMirrorAxis mirroredAxis ) { flags &= ~(1 << static_cast<uint32>(mirroredAxis)); }
	bool test( EMirrorAxis mirroredAxis ) { return (flags & (1 << static_cast<uint32>(mirroredAxis))) > 0; }

	FQuat mirror( FQuat& other )
	{
		if(flags == 0)
			return other;

		FVector normal(
			test( EMirrorAxis::X ) ? 1.0f : 0.0f,
			test( EMirrorAxis::Y ) ? 1.0f : 0.0f,
			test( EMirrorAxis::Z ) ? 1.0f : 0.0f
		);

		FPlane symmetryPlane = FPlane( FVector::ZeroVector, normal );

		FVector axis;
		float angle;
		other.ToAxisAndAngle( axis, angle );

		axis = axis.MirrorByPlane( symmetryPlane );
		angle *= -1.0f;

		return FQuat( axis, angle );
	}

	FVector mirror( FVector& other )
	{
		return FVector(
			other.X * (test( EMirrorAxis::X ) ? -1.0f : 1.0f),
			other.Y * (test( EMirrorAxis::Y ) ? -1.0f : 1.0f),
			other.Z * (test( EMirrorAxis::Z ) ? -1.0f : 1.0f)
		);
	}

	operator EAxis::Type()
	{
		EAxis::Type axis = EAxis::None;
		if(test( EMirrorAxis::X ))
			axis = EAxis::X;
		else if(test( EMirrorAxis::Y ))
			axis = EAxis::Y;
		else if(test( EMirrorAxis::Z ))
			axis = EAxis::Z;

		return axis;
	}
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FMirrorPoints
{
	GENERATED_USTRUCT_BODY()

public:
	FMirrorPoints() {}
	~FMirrorPoints() {}

	UPROPERTY()
	FMirrorAxis mirrorredAxis;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
		FPlane symmetryPlane = FPlane( FVector::ZeroVector, FVector::RightVector );

	int32 numSymmetryPoints()
	{
		int n = 1 + int( mirrorredAxis.test( EMirrorAxis::Y ) );

		return n;
	}

	// Returns a mirrored set of locations about all of the current symmetry planes
	// The returned array will not contain the input location unless includeLocation is true.
	TArray<FVector> mirrorLocation( FVector location, bool includeLocation = false )
	{
		TArray<FVector> locations;

		if(includeLocation)
			locations.Add( location );

		if(mirrorredAxis.test(EMirrorAxis::Y) )
		{
			FVector mirrored = location.MirrorByPlane( symmetryPlane );

			locations.Add( mirrored );
		}

		return locations;
	}

	TArray<FQuat> mirrorRotation( FQuat rotation, bool includeLocation = false )
	{
		TArray<FQuat> rotations;

		if(includeLocation)
			rotations.Add( rotation );

		if(mirrorredAxis.test( EMirrorAxis::Y ))
		{
			FVector axis;
			float angle;
			rotation.ToAxisAndAngle( axis, angle );
			
			axis = axis.MirrorByPlane( symmetryPlane );
			angle *= -1.0f;

			rotations.Add( FQuat(axis, angle) );
		}

		return rotations;
	}

	// plural, returns the set of mirrored axis
	TArray<FMirrorAxis> mirrorAxes( bool includeLocation = false )
	{
		TArray<FMirrorAxis> axis;

		if( includeLocation )
			axis.Add( FMirrorAxis() );
		
		if( mirrorredAxis.flags )
			axis.Add( mirrorredAxis );

		return axis;
	}

	TArray<FTransform> mirrorTransform( FTransform& transform, bool includePassed = false )
	{
		auto mirrors = mirrorAxes( includePassed );

		TArray<FTransform> result;

		for(auto& mirror : mirrors)
		{
			EAxis::Type axis = EAxis::None;
			if(mirror.test( EMirrorAxis::X ))
				axis = EAxis::X;
			else if(mirror.test( EMirrorAxis::Y ))
				axis = EAxis::Y;
			else if(mirror.test( EMirrorAxis::Z ))
				axis = EAxis::Z;

			FTransform newTransform = transform;
			newTransform.Mirror( axis, EAxis::None );

			result.Add( newTransform );
		}

		return result;
	}

	TArray<FPlane> activeSymmetryPlanes()
	{
		TArray<FPlane> planes;

		if(mirrorredAxis.test( EMirrorAxis::Y ))
			planes.Add( symmetryPlane );

		return planes;
	}
};
