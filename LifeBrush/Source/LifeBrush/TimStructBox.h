// Copyright 2018, Timothy Davison. All rights reserved.

// Based on Epic's FTimStructBox in StructBox.h

#pragma once

#include "TimStructBox.generated.h"

USTRUCT( BlueprintType )
struct FTimStructBox
{
	GENERATED_USTRUCT_BODY()

public:
	UPROPERTY()
	UScriptStruct* scriptStruct;

	// we should figure out how to make this a UPROPERTY
	uint8* structMemory;

	FTimStructBox()
		: structMemory( nullptr )
	{}

	bool IsValid() const
	{
		return scriptStruct && structMemory;
	}

	void initStruct();

	~FTimStructBox();

	bool Serialize( FArchive& Ar );

	bool ExportTextItem(FString& ValueStr, FTimStructBox const& DefaultValue, UObject* Parent, int32 PortFlags, UObject* ExportRootScope) const;
	bool ImportTextItem(const TCHAR*& Buffer, int32 PortFlags, UObject* Parent, FOutputDevice* ErrorText);

	bool Identical( const FTimStructBox* Other, uint32 PortFlags ) const;

	void AddStructReferencedObjects( class FReferenceCollector& Collector ) const;

	FTimStructBox& operator=( const FTimStructBox& Other );

protected:
	void Destroy( UScriptStruct* ActualStruct );
	void Create( UScriptStruct* ActualStruct ); 
	
public:
	FTimStructBox( const FTimStructBox& );
};

template<> struct TIsPODType<FTimStructBox> { enum { Value = false }; };
template<> struct TIsTriviallyCopyAssignable<FTimStructBox> { enum { Value = false }; };

template<>
struct TStructOpsTypeTraits<FTimStructBox> : public TStructOpsTypeTraitsBase2<FTimStructBox>
{
	enum
	{
		WithZeroConstructor = true,
		WithCopy = true,
		WithIdentical = true,
		WithAddStructReferencedObjects = true,
		WithSerializer = true,
		WithExportTextItem = true,
		WithImportTextItem = true,
		// TODO.. WithPostSerialize etc..
	};
};

