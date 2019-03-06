// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include "ObjectRecord.generated.h"



// https://answers.unrealengine.com/questions/35618/savingloading-an-array-of-objects.html

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FObjectRecord
{
	GENERATED_USTRUCT_BODY()

public:
	FObjectRecord() {}
	FObjectRecord( UObject * object )
	{
		objectName = object->GetName();
		objectClass = object->GetClass();
	}
	~FObjectRecord() {}

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	FString objectName;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UClass * objectClass;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	TArray<uint8> data;

public:
	template<class ElementType>
	static TArray<FObjectRecord> toObjectRecords( TArray<ElementType*> objects )
	{
		TArray<FObjectRecord> records;

		for(auto object : objects)
		{
			auto index = records.Emplace( object ); 
			auto& record = records[index];

			UObject * asObject = Cast<UObject>( object );

			record.objectName = asObject->GetName();
			record.objectClass = asObject->GetClass();

			FMemoryWriter writer( record.data, true );

			FObjectAndNameAsStringProxyArchive recordArchive( writer, true );

			object->Serialize( recordArchive );
		}

		return records;
	}

	template<class ElementType>
	static TArray<ElementType*> fromObjectRecords( TArray<FObjectRecord> records )
	{
		TArray<ElementType*> objects;
		for(auto& record : records)
		{
			UObjectSimulation * object = NewObject<UObjectSimulation>( (UObjectSimulation*)GetTransientPackage(), record.objectClass, FName( *record.objectName ) );
			objects.Add( object ); 

			FMemoryReader reader( record.data );

			FObjectAndNameAsStringProxyArchive recordArchive( reader, true );

			object->Serialize( recordArchive );
		}

		return objects;
	}
};

FORCEINLINE FArchive& operator<<( FArchive& archive, FObjectRecord * record )
{
	archive << record->objectName;
	archive << record->objectClass;
	archive << record->data;

	return archive;
}

FORCEINLINE FArchive& operator<<( FArchive& archive, FObjectRecord& record )
{
	archive << record.objectName;
	archive << record.objectClass;
	archive << record.data;

	return archive;
}