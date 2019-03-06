// Copyright 2017 Code Monkey Castle, all rights reserved.

#pragma once

#include "Components/ActorComponent.h"
#include "MarchingCubes.h"
#include "RuntimeMeshComponent.h"
#include <atomic>

#include "VolumeComponent.generated.h"

UENUM( Blueprintable )
enum class EBlocksGenerationMethod : uint8
{
	FloatingIsland,
	OriginalFloatingIsland,
	Sphere,
	Empty,
};

UCLASS( ClassGroup = (Custom), meta = (BlueprintSpawnableComponent) )
class LIFEBRUSH_API UVolumeComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	UMaterialInterface * material;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	FVector2D uvScale = FVector2D( 1.0f, 1.0f );

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 size = 64;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	bool flatShading = true;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	EBlocksGenerationMethod method = EBlocksGenerationMethod::FloatingIsland;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float isoLevel = 3.1f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 randomSeed = 426472;

	UPROPERTY( EditAnywhere, Category = "ShipEditor" )
	uint32 numWorms = 10;

	UPROPERTY( EditAnywhere, Category = "ShipEditor" )
	float wormCellRadius = 4.0f;	// the worm's radius is wormCellRadius * the cell radius

	UPROPERTY( EditAnywhere, Category = "ShipEditor" )
	FVector wormNoiseScale = FVector( 1.0f, 1.0f, 1.0f );

public:
	// Sets default values for this component's properties
	UVolumeComponent();

	// Called when the game starts
	virtual void BeginPlay() override;

	void build();

	UniformGrid<float>& grid();

	FIntVector componentToIndex( FVector localPoint );
	FVector indexToComponent( FIntVector index );

	void markDirty( FVector minExtents, FVector maxExtents );

protected:
	UniformGrid<float> _grid;

	URuntimeMeshComponent* runtimeMesh();

	std::vector<int32> _cachedGridIndices; // indexed by triangle indices (first vertex index / 3)

	std::vector<FVector> perlin_worm( FIntVector start, const unsigned int nSegments );
	void decrepify( UniformGrid<float>& grid, size_t wormSegments, uint32 numWorms, float wormRadius, FRandomStream& randomStream );

	std::atomic_bool _marching = false;
};

UCLASS(ClassGroup = (Custom), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UChunkedVolumeComponent : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	UMaterialInterface * material;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	FVector2D uvScale = FVector2D(1.0f, 1.0f);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	int32 size = 64;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	bool flatShading = true;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	EBlocksGenerationMethod method = EBlocksGenerationMethod::FloatingIsland;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float isoLevel = 3.1f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	int32 randomSeed = 426472;

	UPROPERTY(EditAnywhere, Category = "ShipEditor")
	uint32 numWorms = 10;

	UPROPERTY(EditAnywhere, Category = "ShipEditor")
	float wormCellRadius = 4.0f;	// the worm's radius is wormCellRadius * the cell radius

	UPROPERTY(EditAnywhere, Category = "ShipEditor")
	FVector wormNoiseScale = FVector(1.0f, 1.0f, 1.0f);

public:
	// Sets default values for this component's properties
	UChunkedVolumeComponent();

	// Called when the game starts
	virtual void BeginPlay() override;

	void init();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void initializeTestVolume();

	UFUNCTION(BlueprintCallable, Category = LifeBrush)
	void markTestVolumeDirty();

	ChunkGrid<float>& grid();

	FIntVector componentToIndex(FVector localPoint);
	FVector indexToComponent(FIntVector index);

	void markDirty(FVector minExtents, FVector maxExtents);
	void markDirtyIndex(FIntVector minGridIndex, FIntVector maxGridIndex);

protected:
	ChunkGrid<float> _grid;

	FVector _cellSize;
	FVector _invCellSize;

	URuntimeMeshComponent* runtimeMesh();

	struct ChunkInfo
	{
		std::vector<FIntVector> chunkCellIndexForVertex;

		int sectionIndex = -1;
	};

	std::unordered_map<FIntVector, ChunkInfo> _chunkInfo;

	std::unordered_map<FIntVector, std::vector<int32>> _cachedGridIndices; // indexed by triangle indices (first vertex index / 3)


	std::vector<FVector> perlin_worm(FIntVector start, const unsigned int nSegments);
	void decrepify(UniformGrid<float>& grid, size_t wormSegments, uint32 numWorms, float wormRadius, FRandomStream& randomStream);

	std::atomic_bool _marching = false;
};








