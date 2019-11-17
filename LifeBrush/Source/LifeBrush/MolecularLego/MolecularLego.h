// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"

#include "Eigen/Dense"

#include "aabbcc/unrealAABB.h"

#include "InteractionPair.h"

#include "MolecularLego.generated.h"


USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLSpeciesID
{
	GENERATED_BODY()

public:
	FMLSpeciesID() : id(-1) {}
	explicit FMLSpeciesID(int32 i) : id(i) {}

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Molecular Lego")
	int32 id = -1;

	inline explicit operator bool() const
	{
		return id >= 0;
	}

	friend bool operator<(const FMLSpeciesID& a, const FMLSpeciesID& b)
	{
		return a.id < b.id;
	}

	bool operator==(const FMLSpeciesID& other)
	{
		return id == other.id;
	}

	bool operator!=(const FMLSpeciesID& other)
	{
		return id != other.id;
	}

	friend inline uint32 GetTypeHash(const FMLSpeciesID& Key)
	{
		return GetTypeHash(Key.id);
	}

	static const FMLSpeciesID null;
};

inline bool operator==(const FMLSpeciesID& A, const FMLSpeciesID& B)
{
	return A.id == B.id;
}

namespace std
{
	template<> struct hash<FMLSpeciesID>
	{
		std::size_t operator()(const FMLSpeciesID& handle) const
		{
			return std::hash<decltype(handle.id)>()(handle.id);
		}
	};
}

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLParticle : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	int32 speciesID = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	int32 numActiveBonds = 0;
};
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLBrownianMotion : public FGraphObject
{
	GENERATED_BODY()

public:


	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	float time = 0.0f;

	FGraphNodeHandle _rigidBody;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLOrientedParticle : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	int32 speciesID = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	int32 numActiveBonds = 0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	float interactionMultiplier = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush", meta = (UIMin = "0"))
	float radiusMultiplier = 1.0f;

	int32 _subsBound = 0;
	std::array<FGraphNodeHandle, 4> subs;

	FGraphNodeHandle _rigidBody;

	FGraphNodeHandle _bondedPartner;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool active = true;

	bool isSlave = false;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool alignedToBind = false;

	// don't directly manipulate this stuff
	bool _didInit;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLOrientedConnection : public FGraphEdgeObject
{
	GENERATED_BODY()

public:
	float targetDistance = 2.0f;
	float breakingDistance = 100.0f;
};



UCLASS(BlueprintType)
class LIFEBRUSH_API UMLParticleSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * debugBindingSiteMesh = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * debugBindingSiteMaterial = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float debugBindingSiteScale = 0.01f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float baseBindingStrength = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float bindingOffset = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float minTime = 0.05f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxTime = 0.3f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float minSpeed = 1.0f;
	
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxSpeed = 4.0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float minTau = 40.0;
	
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxTau = 100.0;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float flexBindingRadius = 0.3f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float hackNeighbourhoodRadius = 8.0f;


public:
	// Recently bonded FMLOrientedParticles
	TArray<FGraphNodeHandle> recentlyBonded;

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void tick(float deltaT) override;

	static FQuat mirror(FQuat rotation);

protected:
	void tick_mlParticles(float deltaT);
	void tick_bindOrientedParticles(float deltaT);
	void tick_bindOrientedParticles2(float deltaT);
	void tick_bindOrientedParticles_greedyFlexConnections(float deltaT);
	void tick_bindOrientedParticles4(float deltaT);
	void tick_bindOrientedParticles_sortedInteractions(float deltaT);

	// find and cache our helper particles
	void _cacheBindingSubs();

	FVector _torqueSpring(const FQuat &dq_a, float& torqueMagnitude_out);

	void centerOfMass(FGraphNodeHandle rigidHandle, int32 &rigidParticleCount, FVector &com);

	void _updateParticleBVH();

	TArray<FGraphNodeHandle> _createStar(
		FVector position, 
		FQuat orientation,
		float scale, 
		float separation,
		FGraphNodeHandle rigidBody,
		UStaticMesh * optionalMesh= nullptr,
		UMaterialInterface * optionalMaterial = nullptr);

	bool _facing(FGraphNode& a, FGraphNode& b);

	// critically damped positional spring
	// see: http://mathproofs.blogspot.com/2013/07/critically-damped-spring-smoothing.html

	bool _aligned(FGraphNode& node_a, FMLOrientedParticle& particle_a, FGraphNode& node_b, FMLOrientedParticle& particle_b);

	void tick_brownianMotion(float deltaT);

protected:


protected:
	// Attraction (and repulsion) matrix. c_ij is the attraction between particle
	// species i and j.
	// We need to figure out a way to expose this to Unreal in a nice way.
	Eigen::MatrixXf _attractionMatrix;

	Eigen::MatrixXf _bindingRadiusMatrix;



	float _hackRadius = 1.0f;

	typedef unrealAABB::Tree BVH_t;
	
	// indices are raw-indices into the FMLParticle-storage.
	BVH_t _particleBVH;
	BVH_t _orientedParticleBVH;
	BVH_t _bindingParticleBVH;

	// The number of bonds between two rigid body handles
	TMap<InteractionPair, int32> _numBonds;

	// Each chain gets its own id.
	// This maps rigid-body handles to the chain id
	TMap<FGraphNodeHandle, int32> _chainID;
	int32 _nextChain = 0;
};



UCLASS(BlueprintType)
class LIFEBRUSH_API UMolecularLegoSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:

};