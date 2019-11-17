// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"

#include "MolecularLego/MolecularLego_Relaxation.h"
#include "Simulation/ParticleSlots.h"

#include "DNASimulation.generated.h"	

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_tRNA_pickupAminoAcid : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_tRNA_pickupAminoAcid() {}
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_tRNAEjected : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_tRNAEjected() {}
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRNACodon : public FGraphObject
{
	GENERATED_BODY()

public:
	// Default is the start codon, AUG
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName codon = FName(TEXT("AUG")); 

	static const FName StartCodon;
	static const FName ACA;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FtRNA : public FGraphObject
{
	GENERATED_BODY()

public:
	// Default is the start codon, AUG
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName codon = FName(TEXT("AUG"));

	UPROPERTY() FGraphNodeHandle payloadSlot;

	static const FName StartCodon;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FAminoAcid : public FGraphObject
{
	GENERATED_BODY()

public:

};



UENUM(BlueprintType)
enum class ERibosomeSlotType : uint8
{
	entrance UMETA(DisplayName = "entrance"),
	current UMETA(DisplayName = "current"),
	exit UMETA(DisplayName = "exit"),
}; 

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRibosomeMRNASlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ERibosomeSlotType slotType;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRibosomeAASlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY() uint32 group = 0;
	UPROPERTY() uint32 segmentID = 0;

	UPROPERTY() FGraphNodeHandle lastPigeon = FGraphNodeHandle::null;

public:
	FGraphNodeHandle unbind(FGraph& graph);
	void bind(FGraphNodeHandle target, FGraph& graph);
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRibosomeTRNASlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ERibosomeSlotType slotType;

	UPROPERTY() FGraphNodeHandle ignore0 = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle ignore1 = FGraphNodeHandle::null;
};

// Attached to a tRNA, holds onto an amino acid payload
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FtRNAPayloadSlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	// The FMLElement::type of amino acid to hold onto.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName aaName;
};


USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRibosomeObject : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName type;

	// Pigeon holes
	UPROPERTY() FGraphNodeHandle slot_mRNA_entrance = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_mRNA_cur = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_mRNA_exit = FGraphNodeHandle::null;

	UPROPERTY() FGraphNodeHandle slot_tRNA_entrance = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_tRNA_cur = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_tRNA_exit = FGraphNodeHandle::null;

	UPROPERTY() FGraphNodeHandle slot_aa_exit = FGraphNodeHandle::null;
	
	bool didInit();

	UPROPERTY() float time = 1.0f;
};

// -------------------------------------------------------------
// Polymerase
// -------------------------------------------------------------

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRNAPolymeraseObject : public FGraphObject
{
	GENERATED_BODY()

public:
	// Pigeon holes
	UPROPERTY() FGraphNodeHandle slot_DNA_entrance = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_DNA_cur = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle slot_DNA_exit = FGraphNodeHandle::null;

	UPROPERTY() FGraphNodeHandle slot_mRNA_exit = FGraphNodeHandle::null;

	bool didInit();

	UPROPERTY() float time = 1.0f;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FDNAObject : public FGraphObject
{
	GENERATED_BODY()

public:

};


// Convenience object for finding the tata sites
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FTATAObject : public FGraphObject
{
	GENERATED_BODY()

public:
	FGraphNodeHandle find5Prime(FGraph& graph);
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRNAPolymerase_DNASlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ERibosomeSlotType slotType = ERibosomeSlotType::current;
};

// Attached to a tRNA, holds onto an amino acid payload
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FRNAPolymerase_RNASlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ERibosomeSlotType slotType = ERibosomeSlotType::current;

	// When this fills to 3 items, we create an mRNA codon
	UPROPERTY()
	TArray<FName> unfinishedCodon;

	UPROPERTY() FGraphNodeHandle lastHandle;

	uint32_t segmentGroup = 0;
	uint32_t segmentID = 0;

	void clear();
};

// -------------------------------------------------------------
// Simulation
// -------------------------------------------------------------

UCLASS(BlueprintType)
class LIFEBRUSH_API UDNASimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float ribosomeSeekRadius = 2.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float slotRadius = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float advancementRate = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float slotAcceleration = 15.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float aminoAcidSegmentLength = 1.8f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float rnaSegmentLength = 1.8f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float rnaSegmentRadius = 0.4f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * aminoAcidMaterial = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * rnaMaterial = nullptr;

public:

	virtual void tick(float deltaT) override;

protected:
	virtual void attach() override;
	virtual void detach() override;

protected:
	// -------------------------------------------------------------
	// Ribosomes
	// -------------------------------------------------------------
	void _initRibosomes();
	void _tickRibosomes(float dt);

	void _initTRNAs();
	void _tick_tRNAPayloadSlots(float dt);

	void _tickSlots_ribosomeMRNA(float dt);
	void _tickSlots_ribosomeTRNA(float dt);
	void _tickSlots_ribosomeAA(float dt);

	void _attachAminoAcidToPeptideChain(FGraphNodeHandle aminoAcidHandle, FRibosomeAASlot * aminoAcidSlot);

	FGraphNodeHandle _nextInFilament(FGraphNodeHandle handle);

	// -------------------------------------------------------------
	// DNA Polymerase
	// -------------------------------------------------------------
	void _initRNAPolymerase();
	void _tickPolymerase(float dt);

	void _tickSlots_RNAPolymerase_DNA(float dt);
	void _tickSlots_RNAPolymerase_RNA(float dt);

	FGraphNodeHandle _findNextDNA(FGraphNodeHandle start);
	// Finds a TATA near the nodeHandle
	FGraphNodeHandle _findTATA(FGraphNodeHandle nodeHandle);

	void _attachNucleotideToRNAFilament(FDNAObject& dna, FRNAPolymerase_RNASlot& rnaSlot);
};



