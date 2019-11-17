// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/Graph.h"

#include "aabbcc/unrealAABB.h"
#include "Simulation/FlexGraphSimulation_interface.h"
#include <set>
#include "Simulation/ParticleSlots.h"

#include "Synapse.generated.h"


// -------------------------------------------------------------------------------
// Events
// -------------------------------------------------------------------------------

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_VolgateGatedSodiumChannel_passSodium : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_VolgateGatedSodiumChannel_passSodium() {}
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_AChR_bindACh : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_AChR_bindACh() {}
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_LeuT_passSodium : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_LeuT_passSodium() {}
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_AChR_passSodium : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_AChR_passSodium() {}
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UEvent_AChE_breakACh : public USEGraphEvent
{
	GENERATED_BODY()

public:
	virtual ~UEvent_AChE_breakACh() {}
};


// -------------------------------------------------------------------------------
// Synapse Node objects
// -------------------------------------------------------------------------------

UENUM(BlueprintType)
enum class EMolecularChannelType : uint8
{
	BiDirectional UMETA(DisplayName = "BiDirectional"),
	BottomToTop UMETA(DisplayName = "BottomToTop"),
	TopToBottom UMETA(DisplayName = "TopToBottom"),
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMolecularChannel : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY() bool open = true;

	UPROPERTY() float timer = 0.0f;

	// this isn't really polymorphic, but the template code behaves like it is
	static EMolecularChannelType channelType() { return EMolecularChannelType::BiDirectional; }
};


USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSodiumPotassiumPump : public FMolecularChannel
{
	GENERATED_BODY()

public:
	static EMolecularChannelType channelType() { return EMolecularChannelType::BottomToTop; }
};

// Uses 2 Na+ to carry choline across the membrane.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FLeuTSodiumCoupledTransporter : public FMolecularChannel
{
	GENERATED_BODY()

public:
	static EMolecularChannelType channelType() { return EMolecularChannelType::BottomToTop; }

public:
	UPROPERTY() FGraphNodeHandle cholineHandle = FGraphNodeHandle::null;

	UPROPERTY() FGraphNodeHandle naclHandle = FGraphNodeHandle::null;
};


// https://en.wikipedia.org/wiki/Ball_and_chain_inactivation
UENUM(BlueprintType)
enum class EVoltageGatedIonChannelState : uint8
{
	Closed UMETA(DisplayName = "Closed"),
	Open UMETA(DisplayName = "Open"),
	Inactivated UMETA(DisplayName = "Inactivated"),
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMovingAverage
{
	GENERATED_BODY()

public:

	UPROPERTY() TArray<float> values;
	UPROPERTY() int32 maxValues = 100;
	UPROPERTY() int32 nextIndex = 0;

	float average();
	void addValue(float value);
	void clear();
};


USTRUCT(BlueprintType)
struct LIFEBRUSH_API FVoltageGatedIonChannel : public FMolecularChannel
{
	GENERATED_BODY()

public:
	static EMolecularChannelType channelType() { return EMolecularChannelType::BiDirectional; }

public:
	UPROPERTY() FVector voltage;

	// for moving average
	UPROPERTY() FMovingAverage averageValues;

	UPROPERTY() float previousAverage = 0.0f;

	// for logging
	UPROPERTY() float deltaVoltage = 0.0f;

	UPROPERTY() int upCount = 0;
	UPROPERTY() int downCount = 0;

	UPROPERTY() EVoltageGatedIonChannelState state = EVoltageGatedIonChannelState::Closed;
};

// ----------------------------------------------------------
// Acetylcholine stuff
// ----------------------------------------------------------

// https://pdb101.rcsb.org/motm/71 David Goodsell 2005
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FAcetylcholineReceptor : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "LifeBrush" )
	int counter = 0;

	UPROPERTY() FGraphNodeHandle achSlot_a = FGraphNodeHandle::null;
	UPROPERTY() FGraphNodeHandle achSlot_b = FGraphNodeHandle::null;

	UPROPERTY() bool open = false;

	UPROPERTY() float timer = 0.0f;

	bool didInit();

	static EMolecularChannelType channelType() { return EMolecularChannelType::BiDirectional; }
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FCooldownObject : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY() float cooldown = 0.0f;

	bool onCooldown() { return cooldown > 0.0f; }
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FCholine : public FGraphObject
{
	GENERATED_BODY()

public:

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FAcetylcholine : public FCooldownObject
{
	GENERATED_BODY()

public:

};

// https://pdb101.rcsb.org/motm/54 David Goodsell 2004
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FAcetylcholinesterase  : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY() float cooldown = 0.0f;

	bool onCooldown() { return cooldown > 0.0f; }
};

// Attached to a tRNA, holds onto an amino acid payload
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FAChSlot : public FParticleSlot
{
	GENERATED_BODY()

public:
	UPROPERTY() float timeLeft = -1.0f;

	bool canBind(FGraphNodeHandle target, FGraph& graph);
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FSynapseTimerObject : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY()
	float timeLeft = 0.0f;


};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FNaCl : public FCooldownObject
{
	GENERATED_BODY()

public:

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FPotassium : public FCooldownObject
{
	GENERATED_BODY()

public:

};

UCLASS( Blueprintable )
class LIFEBRUSH_API USynapseSimulation : public UObjectSimulation, public IFlexGraphSimulation
{
	GENERATED_BODY()

public:
	virtual void tick(float deltaT) override;

protected:
	virtual void attach() override;
	virtual void detach() override;

public:
	virtual void componentAdded(FGraphNodeHandle node, ComponentType type) override;
	virtual void componentRemoved(FGraphNodeHandle node, ComponentType type) override;


	virtual void flexTick(float deltaT, 
		NvFlexVector<int>& neighbourIndices, 
		NvFlexVector<int>& neighbourCounts,
		NvFlexVector<int>& apiToInternal, 
		NvFlexVector<int>& internalToAPI, 
		int maxParticles) override;


public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FVector bindingPosition_a_ACh = FVector(-1.8f, 3.0f, 0.0f);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FVector bindingPosition_b_ACh = FVector(1.8f, 3.0f, 0.0f);

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FVector forwardDirection = FVector(0.0f, 1.0f, 0.0f);


	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float bindingRadius_ACh_Receptor = 3.0f;

	// The expected value for how long ACh can remain in the bound state
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float expectedBindingTime_ACh = 2.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float generalRefractoryTime = 4.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float refractoryTime_AChSlot = 4.0f;

	// This is a positive number.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelActivation_deltaVoltage = -5.0f;

	// Once the channel is open, how long does it stay open?
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannenActivationPeriod = 5.0f;

	// When the channel becomes inactivated, how long does it stay inactivated before it can open again?
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannenInactivationPeriod = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannenInactivationSampleRadius = 10.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumPumpRefractoryPeriod = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float acetycholinesteraseRefractoryPeriod = 0.5f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInstance * cholineMaterial = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UStaticMesh * cholineMesh = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float cholineScale = 1.0f;





	// When Closed and voltage > sodiumChannelActivationVoltage => Open
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelActivationVoltage = -1.0f;

	// When Open and voltage > sodiumChannelInactivationVoltage => Inactivate
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelInactivationVoltage = -0.5f;

	// When Inactivated and voltage < sodiumChannelClosedVoltage => Closed
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelClosedVoltage = -1.0f;


	// When Closed and voltage > sodiumChannelActivationVoltage => Open
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelActivationRatio = 0.95;

	// When Open and voltage > sodiumChannelInactivationVoltage => Inactivate
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelInactivationRatio = 0.6;

	// When Inactivated and voltage < sodiumChannelClosedVoltage => Closed
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sodiumChannelClosedRatio = 0.75f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * sodiumChannelMaterial_openState = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * sodiumChannelMaterial_closedState = nullptr;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	UMaterialInterface * sodiumChannelMaterial_inactivatedState = nullptr;

protected:
	void _updateReceptorBVH();
	void _updateAChSlotBVH();
	void _updateVoltageGatedChannelsBVH();
	void _updateSodiumPotassiumPumpBVH();
	void _updateAChEBVH();
	void _updateLeuTBVH();

	void _tickTimers(float deltaT);

	void _bindToAChSlots(float dt);
	void _tickAChs(float dt);

	void _expireAChSlots(float dt);

	void _initReceptors(float dt);
	void _openCloseReceptors(float dt);

	void _openCloseVoltageGatedChannels(float dt);

	void _tickNaCl_and_potassium(float dt);
	void _passNaCl(float dt);
	void _tickLeuT(float dt);

	void _pumpSodium(float dt);
	void _pumpPotassium(float dt);

	void _tickAcetycholinesterase(float dt);

	template<class TSlotType, typename TWork>
	void tickSlots(float dt, TWork operation)
	{
		auto& slots = graph->componentStorage<TSlotType>();
		auto& rigids = graph->componentStorage<FFlexRigidBodyObject>();

		for (TSlotType& slot : slots)
		{
			if( !slot.isValid() ) continue;

			// hold onto the pigeon
			if (slot.occupied())
			{
				FGraphNodeHandle pigeon = slot.pigeon();
				FVector slotPosition = graph->node(slot.nodeHandle()).position;
				FQuat slotOrientation = graph->node(slot.nodeHandle()).orientation;

				// teleport
				FGraphNodeHandle rigidHandle = FFlexRigidBodyObject::getRigidBodyHandle(*graph, pigeon);

				if (rigidHandle)
				{
					rigids.componentPtrForNode(rigidHandle)->setRotationPosition(slotOrientation, slotPosition, *graph);
				}
				else
				{
					graph->node(pigeon).position = slotPosition;
					graph->node(pigeon).orientation = slotOrientation;
				}
			}

			operation(slot);
		}
	}

	typedef unrealAABB::Tree BVH_t;

	template<class TFCooldownObject, class TChannel>
	void _tickChannels(float dt, BVH_t& bvh, float refactoryTime, float pumpRefactoryTime)
	{
		auto nothing = [](TFCooldownObject& element, TChannel& channel) {};

		_tickChannels_withCallback<TFCooldownObject, TChannel>(dt, bvh, refactoryTime, pumpRefactoryTime, nothing);
	}

	template<class TFCooldownObject, class TChannel, class TOperation>
	void _tickChannels_withCallback(float dt, BVH_t& bvh, float refactoryTime, float pumpRefactoryTime, TOperation pumpEvent)
	{
		auto& elements = graph->componentStorage<TFCooldownObject>();
		auto& channels = graph->componentStorage<TChannel>();
		auto& velocities = graph->componentStorage<FVelocityGraphObject>();

		for (TChannel& channel : channels)
		{
			channel.timer -= dt;
		}

		for (TFCooldownObject& element : elements)
		{
			if( !element.isValid() ) continue;
			if( element.onCooldown() ) continue;

			FGraphNode& elementNode = graph->node(element.nodeHandle());

			// find a receptor 
			unrealAABB::AABB query(elementNode.position, bindingRadius_ACh_Receptor);

			FGraphNodeHandle channelHandle = FGraphNodeHandle::null;

			FVelocityGraphObject * velocity = nullptr;
			if (TChannel::channelType() == EMolecularChannelType::BiDirectional)
			{
				velocity = velocities.componentPtrForNode(element.nodeHandle());
			}




			bvh.query(query, [&](unsigned int index_channel) {
				FGraphNodeHandle handle = FGraphNodeHandle(index_channel);

				if (velocity && TChannel::channelType() == EMolecularChannelType::BiDirectional)
				{
					FGraphNode& channelNode = graph->node(handle);

					// check the direction of the velocity, we'll transport if we are moving towards the channel
					FVector dir = channelNode.position - elementNode.position;

					if (FVector::DotProduct(velocity->linearVelocity, dir) >= 0.0f)
					{
						// we found it, we're done the query
						channelHandle = handle;
						return false;
					}

					// keep searching
					return true;
				}

				channelHandle = handle;

				// done the query, exit it
				return false;
			});

			if (!channelHandle) continue;

			TChannel * pump = channels.componentPtrForNode(channelHandle);

			if( !pump ) continue;
			if( !pump->open ) continue;
			if( pump->timer > 0.0f ) continue;

			element.cooldown = refactoryTime;

			// and teleport
			FGraphNode& channelNode = graph->node(channelHandle);

			FVector dir = (elementNode.position - channelNode.position).GetSafeNormal();

			FVector pumpDir = channelNode.orientation.RotateVector(FVector::UpVector);

			bool shouldPump = false;

			if (TChannel::channelType() == EMolecularChannelType::BiDirectional)
				shouldPump = true;
			else if (TChannel::channelType() == EMolecularChannelType::BottomToTop)
				shouldPump = FVector::DotProduct(dir, pumpDir) < 0.0f;
			else if (TChannel::channelType() == EMolecularChannelType::TopToBottom)
				shouldPump = FVector::DotProduct(dir, pumpDir) > 0.0f;

			if (shouldPump)
			{
				elementNode.position = channelNode.position - dir * bindingRadius_ACh_Receptor;

				if (FVelocityGraphObject * velocioty = velocities.componentPtrForNode(elementNode.handle()))
					velocioty->linearVelocity = -dir * 8.0f;

				pump->timer = pumpRefactoryTime;

				pumpEvent(element, *pump);
			}
		}
	}


	template<class TElementType>
	void _updateBVH(BVH_t& bvh);

protected:
	BVH_t _receptorBVH;
	BVH_t _receptorSlotsBVH;
	BVH_t _voltageGatedChannelsBVH;
	BVH_t _sodiumPotassiumPumpBVH;
	BVH_t _AChEBVH;
	BVH_t _LeuTBVH;

	BVH_t _achBVH;

	TSet<FGraphNodeHandle> _receptorsToRemove;

	float _hackRadius = 2.0f;
};

template<class TElementType>
void USynapseSimulation::_updateBVH(BVH_t& bvh)
{
	auto& elements = graph->componentStorage<TElementType>();

	// update or insert into the BVH
	for (auto& e : elements)
	{
		auto particleIndex = e.nodeHandle().index;

		if (!e.isValid() && bvh.containsParticle(particleIndex))
		{
			bvh.removeParticle(particleIndex);
		}
		else
		{
			FVector p = graph->node(e.nodeHandle()).position;

			if (bvh.containsParticle(particleIndex))
				bvh.updateParticle(particleIndex, p, _hackRadius);
			else
				bvh.insertParticle(particleIndex, p, _hackRadius);
		}
	}
}
