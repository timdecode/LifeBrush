// Copyright (c) 2019 Timothy Davison. All rights reserved.

#pragma once

#include "ShipEditorSimulation/ObjectSimulation.h"
#include "Algorithm/Algorithm.h"
#include "Algorithm/Element.h"
#include "aabbcc/unrealAABB.h"

#include "MolecularLego/InteractionPair.h"

#include <vector>

#include "MolecularLego_Relaxation.generated.h"	

// There are different types of rule positions. Empty, Occupied or Optional. When the conditions
// of a rule position are satisfied and all positions in a rule are satisfied, an interaction
// with the occupying element will be generated.
UENUM(BlueprintType)
enum class ERulePositionType : uint8
{
	// The rule position must be empty (no overlapping elements).
	Empty UMETA(DisplayName = "Empty"), 
	// The rule position must be occupied by one or more elements. The occupied element will get a prediction from the optimization.
	// All occupied positions in a rule must be satisfied to apply the rule.
	Occupied UMETA(DisplayName = "Occupied"),
	// Can be a base position. Must be occupied.
    OccoupiedBase UMETA(DisplayName = "OccoupiedBase"),
	// The rule position may be occupied, if it is, a prediction is generated for the occupying element.
	// Any or none of the rule positions in a rule may be satisfied to apply the rule.
	Optional UMETA(DisplayName = "Optional"),
};

// Lower priority rules should not not overlap with elements that are occluders. The shape of the occluder can be controlled
// by adding FMLOccluderEdgeObjects between the occluder and adjacent nodes in the aggregate.
UENUM(BlueprintType)
enum class EOccluderType : uint8
{
	NoOcclusion UMETA(DisplayName = "NoOcclusion"),

	Occluder UMETA(DisplayName = "Occluder"),

};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLElement : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FName type;	

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	ERulePositionType positionType;

	FGraphNodeHandle rigid;

	bool bound = false;
};

// The node should also have a FMLElement.
USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLElementRule : public FGraphObject
{
	GENERATED_BODY()
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLOccluder : public FGraphObject
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float radius = 1.0f;

	// The main node of this occluder rigid-body, the one that contains the FMLElement.
	FGraphNodeHandle transientHost;
};

USTRUCT(BlueprintType)
struct LIFEBRUSH_API FMLElementRuleConnection : public FGraphEdgeObject
{
	GENERATED_BODY()
};

// ---------------------------------------------------------------------
// Actor Interface
// ---------------------------------------------------------------------
UCLASS()
class LIFEBRUSH_API AMLRuleActor : public AStaticMeshActor
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	FMLElementRule rule;

	// This actors element will be added to the rule if this is true. If all the elements
	// are in child actors, set this to false.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	bool hasElement = true;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	FMLElement element; // a rule node also has an element

	void writeToGraph(FGraph& ruleGraph);
};

UCLASS()
class LIFEBRUSH_API AMLElementActor : public AStaticMeshActor
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "Mitochondria")
	FMLElement element;
};

UCLASS(BlueprintType)
class LIFEBRUSH_API UMLElementSimulation : public UObjectSimulation
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	FGraph ruleGraph;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float minAssignmentDistance = 10.0f;

	// Rules will not be applied to elements when the density at their position is greater than this value.
	// The base density is 1, therefore, this value must be greater than 1.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxAssignmentDensity = 1.8f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float typeCost = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float hackNeighbourhoodRadius = 8.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float maxSpeed = 5.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float sigma = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	float torqueMultiplier = 1.0f; 

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "LifeBrush")
	bool disabled = false;

	float _hackRadius = 1.0f;

public:

	struct RulePosition
	{
		FVector position;
		FQuat orientation;
		FGraphNodeHandle handle;
		int32 reciprocal = -1;
		FName type;
		ERulePositionType positionType = ERulePositionType::Empty;
	};

	struct Rule
	{
		// Positions relative to the subject. The subject must be at zero.
		std::vector<RulePosition> positions;

		int32 priority = 0;

		FGraphNodeHandle ruleHandle;
	};

	struct PredictedPosition
	{
		PredictedPosition(size_t elementIndex, Eigen::Vector3f position, FQuat rotation = FQuat::Identity, float weight = 1.0f) 
			: position(position), rotation(rotation), elementIndex(elementIndex), weight(weight) {}

		Eigen::Vector3f position;
		FQuat rotation;

		int32_t elementIndex;

		float weight = 1.0f;
	};

	struct PositionQuat
	{
		FVector position;
		FQuat rotation;
	};

	struct InteractionPairRuleIndices
	{
		// first corresponds to the min index of an interaction pair \see RuleInteractions
		// and second to the max index of an interaction pair
		int8_t first_base = 0; // assume base position of the rule
		int8_t first_a = -1; // -1 means no assignment
		float first_weight = std::numeric_limits<float>::max();

		int8_t second_base = 0; // assume base position of the rule
		int8_t second_a = -1;
		float second_weight = std::numeric_limits<float>::max();
	};

	struct RuleInteractions
	{
		Rule * rule;

		TMap<FGraphNodeHandle, float> weights;

		// Looks like this:
		// <h_0,h_2> : <r_1,r_3> // neighbours h_0 and h_2 interact, h_0 interacts using the sub-rule index, r_1, and h_2 interacts with h_0 using r_3
		// <h_0,h_5> : <r_2,r_5>
		TMap<InteractionPair, InteractionPairRuleIndices> interactionPairs;

		// Lower priority rules should not not overlap with the handles in this set. The shape of the occluder can be controlled
		// by adding FMLOccluderEdgeObjects between the occluder and adjacent nodes in the aggregate.
		TSet<FGraphNodeHandle> occluders;
	};

	std::vector<Rule> cachedRules;

	typedef unrealAABB::Tree BVH_t;
	BVH_t elementBVH;
	BVH_t occluderBVH;

	int32 _ticksLeftBeforeRebuild;

protected:
	std::vector<RuleInteractions> _previousInteractions;

	std::vector< std::function<void()> > _tickWork;


protected:
	virtual void attach() override;

public:
	virtual void begin() override;
	virtual void detach() override;

	virtual void tick(float deltaT) override;
	virtual void tick_paused(float deltaT) override;

	virtual void componentAdded(FGraphNodeHandle node, ComponentType type) override;
	virtual void componentRemoved(FGraphNodeHandle node, ComponentType type) override;

	auto addTickWork(std::function<void()> work) -> void;

	void _updateBVHs();

protected:
	auto _initCachedRules() -> void;

	auto _initOccluders() -> void;

	auto _initElements() -> void;

	std::vector<RuleInteractions> _interactions();
	void _filterInteractions(RuleInteractions& interactions);
	void _normalizeInteractions(std::vector<RuleInteractions>& interactions);

	void _buildPredictions(
		RuleInteractions& interactions, 
		std::map<FGraphNodeHandle, unsigned int>& problem, 
		TMap<FGraphNodeHandle, FGraphNodeHandle>& rigidHandles,
		std::vector<std::vector<PredictedPosition>>& predictions_in_out);

protected:

	// Finds the greedy minimum assignment of pairs between two aligned neighborhoods
	auto _greedyPairs(
		const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
		const Rule& rule,
		std::vector<std::pair<int16, int16>>& pairings_out) -> float;

	auto _density(const FVector p_a) -> float;

	auto _optimalPairs(
		FGraphNode& aElement,
		const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph,
		const Rule& rule, size_t positionIndex, const RuleInteractions& previousInteractions,
		std::vector<std::pair<int16, int16>>& pairings_out) -> float;

	auto _optimalPairs(FGraphNode& aNode, 
		const std::vector<PositionQuat>& aPositions, const std::vector<FGraphNodeHandle>& aNeighbours, FGraph& aGraph, 
		FGraphNodeHandle ruleHandle, FGraphNodeHandle ruleBase, 
		const RuleInteractions& previousInteractions, 
		std::vector<std::pair<FGraphNodeHandle, FGraphNodeHandle>>& pairings_out) -> float;
	
	auto _ruleSatisfied(Rule& rule, std::vector<std::pair<int16, int16>> pairs) -> bool;
	auto _ruleOccluded(int i_rule, int i_base, FGraphNode& node, std::vector<FGraphNodeHandle>& handles, const std::vector<std::pair<int16, int16>>& pairings) -> bool;

	void _rotationalOptimization(float deltaT, std::vector<RuleInteractions>& interactions);
	void _positionalOptimization(float deltaT, std::vector<RuleInteractions>& interactions);

	FVector4 _vec4(FQuat quat) { return FVector4(quat.X, quat.Y, quat.Z, quat.W); }
	FQuat _quat(FVector4 vec) { return FQuat(vec.X, vec.Y, vec.Z, vec.W); }

	template<class TElementType>
	void _updateBVH(BVH_t& bvh)
	{
		auto& elements = graph->componentStorage<TElementType>();

		// update or insert into the BVH
		for (auto& e : elements)
		{
			auto particleIndex = e.nodeHandle().index;

			if (!e.isValid())
			{
				if (bvh.containsParticle(particleIndex))
					bvh.removeParticle(particleIndex);
			}

			FVector p = graph->node(e.nodeHandle()).position;

			if (!e.isValid()) continue;

			if (bvh.containsParticle(particleIndex))
				bvh.updateParticle(particleIndex, p, _hackRadius);
			else
				bvh.insertParticle(particleIndex, p, _hackRadius);
		}
	}

protected:
	struct Neighbourhood
	{
		FGraphNodeHandle element;

		FGraphNodeHandle ruleElement;

		std::vector< std::pair<FGraphNodeHandle, FGraphNodeHandle> > pairs;

		float weight = 1.0f;
		bool satisfied = false;
	};

	void _createBonds(float deltaT, std::vector< std::vector<Neighbourhood> >& predictionNeighbourhoods);
	void _bind(FGraphNodeHandle rigid_a, FGraphNodeHandle rigid_b, FVector position_b, FQuat rotation_b);
	auto _bindingHandles(FGraphNodeHandle rigidHandle, FVector target, int numHandles)->std::vector<FGraphNodeHandle>;

	void _repulsion(
		std::vector<Neighbourhood> predictionNeighbourhoods,
		std::map<FGraphNodeHandle, unsigned int>& problem,
		std::vector<FGraphNodeHandle>& indexToElement,
		std::vector<std::vector<PredictedPosition>>& predictions_in_out);
};