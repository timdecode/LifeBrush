// Copyright 2016, Timothy Davison. All rights reserved.

#include "LifeBrush.h"

#include "Utility.h"

#include "ElementActor.h"
#include "AggregateProxyComponent.h"
#include "ShipEditorSimulation/MeshSimulation.h"
#include "Simulation/Aggregates.h"
#include "Simulation/FlexElements.h"
#include "ShipEditorSimulation/ObjectSimulation.h"


// Sets default values
AElementActor::AElementActor()
{
 	// Set this actor to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

	SetMobility(EComponentMobility::Movable);
}

// Called when the game starts or when spawned
void AElementActor::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AElementActor::Tick( float DeltaTime )
{
	Super::Tick( DeltaTime );

}

void AElementActor::_loadAggregate(FGraphNodeHandle elementNode, FGraph& graph)
{
	USceneComponent * sceneComponent = this->GetRootComponent();

	// load any aggregate proxy children
	TArray<USceneComponent*> children;
	TArray<FGraphNodeHandle> particlesInBody;
	sceneComponent->GetChildrenComponents(true, children);
	if (children.Num())
	{
		for (USceneComponent * child : children)
		{
			UAggregateProxyComponent * aggregateProxy = Cast<UAggregateProxyComponent>(child);

			if (!aggregateProxy) continue;

			FTransform transform = aggregateProxy->GetComponentTransform();

			FGraphNodeHandle aggregateHandle = FGraphNodeHandle(graph.addNode(
				transform.GetLocation(),
				transform.GetRotation(),
				transform.GetMaximumAxisScale()
			));

			Utility::unboxComponents(aggregateProxy->graphObjects, aggregateHandle(graph), graph);

			if (aggregateHandle.node(graph).hasComponent<FGraphMesh>())
			{
				FGraphMesh& mesh = aggregateHandle.node(graph).component<FGraphMesh>(graph);

				mesh.staticMesh = aggregateProxy->GetStaticMesh();
				mesh.material = aggregateProxy->GetMaterial(0);
			}

			particlesInBody.Add(aggregateHandle);
		}

		if (!elementNode(graph).hasComponent<FFlexParticleObject>())
		{
			UE_LOG(LogTemp, Warning, TEXT("Adding a FFlexParticleObject to the element node, so it's position will move with the rigid body aggregate."));

			elementNode(graph).addComponent<FFlexParticleObject>(graph);
		}
		
		// link them together, as a rigid and an aggregate
		particlesInBody.Add(elementNode);
		auto& rigidBody = FFlexRigidBodyObject::createRigidBody(graph, particlesInBody, FGraphNodeHandle::null /* we want a new node created */ );

		// remove the elementNode
		particlesInBody.Remove(elementNode);
		particlesInBody.Add(rigidBody.nodeHandle());
		FMLAggregateNO::aggregateNodes(particlesInBody, graph, elementNode);
	}
}

void AElementActor::writeToElement(ElementTuple& tuple, FGraph& graph )
{
	tuple.element.radius = overrideRadius;

	USceneComponent * sceneComponent = this->GetRootComponent();

	if(scaleOverrideRadius)
		tuple.element.radius *= sceneComponent->GetComponentTransform().GetMaximumAxisScale();

	tuple.element.generative = generative;

	tuple.element.generationParameters = generationParameters;
	tuple.element.optimizationParameters = optimizationParameters;

	tuple.element.minAssignmentDistance = minAssignmentDistance;
	tuple.element.freespaceRadius = freespaceRadius;

	tuple.element.generationInnerRadius = generationInnerRadius;

	tuple.node.scale = sceneComponent->GetComponentScale().X;

	// give it a mesh
	if(UStaticMeshComponent * staticMeshComponent = FindComponentByClass<UStaticMeshComponent>())
	{
		FGraphMesh& graphMesh = tuple.node.addComponent<FGraphMesh>(graph);

		graphMesh.staticMesh = staticMeshComponent->GetStaticMesh();
		graphMesh.material = staticMeshComponent->GetMaterial(0);
	}


	// unbox the attached components
	Utility::unboxComponents(graphObjects, tuple.node, graph);

	_loadAggregate(tuple.handle(), graph);

	// don't touch the tuple past this point, the reference could be invalidated (retrieve it from the graph again with
	// graph.node(tuple.handle()).
}

void AElementActor::readFromElement(ElementTuple& tuple, FGraph& graph)
{
	FElementObject& element = tuple.element;
	FGraphNode& node = tuple.node;

	overrideRadius = element.radius;
	scaleOverrideRadius = true;
	generative = element.generative;

	generationParameters = element.generationParameters;
	optimizationParameters = element.optimizationParameters;

	minAssignmentDistance = element.minAssignmentDistance;
	freespaceRadius = element.freespaceRadius;

	generationInnerRadius = element.generationInnerRadius;

	USceneComponent * sceneComponent = GetRootComponent();
	sceneComponent->SetWorldScale3D( FVector( node.scale ) );

	graphObjects = Utility::boxComponents(tuple.node, graph);
}

void AElementActor::readFromActor( AElementActor* elementActor )
{
	overrideRadius = elementActor->overrideRadius;
	scaleOverrideRadius = elementActor->scaleOverrideRadius;
	gradient = elementActor->gradient;
	generative = elementActor->generative;

	generationParameters = elementActor->generationParameters;
	optimizationParameters = elementActor->optimizationParameters;

	minAssignmentDistance = elementActor->minAssignmentDistance;
	freespaceRadius = elementActor->freespaceRadius;
	generationInnerRadius = elementActor->generationInnerRadius;

	GetRootComponent()->SetWorldScale3D( elementActor->GetRootComponent()->GetComponentScale() );

	graphObjects = elementActor->graphObjects;
}

