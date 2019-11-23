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
	TArray<FGraphNodeHandle> membersInBody;
	TArray<FGraphNodeHandle> aggregatesInBody;

	TSet<AActor*> visitedActors;
	visitedActors.Add(this);

	sceneComponent->GetChildrenComponents(true, children);

	for (USceneComponent * child : children)
	{
		AElementActor * elementActor = Cast<AElementActor>(child->GetOwner());

		if (elementActor && !visitedActors.Contains(elementActor))
		{
			visitedActors.Add(elementActor);

			_loadSubElementActor(elementActor, graph, particlesInBody, membersInBody, aggregatesInBody);
		}

		if (UAggregateProxyComponent * aggregateProxy = Cast<UAggregateProxyComponent>(child)) 
		{
			_loadAggregateProxy(aggregateProxy, graph, particlesInBody, membersInBody, aggregatesInBody);
		}
	}

	// only create an aggregate (adding ourselves) when we loaded some subunits above
	if (membersInBody.Num() > 0)
	{
		// insert at the front, because we want this actor to be the root
		if (this->particleOrMember == EAggregateProxyParticleOrMember::Particle)
			particlesInBody.Insert(elementNode, 0);
		else if (this->particleOrMember == EAggregateProxyParticleOrMember::Member)
			membersInBody.Insert(elementNode, 0);

		bool isFlexRigid = particlesInBody.Num() > 0;

		if (isFlexRigid && !elementNode(graph).hasComponent<FFlexParticleObject>())
		{
			UE_LOG(LogTemp, Warning, TEXT("Adding a FFlexParticleObject to the element node, so it's position will move with the rigid body aggregate."));

			elementNode(graph).addComponent<FFlexParticleObject>(graph);
		}

		// link them together, as a rigid and an aggregate
		if (isFlexRigid)
		{
			auto& rigidBody = FFlexRigidBodyObject::createRigidBody(graph, particlesInBody, FGraphNodeHandle::null /* we want a new node created */);
			rigidBody.node(graph).orientation = elementNode(graph).orientation;
			aggregatesInBody.Add(rigidBody.nodeHandle());

			FMLAggregateNO::aggregateNodes(aggregatesInBody, graph, elementNode);

			// finally, connect the members
			for (FGraphNodeHandle member : membersInBody)
			{
				graph.connectNodes<FFlexRigidBodyConnection>(rigidBody.nodeHandle(), member);
				member(graph).addComponent<FFlexRigidMember>(graph);
			}
		}
		else
			FMLAggregateNO::aggregateNodes(aggregatesInBody, graph, elementNode);
	}
}


void AElementActor::_loadSubElementActor(
	AElementActor * elementActor, 
	FGraph &graph, 
	TArray<FGraphNodeHandle> &particlesInBody,
	TArray<FGraphNodeHandle> &membersInBody,
	TArray<FGraphNodeHandle> &aggregatesInBody)
{
	FTransform transform = elementActor->GetRootComponent()->GetComponentTransform();

	FGraphNodeHandle elementActorHandle = FGraphNodeHandle(graph.addNode(
		transform.GetLocation(),
		transform.GetRotation(),
		transform.GetMaximumAxisScale()
	));

	FGraphNode& elementActorNode = graph.node(elementActorHandle);

	// unbox the attached components
	Utility::unboxComponents(elementActor->graphObjects, elementActorNode, graph);

	// give it a mesh
	UStaticMeshComponent * staticMeshComponent = elementActor->FindComponentByClass<UStaticMeshComponent>();
	if (staticMeshComponent && elementActor->shouldGenerateMeshObject)
	{
		FGraphMesh& graphMesh = elementActorNode.addComponent<FGraphMesh>(graph);

		graphMesh.staticMesh = staticMeshComponent->GetStaticMesh();
		graphMesh.material = staticMeshComponent->GetMaterial(0);
	}

	if (elementActor->particleOrMember == EAggregateProxyParticleOrMember::Particle)
		particlesInBody.Add(elementActorHandle);
	else if (elementActor->particleOrMember == EAggregateProxyParticleOrMember::Member)
		membersInBody.Add(elementActorHandle);

	aggregatesInBody.Add(elementActorHandle);
}

void AElementActor::_loadAggregateProxy(UAggregateProxyComponent * aggregateProxy, FGraph &graph, TArray<FGraphNodeHandle> &particlesInBody, TArray<FGraphNodeHandle> &membersInBody, TArray<FGraphNodeHandle> &aggregatesInBody)
{


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

	if (aggregateProxy->particleOrMember == EAggregateProxyParticleOrMember::Particle)
		particlesInBody.Add(aggregateHandle);
	else if (aggregateProxy->particleOrMember == EAggregateProxyParticleOrMember::Member)
		membersInBody.Add(aggregateHandle);

	aggregatesInBody.Add(aggregateHandle);
}

void AElementActor::_loadSubElements(FGraphNodeHandle elementNode, TArray<FGraphNodeHandle>& particlesInBody, TArray<FGraphNodeHandle>& membersInBody, TArray<FGraphNodeHandle>& aggregatesInBody)
{

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
	UStaticMeshComponent * staticMeshComponent = FindComponentByClass<UStaticMeshComponent>();
	if( staticMeshComponent && shouldGenerateMeshObject )
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

void AElementActor::_setSelectionOutlineVisibility(bool visibility)
{
	UStaticMeshComponent * mesh = FindComponentByClass<UStaticMeshComponent>();

	// we visualize selection through a post-process effect on the custom depth
	if (mesh)
	{
		mesh->SetRenderCustomDepth(visibility);
	}

	// show the children
	USceneComponent * sceneComponent = this->GetRootComponent();
	TArray<USceneComponent*> children;

	sceneComponent->GetChildrenComponents(true, children);

	for (USceneComponent * child : children)
	{
		if (UStaticMeshComponent * childMesh = Cast<UStaticMeshComponent>(child))
			childMesh->SetRenderCustomDepth(visibility);
	}
}

void AElementActor::showSelectionOutline()
{
	_setSelectionOutlineVisibility(true);
}

void AElementActor::hideSelectionOutline()
{
	_setSelectionOutlineVisibility(false);
}