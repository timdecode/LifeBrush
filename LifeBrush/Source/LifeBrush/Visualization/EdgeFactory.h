// Copyright 2016 Timothy Davison, all rights reserved. Written for Ship Editor.

#pragma once

#include "Components/ActorComponent.h"


#include <utility>
#include <unordered_map>
#include <vector>
#include <array>

#include "Utility.h"
#include "MeshFactory.h"

#include "RuntimeMeshComponent.h"

#include "EdgeFactory.generated.h"

UCLASS(ClassGroup = (Custom), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API ULineFactory : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	class UMaterialInterface * material;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvBottomY = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvTopY = 1.0f;

	// This will add additional scaling of the texture between two nodes of an edge. If the section specified by uvBottomLeft 
	// and uvTopRight does not horizontally span the texture, there will be artifacts.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvXScale = 1.0f;

	UPROPERTY(Instanced, EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	URuntimeMeshComponent * _runtimeMeshComponent = nullptr;

public:

	virtual void OnComponentCreated()
	{
		Super::OnComponentCreated();

		_runtimeMeshComponent = NewObject<URuntimeMeshComponent>(GetOwner());

		_runtimeMeshComponent->AttachToComponent(this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		_runtimeMeshComponent->bCastDynamicShadow = true;

		_runtimeMeshComponent->SetCollisionProfileName(TEXT("NoCollision"));

		_runtimeMeshComponent->RegisterComponent();
	}

	virtual void OnComponentDestroyed(bool bDestroyingHierarchy)
	{
		Super::OnComponentDestroyed(bDestroyingHierarchy);

		if (_runtimeMeshComponent)
			_runtimeMeshComponent->DestroyComponent();
	}

	struct LineElement
	{
		enum class SegmentType
		{
			Start,
			Segment,
			End
		};

		FVector point;

		float radius;

		SegmentType type;

		LineElement(FVector point, float radius) : point(point), radius(radius), type(SegmentType::Segment) {}
		LineElement(FVector point, float radius, SegmentType type) : point(point), radius(radius), type(type) {}
	};

	virtual void processElements(TArray<LineElement> lineSegments, FVector upNormal, int32 section, UMaterialInterface * material_in = nullptr);

	virtual void clearSection(int32 section)
	{
		_runtimeMeshComponent->ClearMeshSection(section);
	}

	virtual int32 numSections()
	{
		return _runtimeMeshComponent->GetNumSections();
	}

protected:
	// Unreal seems to use a flipped uv coordinate system? bottom-left is top-right.
	std::array<FVector2D, 4> getUVs()
	{
		// This isn't a static variable, because we need to access uvBottomY and uvTopY
		std::array<FVector2D, 4> uvs =
		{
			FVector2D(1.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvTopY),
			FVector2D(1.0f, 1.0f - uvTopY)
		};

		return uvs;
	}

	bool _didInit = false;
};

struct FColoredLineBuilder
{
public:
	struct LineElement
	{
		enum class SegmentType
		{
			Begin,
			Middle,
			End
		};

		FVector point;
		FColor color;

		float radius;

		SegmentType type;

		LineElement(FVector point, float radius) : point(point), radius(radius), type(SegmentType::Middle) {}
		LineElement(FVector point, float radius, SegmentType type) : point(point), radius(radius), type(type) {}
	};

	size_t numCircleComponents = 4;

	TArray<LineElement> lineSegments;

	void begin(FVector& position, float radius);
	void addPoint(FVector& position, float radius);
	void end(FVector& position, float radius);

	void clear();


	template<typename Func>
	void processWithFunctionOld(Func func)
	{
		using SegmentType = LineElement::SegmentType;

		const auto uvs = getUVs(0.0f, 1.0f);

		// pre-compute the unit circle
		TArray<FVector> unitCircle;
		{
			unitCircle.SetNum(numCircleComponents);

			float piStep = M_PI * 2.0 / float(numCircleComponents);
			for (int i = 0; i < numCircleComponents; ++i)
			{
				float t = float(i) * piStep;

				float x = FMath::Cos(t);
				float z = FMath::Sin(t);

				unitCircle[i] = FVector(x, 0.0f, z);
			}
		}

		// find the vertex positions
		{
			const FVector up = FVector::UpVector;

			const int32 m = numCircleComponents * 2;

			int32 pi = 0;

			TArray<FVector> circleVertices;
			circleVertices.SetNum(m);

			FQuat lastQuat = FQuat::Identity;
			FVector lastPosition = FVector::ZeroVector;

			for (int32 i = 1; i < lineSegments.Num(); ++i)
			{
				auto& element_last = lineSegments[i - 1];
				auto& element_cur = lineSegments[i + 0];

				const FVector dirCur = (element_cur.point - element_last.point).GetSafeNormal();

				if (element_cur.type == SegmentType::Begin)
				{
					continue;
				}

				if (element_last.type == SegmentType::Begin)
				{
					lastQuat = FQuat::FindBetween(FVector::RightVector, dirCur);
				}

				// Find rotation_a
				FQuat quatA = lastQuat;

				// Find rotation_b
				FQuat quatB = FQuat::Identity;

				if (element_cur.type == SegmentType::Middle)
				{
					quatB = FQuat::FindBetweenNormals(FVector::RightVector, dirCur);
				}
				else if (element_cur.type == SegmentType::End)
				{
					quatB = FQuat::FindBetweenNormals(FVector::RightVector, dirCur);
				}

				lastQuat = quatB;

				const FVector& a = element_last.point;
				const FVector& b = element_cur.point;

				// compute the transformed unit circle vertices
				for (int ci = 0; ci < m; ci += 2)
				{
					FVector va = quatA.RotateVector(unitCircle[ci / 2]);
					FVector vb = quatB.RotateVector(unitCircle[ci / 2]);

					circleVertices[ci + 1] = va * element_cur.radius + a;
					circleVertices[ci] = vb * element_cur.radius + b;
				}

				const LineElement& evenElement = element_last;
				const LineElement& oddElement = element_cur;


				for (int j = 0; j < m; j += 2)
				{

					func(circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0],
						uvs[(j + 3) % m], uvs[j + 1], uvs[j + 0],
						oddElement.color, oddElement.color, evenElement.color);

					func(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0],
						uvs[(j + 2) % m], uvs[(j + 3) % m], uvs[j + 0],
						evenElement.color, oddElement.color, evenElement.color
					);
				}
			}
		}
	}

	// Frenet-Serret frames: https://en.wikipedia.org/wiki/Frenet%E2%80%93Serret_formulas
	template<typename Func>
	void processWithFunction(Func func)
	{
		using SegmentType = LineElement::SegmentType;

		const auto uvs = getUVs(0.0f, 1.0f);

		// pre-compute the unit circle
		TArray<FVector> unitCircle;
		{
			unitCircle.SetNum(numCircleComponents);

			float piStep = M_PI * 2.0 / float(numCircleComponents);
			for (int i = 0; i < numCircleComponents; ++i)
			{
				float t = float(i) * piStep;

				float x = FMath::Cos(t);
				float z = FMath::Sin(t);

				unitCircle[i] = FVector(0.0, x, z);
			}
		}

		// find the vertex positions
		{
			const FVector up = FVector::UpVector;

			const int32 m = numCircleComponents * 2;

			int32 pi = 0;

			TArray<FVector> circleVertices;
			circleVertices.SetNum(m);

			FQuat lastQuat = FQuat::Identity;
			FVector dirLast;

			FVector tangent;
			FVector normal;
			FVector binormal;

			// draw segments
			for (int32 i = 0; i < lineSegments.Num() - 1; ++i)
			{
				auto& element_cur = lineSegments[i + 0];
				auto& element_next = lineSegments[i + 1];

				FVector dirCur = (element_next.point - element_cur.point).GetSafeNormal();

				// protect against no direction
				if (dirCur.Equals(FVector::ZeroVector))
					dirCur = dirLast;

				FQuat curQuat = FQuat::Identity;

				if (element_cur.type == SegmentType::Begin)
				{
					tangent = dirCur;

					if (tangent.Equals(FVector::ForwardVector))
						normal = FVector::UpVector;
					else
						normal = FVector::CrossProduct(FVector::ForwardVector, tangent);

					binormal = FVector::CrossProduct(tangent, normal);

					curQuat = FRotationMatrix::MakeFromXY(tangent, binormal).ToQuat();
					lastQuat = curQuat;
				}
				else if (element_cur.type == SegmentType::Middle)
				{
					tangent = dirCur;

					normal = FVector::CrossProduct(binormal, tangent);
					binormal = FVector::CrossProduct(tangent, normal);

					curQuat = FRotationMatrix::MakeFromXY(tangent, binormal).ToQuat();
				}
				else if (element_cur.type == SegmentType::End)
					continue;


				auto drawSegment = [&](LineElement& e_a, LineElement& e_b, FQuat quat_a, FQuat quat_b) {
					const FVector& a = e_a.point;
					const FVector& b = e_b.point;

					// compute the transformed unit circle vertices
					for (int ci = 0; ci < m; ci += 2)
					{
						FVector va = quat_a.RotateVector(unitCircle[ci / 2]);
						FVector vb = quat_b.RotateVector(unitCircle[ci / 2]);

						circleVertices[ci + 1] = va * element_cur.radius + a;
						circleVertices[ci] = vb * element_cur.radius + b;
					}

					const LineElement& evenElement = e_a;
					const LineElement& oddElement = e_b;

					for (int j = 0; j < m; j += 2)
					{

						func(circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0],
							uvs[(j + 3) % m], uvs[j + 1], uvs[j + 0],
							oddElement.color, oddElement.color, evenElement.color);

						func(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0],
							uvs[(j + 2) % m], uvs[(j + 3) % m], uvs[j + 0],
							evenElement.color, oddElement.color, evenElement.color
						);
					}
				};

				if (i >= 1)
				{
					auto& element_last = lineSegments[i - 1];

					if( element_last.type != SegmentType::End && element_last.type != SegmentType::Begin )
						drawSegment(element_last, element_cur, lastQuat, curQuat);
				}

				drawSegment(element_cur, element_next, curQuat, curQuat);


				lastQuat = curQuat;
				dirLast = dirCur;
			}

			// draw caps
			for (int32 i = 0; i < lineSegments.Num() - 1; ++i)
			{
				auto& element_cur = lineSegments[i + 0];
				auto& element_next = lineSegments[i + 1];

				FVector dirCur = (element_next.point - element_cur.point).GetSafeNormal();

				// protect against no direction
				if (dirCur.Equals(FVector::ZeroVector))
					dirCur = dirLast;

				FQuat curQuat = FQuat::Identity;

				if (element_cur.type == SegmentType::Begin)
				{
					tangent = dirCur;

					if (tangent.Equals(FVector::ForwardVector))
						normal = FVector::UpVector;
					else
						normal = FVector::CrossProduct(FVector::ForwardVector, tangent);

					binormal = FVector::CrossProduct(tangent, normal);

					curQuat = FRotationMatrix::MakeFromXY(tangent, binormal).ToQuat();
					lastQuat = curQuat;
				}
				else if (element_cur.type == SegmentType::Middle)
				{
					tangent = dirCur;

					normal = FVector::CrossProduct(binormal, tangent);
					binormal = FVector::CrossProduct(tangent, normal);

					curQuat = FRotationMatrix::MakeFromXY(tangent, binormal).ToQuat();
				}
				else if (element_cur.type == SegmentType::End)
					continue;



				//auto drawCap = [&](LineElement& e_a, FQuat quat_a, bool flip = false) {
				//	const FVector& a = e_a.point;

				//	// compute the transformed unit circle vertices
				//	for (int ci = 0; ci < m; ci += 2)
				//	{
				//		FVector va = quat_a.RotateVector(unitCircle[ci / 2]);
				//		
				//		int ia = !flip ? 1 : 0;
				//		int ib = !flip ? 0 : 1;

				//		circleVertices[ci + ia] = a;
				//		circleVertices[ci + ib] = va * e_a.radius + a;
				//	}

				//	const LineElement& evenElement = e_a;
				//	const LineElement& oddElement = e_a;

				//	for (int j = 0; j < m; j += 2)
				//	{
				//		func(circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0],
				//			uvs[(j + 3) % m], uvs[j + 1], uvs[j + 0],
				//			oddElement.color, oddElement.color, evenElement.color);

				//		func(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0],
				//			uvs[(j + 2) % m], uvs[(j + 3) % m], uvs[j + 0],
				//			evenElement.color, oddElement.color, evenElement.color
				//		);
				//	}
				//};

				auto drawCap = [&](FVector point_a, float radius_a, FColor color_a, FVector point_b, float radius_b, FColor color_b, FQuat quat, bool flip = false) {
					// compute the transformed unit circle vertices
					for (int ci = 0; ci < m; ci += 2)
					{
						FVector va = quat.RotateVector(unitCircle[ci / 2]);

						int ia = !flip ? 1 : 0;
						int ib = !flip ? 0 : 1;

						circleVertices[ci + ia] = va * radius_a + point_a;
						circleVertices[ci + ib] = va * radius_b + point_b;
					}

					for (int j = 0; j < m; j += 2)
					{
						func(circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0],
							uvs[(j + 3) % m], uvs[j + 1], uvs[j + 0],
							color_b, color_b, color_a);

						func(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0],
							uvs[(j + 2) % m], uvs[(j + 3) % m], uvs[j + 0],
							color_a, color_b, color_a
						);
					}
				};

				auto drawHemisphere = [&](LineElement& element, bool flip = false)
				{
					FVector forwardDir = curQuat  * FVector::ForwardVector * (flip ? 1.0f : -1.0f);

					FVector p = element.point;
					float r = element.radius;

					const int numLats = 2;
					for (int lat = 0; lat < numLats; ++lat)
					{
						float frac_a = (float(lat + 1) / float(numLats));
						float frac_b = (float(lat) / float(numLats));

						FVector a = p + forwardDir * frac_a;
						FVector b = p + forwardDir * frac_b;

						float r_a = FMath::Sqrt(r * r - FMath::Pow(frac_a * r, 2.0f));
						float r_b = FMath::Sqrt(r * r - FMath::Pow(frac_b * r, 2.0f));

						drawCap(a, r_a, element.color, b, r_b, element.color, curQuat, flip);
					}
				};

				if (element_cur.type == SegmentType::Begin)
				{
					drawHemisphere(element_cur, false);
				}

				if (element_next.type == SegmentType::End)
				{
					drawHemisphere(element_next, true);
				}

				lastQuat = curQuat;
				dirLast = dirCur;
			}

		}
	}

	void createToBuilder(TSharedRef<FRuntimeMeshBuilder> builderRef)
	{
		FRuntimeMeshBuilder& builder = builderRef.Get();

		FRuntimeMeshVertexSimpleFactory factory(builder);

		processWithFunction([&](
			const FVector& v0, const FVector& v1, const FVector& v2,
			const FVector2D& u0, const FVector2D &u1, const FVector2D& u2,
			const FColor& c0, const FColor& c1, const FColor& c2)
		{
			factory.pushTriangle(
				v0, v1, v2,
				u0, u1, u2,
				c0, c1, c2);
		});
	}

	void updateBuilder(TSharedRef<FRuntimeMeshBuilder> builderRef)
	{
		FRuntimeMeshBuilder& builder = builderRef.Get();

		FRuntimeMeshVertexSimpleFactory factory(builder);

		int32 index = 0;

		processWithFunction([&](
			const FVector& v0, const FVector& v1, const FVector& v2,
			const FVector2D& u0, const FVector2D &u1, const FVector2D& u2,
			const FColor& c0, const FColor& c1, const FColor& c2)
		{
			factory.updateTriangle(
				index,
				v0, v1, v2,
				u0, u1, u2,
				c0, c1, c2);

			index += 3;
		});
	}

	void updateRuntimeMeshData(FRuntimeMeshScopedUpdater& meshData)
	{

		int32 index = 0;

		processWithFunction([&](
			const FVector& v0, const FVector& v1, const FVector& v2,
			const FVector2D& u0, const FVector2D &u1, const FVector2D& u2,
			const FColor& c0, const FColor& c1, const FColor& c2)
		{
			FVector n3 = FVector::CrossProduct(v2 - v0, v1 - v0);
			n3.Normalize();

			FVector4 n4 = FVector4(n3, 1.0f);

			FVector tangentV = FVector::CrossProduct(n3, v2 - v0);
			tangentV.Normalize();

			int32 i0 = index + 0;
			int32 i1 = index + 1;
			int32 i2 = index + 2;

			meshData.SetPosition(i0, v0);
			meshData.SetPosition(i1, v1);
			meshData.SetPosition(i2, v2);

			meshData.SetNormal(i0, n4);
			meshData.SetNormal(i1, n4);
			meshData.SetNormal(i2, n4);

			meshData.SetTangent(i0, tangentV);
			meshData.SetTangent(i1, tangentV);
			meshData.SetTangent(i2, tangentV);

			index += 3;
		});
	}

	/*
	\param uvXScale This will add additional scaling of the texture between two nodes of an edge. If the section specified by uvBottomLeft 
	  and uvTopRight does not horizontally span the texture, there will be artifacts.
	*/
	ColoredQuadFactory createQuadFactory(float uvBottomY = 0.0f, float uvTopY = 1.0f, float uvXScale = 1.0f);

	void appendToQuadFactory(ColoredQuadFactory& factory, float uvBottomY = 0.0f, float uvTopY = 1.0f, float uvXScale = 1.0f);

protected:
	// Unreal seems to use a flipped uv coordinate system? bottom-left is top-right.
	std::array<FVector2D, 4> getUVs(float uvBottomY, float uvTopY)
	{
		// This isn't a static variable, because we need to access uvBottomY and uvTopY
		std::array<FVector2D, 4> uvs =
		{
			FVector2D(1.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvTopY),
			FVector2D(1.0f, 1.0f - uvTopY)
		};

		return uvs;
	}
};

UCLASS(ClassGroup = (Custom), meta = (BlueprintSpawnableComponent))
class LIFEBRUSH_API UColoredLineFactory : public UActorComponent
{
	GENERATED_BODY()

public:
	// The material should use vertex colors.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	class UMaterialInterface * material;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvBottomY = 0.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvTopY = 1.0f;

	// This will add additional scaling of the texture between two nodes of an edge. If the section specified by uvBottomLeft 
	// and uvTopRight does not horizontally span the texture, there will be artifacts.
	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	float uvXScale = 1.0f;

	UPROPERTY(Instanced, EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	URuntimeMeshComponent * _runtimeMeshComponent = nullptr;

public:

	virtual void OnComponentCreated()
	{
		Super::OnComponentCreated();

		_runtimeMeshComponent = NewObject<URuntimeMeshComponent>(GetOwner());

		_runtimeMeshComponent->AttachToComponent(this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		_runtimeMeshComponent->bCastDynamicShadow = true;

		_runtimeMeshComponent->SetCollisionProfileName(TEXT("NoCollision"));

		_runtimeMeshComponent->RegisterComponent();
	}

	virtual void OnComponentDestroyed(bool bDestroyingHierarchy)
	{
		Super::OnComponentDestroyed(bDestroyingHierarchy);

		if (_runtimeMeshComponent)
			_runtimeMeshComponent->DestroyComponent();
	}

	virtual void commitWithFastPathOption(FColoredLineBuilder& builder, int32 section, UMaterialInterface * material_in = nullptr, bool topologyDidChange = true);


	virtual void commitSection(FColoredLineBuilder& builder, int32 section, UMaterialInterface * material_in = nullptr, bool topologyDidChange = true);

	virtual void clearSection(int32 section)
	{
		_runtimeMeshComponent->ClearMeshSection(section);
	}

	virtual int32 numSections()
	{
		return _runtimeMeshComponent->GetNumSections();
	}

	TArray<FRuntimeMeshVertexSimple> _simple(FColoredLineBuilder& builder);

protected:
	TSharedRef<FRuntimeMeshBuilder> meshBuilder = MakeRuntimeMeshBuilder(false, false, 1, true);

	bool _didInit = false;
};

UCLASS( ClassGroup=(Custom), meta=(BlueprintSpawnableComponent) )
class LIFEBRUSH_API UEdgeFactory : public UActorComponent
{
	GENERATED_BODY()

public:
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	class UMaterialInterface * material;
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float uvBottomY = 0.0f;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float uvTopY = 1.0f;

	// This will add additional scaling of the texture between two nodes of an edge. If the section specified by uvBottomLeft 
	// and uvTopRight does not horizontally span the texture, there will be artifacts.
	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	float uvXScale = 1.0f;

	UPROPERTY(EditAnywhere, BlueprintReadWrite, Category = "ShipEditor")
	bool enableCollision = false;

	UPROPERTY( Instanced, EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	URuntimeMeshComponent * _runtimeMeshComponent = nullptr;

	UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
	int32 edgeType = 0;

public:	
	// Sets default values for this component's properties
	UEdgeFactory() {};

	virtual void OnComponentCreated()
	{
		Super::OnComponentCreated();

		_runtimeMeshComponent = NewObject<URuntimeMeshComponent>( GetOwner() );

		_runtimeMeshComponent->AttachToComponent( this->GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );

		_initRuntimeMeshComponent();

		_runtimeMeshComponent->RegisterComponent();
	}

	virtual void OnComponentDestroyed( bool bDestroyingHierarchy )
	{
		Super::OnComponentDestroyed( bDestroyingHierarchy );

		if(_runtimeMeshComponent)
			_runtimeMeshComponent->DestroyComponent();
	}

	void _initRuntimeMeshComponent()
	{
		_runtimeMeshComponent->bCastDynamicShadow = true;

		if( enableCollision )
			_runtimeMeshComponent->SetCollisionProfileName( TEXT( "OverlapAllDynamic" ) );
	}

	struct EdgeParameter
	{
		FVector a;
		FVector b;
		float halfExtents;
		bool isCapped = false;

		EdgeParameter( FVector a_, FVector b_, float halfExtents_ ) : a( a_ ), b( b_ ), halfExtents( halfExtents_ ) {}
	};

	virtual void processEdges( TArray<EdgeParameter> edges, FVector upNormal, int32 section, UMaterialInterface * material_in = nullptr ) {};

	virtual void clearSection( int32 section )
	{
		_runtimeMeshComponent->ClearMeshSection( section );
	}

	virtual int32 numSections()
	{
		return _runtimeMeshComponent->GetNumSections();
	}

protected:
	// Unreal seems to use a flipped uv coordinate system? bottom-left is top-right.
	std::array<FVector2D, 4> getUVs()
	{
		// This isn't a static variable, because we need to access uvBottomY and uvTopY
		std::array<FVector2D, 4> uvs =
		{
			FVector2D(1.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvBottomY),
			FVector2D(0.0f, 1.0f - uvTopY),
			FVector2D(1.0f, 1.0f - uvTopY)
		};

		return uvs;
	}

	bool _didInit = false;
};

struct FBlockEdgeFactory
{
public:
	static QuadFactory quadFactoryForEdges(TArray<UEdgeFactory::EdgeParameter>& edges, FVector upNormal)
	{
		QuadFactory quadFactory;


		FVector up = upNormal;

		// fight depth for edges that share the same nodes by giving the generated meshes very slightly different half extents
		std::unordered_map<FVector, int> depthFighter;

		UEdgeFactory::EdgeParameter * lastPair = nullptr;

		for (auto& pair : edges)
		{
			FVector& a = pair.a;
			FVector& b = pair.b;
			float halfExtents = pair.halfExtents;

			FVector dir = b - a;
			float length = dir.Size();
			dir /= length;

			FVector2D uvLengthX = FVector2D(length / (halfExtents * 2.0f), 1.0f);
			FVector2D uvLenghtY = FVector2D(1.0f, length / (halfExtents * 2.0f));

			FVector side = FVector::CrossProduct(dir, up).GetSafeNormal();
			FVector relativeUp = FVector::CrossProduct(side, dir).GetSafeNormal();

			if (dir.Equals(up, 0.0001f))
			{
				side = FVector::RightVector;
				relativeUp = FVector::ForwardVector;
			}
			else if ((-dir).Equals(up, 0.0001f))
			{
				side = -FVector::RightVector; // left-vector
				relativeUp = FVector::ForwardVector;
			}

			int& aCount_ = depthFighter[a];
			int& bCount_ = depthFighter[b];

			aCount_++;
			bCount_++;

			int aCount = aCount_ - 1;
			int bCount = bCount_ - 1;

			float aHalfExtents = halfExtents * (1.0f + float(aCount) * 0.001f);
			float bHalfExtents = halfExtents * (1.0f + float(bCount) * 0.001f);

			FVector za = relativeUp * aHalfExtents;
			FVector xa = side * aHalfExtents;

			FVector zb = relativeUp * bHalfExtents;
			FVector xb = side * bHalfExtents;

			//    v5 --- v6
			//   /|     /|
			// v1 --- v2 |     y 
			// |  v4 -|- v7    | z
			// | /    | /      |/
			// v0 --- v3       o --- x
			FVector vertices[] = {
				a - xa - za, // 0
				a - xa + za, // 1
				a + xa + za, // 2
				a + xa - za, // 3

				b - xb - zb, // 4
				b - xb + zb, // 5
				b + xb + zb, // 6
				b + xb - zb  // 7
			};

			FVector2D uvs[] =
			{
				FVector2D(0,0),
				FVector2D(0,1),
				FVector2D(1,1),
				FVector2D(1,0)
			};

			if (pair.isCapped)
			{
				quadFactory.pushQuad(
					vertices[0],
					vertices[1],
					vertices[2],
					vertices[3],
					uvs[0],
					uvs[1],
					uvs[2],
					uvs[3]);

				quadFactory.pushQuad(
					vertices[4],
					vertices[7],
					vertices[6],
					vertices[5],
					uvs[0],
					uvs[1],
					uvs[2],
					uvs[3]);
			}

			quadFactory.pushQuad(
				vertices[1],
				vertices[0],
				vertices[4],
				vertices[5],
				uvs[0],
				uvs[1],
				uvs[2] * uvLengthX,
				uvs[3] * uvLengthX);

			quadFactory.pushQuad(
				vertices[3],
				vertices[2],
				vertices[6],
				vertices[7],
				uvs[0],
				uvs[1],
				uvs[2] * uvLengthX,
				uvs[3] * uvLengthX);

			quadFactory.pushQuad(
				vertices[2],
				vertices[1],
				vertices[5],
				vertices[6],
				uvs[0],
				uvs[1],
				uvs[2] * uvLengthX,
				uvs[3] * uvLengthX);

			quadFactory.pushQuad(
				vertices[0],
				vertices[3],
				vertices[7],
				vertices[4],
				uvs[0],
				uvs[1],
				uvs[2] * uvLengthX,
				uvs[3] * uvLengthX);

			if (lastPair)
			{

			}

			lastPair = &pair;
		}

		return quadFactory;
	}
};

UCLASS( ClassGroup = (Custom), meta = (BlueprintSpawnableComponent) )
class LIFEBRUSH_API UBlockyEdgeFactory : public UEdgeFactory
{
	GENERATED_BODY()

public:
	UBlockyEdgeFactory() {}

	// Positions are in local space.
	virtual void processEdges( TArray<EdgeParameter> edges, FVector upNormal, int32 section, UMaterialInterface * material_in = nullptr )
	{
		_processEdges( edges, upNormal, section, false, material_in );
	}

	virtual void _processEdges( TArray<EdgeParameter>& edges, FVector& upNormal, int32 meshSection, bool buildPhysics = false, UMaterialInterface * material_in = nullptr )
	{
		QuadFactory quadFactory = FBlockEdgeFactory::quadFactoryForEdges(edges, upNormal);

		_runtimeMeshComponent->ClearMeshSection( meshSection );

		_runtimeMeshComponent->CreateMeshSection( meshSection, quadFactory.vertices, quadFactory.indices, quadFactory.normals, quadFactory.uvs /* uvs*/, TArray<FColor>() /* colors */, quadFactory.tangents, enableCollision, EUpdateFrequency::Frequent );

		_runtimeMeshComponent->SetMeshSectionCastsShadow( meshSection, true );
		_runtimeMeshComponent->SetMeshSectionCollisionEnabled( meshSection, enableCollision );
		_runtimeMeshComponent->SetMaterial( meshSection, material_in ? material_in : material );
	}
};

UCLASS( ClassGroup = (Custom), meta = (BlueprintSpawnableComponent) )
class LIFEBRUSH_API UBeveledEdgeFactory : public UEdgeFactory
{
	GENERATED_BODY()

public:
	UBeveledEdgeFactory() {}

	virtual void processEdges( TArray<EdgeParameter> edges, FVector upNormal, int32 section, UMaterialInterface * material_in = nullptr )
	{
		_processEdges( edges, upNormal, section, true, material_in );
	}

	virtual void _processEdges( TArray<EdgeParameter>& edges, FVector& upNormal, int32 meshSection, bool buildPhysics = false, UMaterialInterface * material_in = nullptr )
	{
		SmoothMeshFactory meshFactor;

		float bevelWidth = 0.5f;

		const int n = 16;

		FVector up = FVector::UpVector;

		// fight depth for edges that share the same nodes by giving the generated meshes very slightly different half extents
		std::unordered_map<FVector, int> depthFighter;

		TArray< TArray<FVector> > convexMeshes;

		auto uvs = getUVs();

		for(auto& pair : edges)
		{
			FVector& a = pair.a;
			FVector& b = pair.b;
			float halfExtents = pair.halfExtents;
			FVector dir = a - b;
			float length = dir.Size();
			dir /= length;

			int& aCount_ = depthFighter[a];
			int& bCount_ = depthFighter[b];

			aCount_++;
			bCount_++;

			int aCount = aCount_ - 1;
			int bCount = bCount_ - 1;

			FVector2D uvLengthX = FVector2D( length / (halfExtents * 2.0f), 1.0f );

			FVector side = FVector::CrossProduct( dir, up ).GetSafeNormal();
			FVector relativeUp = FVector::CrossProduct( side, dir ).GetSafeNormal();

			if(dir.Equals( up, 0.0001f ))
			{
				side = FVector::RightVector;
				relativeUp = FVector::ForwardVector;
			}
			else if((-dir).Equals( up, 0.0001f ))
			{
				side = -FVector::RightVector; // left-vector
				relativeUp = FVector::ForwardVector;
			}

			float aHalfExtents = halfExtents * (1.0f + float( aCount ) * 0.001f);
			float bHalfExtents = halfExtents * (1.0f + float( bCount ) * 0.001f);


			FVector za = relativeUp * aHalfExtents;
			FVector xa = side * aHalfExtents;

			FVector zb = relativeUp * bHalfExtents;
			FVector xb = side * bHalfExtents;

			// bevel x and z offsets

			FVector bev_xa = side * (aHalfExtents - bevelWidth);
			FVector bev_za = relativeUp * (aHalfExtents - bevelWidth);

			FVector bev_xb = side * (bHalfExtents - bevelWidth);
			FVector bev_zb = relativeUp * (bHalfExtents - bevelWidth);

			FVector vertices[] =
			{
				a - bev_xa - za, // 0
				b - bev_xb - zb, // 1
				a + bev_xa - za, // 2
				b + bev_xb - zb, // 3
				a + xa - bev_za, // 4
				b + xb - bev_zb, // 5
				a + xa + bev_za, // 6
				b + xb + bev_zb, // 6
				a + bev_xa + za, // 8
				b + bev_xb + zb, // 8
				a - bev_xa + za, // 10
				b - bev_xb + zb, // 11
				a - xa + bev_za, // 12
				b - xb + bev_zb, // 13
				a - xa - bev_za, // 14
				b - xb - bev_zb  // 15
			};

			// sides
			for( int i = 0; i < n; i += 2)
			{
				meshFactor.pushQuad(
					vertices[i + 0],
					vertices[i + 1],
					vertices[(i + 3) % n],
					vertices[(i + 2) % n],
					uvs[0] * uvLengthX, // why this guy? cus we're doing (1,1) - uv when creating the uvs array
					uvs[1],
					uvs[2],
					uvs[3] * uvLengthX );
			}

			// caps
			for(int i = 0; i < n; i += 2)
			{
				meshFactor.pushTriangle(
					a,
					vertices[i],
					vertices[(i + 2) % n],
					dir,
					dir,
					dir
				);

				meshFactor.pushTriangle(
					b,
					vertices[(i + 3) % n],
					vertices[i + 1],
					-dir,
					-dir,
					-dir
				);
			}

			// add a convex mesh
			FVector z = relativeUp * halfExtents;
			FVector x = side * halfExtents;

			//    v5 --- v6
			//   /|     /|
			// v1 --- v2 |
			// |  v4 -|- v7
			// | /    | /
			// v0 --- v3
			TArray<FVector> convexMesh;
			convexMesh.AddUninitialized( 8 ); 
			
			convexMesh[0] = a - x - z; // 0
			convexMesh[1] = a - x + z; // 1
			convexMesh[2] = a + x + z; // 2
			convexMesh[3] = a + x - z; // 3

			convexMesh[4] = b - x - z; // 4
			convexMesh[5] = b - x + z; // 5
			convexMesh[6] = b + x + z; // 6
			convexMesh[7] = b + x - z; // 7

			convexMeshes.Add( convexMesh );
		}

		_runtimeMeshComponent->CreateMeshSection( meshSection, meshFactor.vertices, meshFactor.indices, meshFactor.normals, meshFactor.uvs /* uvs*/, TArray<FColor>() /* colors */, meshFactor.tangents, enableCollision /* create collision */, EUpdateFrequency::Infrequent );
		_runtimeMeshComponent->SetMeshSectionCastsShadow( meshSection,true);
		_runtimeMeshComponent->SetMeshSectionCollisionEnabled( meshSection, enableCollision);
		_runtimeMeshComponent->SetMaterial( meshSection, material_in ? material_in : material );
		_runtimeMeshComponent->bMultiBodyOverlap = true;

		if(!buildPhysics)
			return;

		URuntimeMeshComponent * physicsMesh = Cast<URuntimeMeshComponent>( GetOwner()->GetRootComponent() );
		if(!physicsMesh)
			return;
		
		physicsMesh->SetCollisionConvexMeshes( convexMeshes );

	}
};

UCLASS( ClassGroup = (Custom), meta = (BlueprintSpawnableComponent) )
class LIFEBRUSH_API UPipeEdgeFactory : public UEdgeFactory
{
	GENERATED_BODY()

public:
	UPipeEdgeFactory() {}

	virtual void processEdges( TArray<EdgeParameter> edges, FVector upNormal, int32 section, UMaterialInterface * material_in = nullptr )
	{
		_processEdges( edges, upNormal, section, false, material_in );
	}

	virtual void _processEdges( TArray<EdgeParameter>& edges, FVector& upNormal, int32 meshSection, bool buildPhysics = false, UMaterialInterface * material_in = nullptr )
	{
		SmoothMeshFactory meshFactory;

		float bevelWidth = 0.5f;

		const int n = 16;
		const int m = n * 2;

		FVector up = FVector::UpVector;

		// fight depth for edges that share the same nodes by giving the generated meshes very slightly different half extents
		std::unordered_map<FVector, int> depthFighter;

		TArray< TArray<FVector> > convexMeshes;

		// unreal seems to use a flipped uv coordinate system? bottom-left is top-right
		std::array<FVector2D, 4> uvs = getUVs();

		// build our unit circle
		FVector unitCircle[n];

		float piStep = M_PI * 2 / float( n );
		for(int i = 0; i < n; ++i)
		{
			float t = float( i ) * piStep;

			float x = FMath::Cos( t );
			float z = FMath::Sin( t );

			unitCircle[i] = FVector( x, 0.0f, z );
		}


		for(auto& pair : edges)
		{
			FVector& a = pair.a;
			FVector& b = pair.b;
			float halfExtents = pair.halfExtents;
			FVector dir = a - b;
			float length = dir.Size();
			dir /= length;

			int& aCount_ = depthFighter[a];
			int& bCount_ = depthFighter[b];

			aCount_++;
			bCount_++;

			int aCount = aCount_ - 1;
			int bCount = bCount_ - 1;

			FVector2D uvLengthX = FVector2D( length * uvXScale, 1.0f ) ;

			FVector side = FVector::CrossProduct( dir, up ).GetSafeNormal();
			FVector relativeUp = FVector::CrossProduct( side, dir ).GetSafeNormal();

			if(dir.Equals( up, 0.0001f ))
			{
				side = FVector::RightVector;
				relativeUp = FVector::ForwardVector;
			}
			else if((-dir).Equals( up, 0.0001f ))
			{
				side = -FVector::RightVector; // left-vector
				relativeUp = FVector::ForwardVector;
			}

			FQuat quat = FQuat::FindBetweenNormals( FVector::RightVector, dir );

			FVector vertices[m];
			FVector normals[m];
			for(int i = 0; i < m; i += 2) 
			{
				FVector v = quat.RotateVector( unitCircle[i / 2] );  

				vertices[i] = v * halfExtents + a;
				vertices[i + 1] = v * halfExtents + b;

				normals[i] = v;
				normals[i + 1] = v;
			}

			for(int i = 0; i < m; i += 2)
			{
				// v0,v1,v2
				// v0,v2,v3
				meshFactory.pushQuad(
					vertices[i + 0],
					vertices[i + 1],
					vertices[(i + 3) % m],
					vertices[(i + 2) % m],
					normals[i + 0],
					normals[i + 1],
					normals[(i + 3) % m],
					normals[(i + 2) % m],
					uvs[0] * uvLengthX,
					uvs[1],
					uvs[2],
					uvs[3] * uvLengthX
				);
			}

			// caps
			for(int i = 0; i < m; i += 2)
			{
				meshFactory.pushTriangle(
					a,
					vertices[i],
					vertices[(i + 2) % m],
					dir,
					dir,
					dir 
				);

				meshFactory.pushTriangle(
					b,
					vertices[(i + 3) % m],
					vertices[i + 1],
					-dir,
					-dir,
					-dir
				);
			}

			// add a convex mesh
			FVector z = relativeUp * halfExtents;
			FVector x = side * halfExtents;

			//    v5 --- v6
			//   /|     /|
			// v1 --- v2 |     y 
			// |  v4 -|- v7    | z
			// | /    | /      |/
			// v0 --- v3       o --- x
			TArray<FVector> convexMesh;
			convexMesh.AddUninitialized( 8 );

			convexMesh[0] = a - x - z; // 0
			convexMesh[1] = a - x + z; // 1
			convexMesh[2] = a + x + z; // 2
			convexMesh[3] = a + x - z; // 3

			convexMesh[4] = b - x - z; // 4
			convexMesh[5] = b - x + z; // 5
			convexMesh[6] = b + x + z; // 6
			convexMesh[7] = b + x - z; // 7

			convexMeshes.Add( convexMesh );
		}

		_runtimeMeshComponent->CreateMeshSection( 
			meshSection, 
			meshFactory.vertices, 
			meshFactory.indices, 
			meshFactory.normals, 
			meshFactory.uvs /* uvs*/, 
			TArray<FColor>() /* colors */, 
			meshFactory.tangents, 
			enableCollision /* create collision */,
			EUpdateFrequency::Infrequent 
		);

		_runtimeMeshComponent->SetMeshSectionCastsShadow( meshSection, true );
		_runtimeMeshComponent->SetMeshSectionCollisionEnabled( meshSection, enableCollision);
		_runtimeMeshComponent->SetMaterial( meshSection, material_in ? material_in : material );
		_runtimeMeshComponent->bMultiBodyOverlap = true;

		if(!buildPhysics)
			return;

		URuntimeMeshComponent * physicsMesh = Cast<URuntimeMeshComponent>( GetOwner()->GetRootComponent() );
		if(!physicsMesh)
			return;

		physicsMesh->SetCollisionConvexMeshes( convexMeshes );
	}

private:
	void _initRuntimeMeshComponent()
	{
		_runtimeMeshComponent->SetShouldSerializeMeshData(false);
		_runtimeMeshComponent->bCastDynamicShadow = true;
		_runtimeMeshComponent->SetCollisionProfileName( TEXT( "OverlapAllDynamic" ) );
	}
};