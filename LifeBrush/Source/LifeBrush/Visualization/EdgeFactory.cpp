// Copyright 2018 Timothy Davison, all rights reserved. Written for Ship Editor.

#include "LifeBrush.h"
#include "EdgeFactory.h"


void ULineFactory::processElements(TArray<LineElement> lineSegments, FVector upNormal, int32 section, UMaterialInterface * material_in /*= nullptr*/)
{
	using SegmentType = LineElement::SegmentType;

	if (lineSegments.Num() < 2)
		return;
	
	QuadFactory quadFactory;

	const int32 numCircleComponents = 4;

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

		for (int32 i = 1; i < lineSegments.Num() - 1; ++i)
		{
			auto& element_last = lineSegments[i - 1];
			auto& element_cur = lineSegments[i + 0];
			auto& element_next = lineSegments[i + 1];



			const FVector dirLast = (element_cur.point - element_last.point).GetSafeNormal();
			const FVector dirNext = (element_next.point - element_cur.point).GetSafeNormal();

			const FVector dirCur = (dirNext + dirLast).GetSafeNormal();

			if (element_cur.type == SegmentType::Start)
			{
				lastQuat = FQuat::FindBetweenNormals(FVector::RightVector, dirNext);

				continue;
			}

			// Find rotation_a
			FQuat quatA = lastQuat;

			// Find rotation_b
			FQuat quatB = FQuat::Identity;

			if (element_cur.type == SegmentType::Segment)
			{
				quatB = FQuat::FindBetweenNormals(FVector::RightVector, dirCur);
			}
			else if (element_cur.type == SegmentType::End)
			{
				quatB = FQuat::FindBetweenNormals(FVector::RightVector, dirLast);
			}

			lastQuat = quatB;

			const FVector& a = element_last.point;
			const FVector& b = element_cur.point;

			// compute the transformed unit circle vertices
			for (int ci = 0; ci < m; ci += 2)
			{
				FVector va = quatA.RotateVector(unitCircle[ci / 2]);
				FVector vb = quatB.RotateVector(unitCircle[ci / 2]);

				circleVertices[ci] = va * element_cur.radius + a;
				circleVertices[ci + 1] = vb * element_cur.radius + b;
			}

			for (int j = 0; j < m; j += 2)
			{
				quadFactory.pushTriangle(circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0]);
				quadFactory.pushTriangle(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0]);
			}
		}
	}

	_runtimeMeshComponent->ClearMeshSection(section);

	_runtimeMeshComponent->CreateMeshSection(section, quadFactory.vertices, quadFactory.indices, quadFactory.normals, quadFactory.uvs /* uvs*/, TArray<FColor>() /* colors */, quadFactory.tangents, false, EUpdateFrequency::Frequent);

	_runtimeMeshComponent->SetMeshSectionCastsShadow(section, true);
	_runtimeMeshComponent->SetMeshSectionCollisionEnabled(section, false);
	_runtimeMeshComponent->SetMaterial(section, material_in ? material_in : material);
}

void UColoredLineFactory::commitWithFastPathOption(FColoredLineBuilder& builder, int32 section, UMaterialInterface * material_in /*= nullptr*/, bool topologyDidChange /*= true*/)
{
	if (topologyDidChange || !_runtimeMeshComponent->DoesSectionExist(section) )
	{
		meshBuilder = MakeRuntimeMeshBuilder(false, false, 1, true);
		

		builder.createToBuilder(meshBuilder);

		_runtimeMeshComponent->CreateMeshSection(section, meshBuilder, false, EUpdateFrequency::Frequent, ESectionUpdateFlags::None);

		_runtimeMeshComponent->SetMeshSectionCollisionEnabled(section, false);
		_runtimeMeshComponent->SetMeshSectionCastsShadow(section, true);
		_runtimeMeshComponent->SetMaterial(section, material_in ? material_in : material);
	}
	else
	{
		FRuntimeMeshData& meshData = _runtimeMeshComponent->GetOrCreateRuntimeMesh()->GetRuntimeMeshData().Get();

		auto meshThing = meshData.BeginSectionUpdate(section);

		builder.updateRuntimeMeshData(*meshThing.Get());

		meshThing->Commit(true, true, false, false, false);
	}
}

TArray<FRuntimeMeshVertexSimple> UColoredLineFactory::_simple(FColoredLineBuilder& builder)
{
	TArray<FRuntimeMeshVertexSimple> result;

	return result;
}

void FColoredLineBuilder::begin(FVector& position, float radius)
{
	lineSegments.Emplace(position, radius, LineElement::SegmentType::Begin);
}

void FColoredLineBuilder::addPoint(FVector& position, float radius)
{
	lineSegments.Emplace(position, radius, LineElement::SegmentType::Middle);
}

void FColoredLineBuilder::end(FVector& position, float radius)
{
	lineSegments.Emplace(position, radius, LineElement::SegmentType::End);
}

void FColoredLineBuilder::clear()
{
	lineSegments.Empty();
}


ColoredQuadFactory FColoredLineBuilder::createQuadFactory(float uvBottomY, float uvTopYf, float uvXScale)
{
	ColoredQuadFactory quadFactory;

	if (lineSegments.Num() < 2)
		return quadFactory;

	appendToQuadFactory(quadFactory, uvBottomY, uvTopYf, uvXScale);

	return quadFactory;
}

void FColoredLineBuilder::appendToQuadFactory(ColoredQuadFactory& quadFactory, float uvBottomY /*= 0.0f*/, float uvTopY /*= 1.0f*/, float uvXScale /*= 1.0f*/)
{
	using SegmentType = LineElement::SegmentType;

	const auto uvs = getUVs(uvBottomY, uvTopY);

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

		for (int32 i = 1; i < lineSegments.Num() - 1; ++i)
		{
			auto& element_last = lineSegments[i - 1];
			auto& element_cur = lineSegments[i + 0];
			auto& element_next = lineSegments[i + 1];



			const FVector dirLast = (element_cur.point - element_last.point).GetSafeNormal();
			const FVector dirNext = (element_next.point - element_cur.point).GetSafeNormal();

			const FVector dirCur = (dirNext + dirLast).GetSafeNormal();

			if (element_cur.type == SegmentType::Begin)
			{
				lastQuat = FQuat::FindBetweenNormals(FVector::RightVector, dirNext);

				continue;
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
				quatB = FQuat::FindBetweenNormals(FVector::RightVector, dirLast);
			}

			lastQuat = quatB;

			const FVector& a = element_last.point;
			const FVector& b = element_cur.point;

			// compute the transformed unit circle vertices
			for (int ci = 0; ci < m; ci += 2)
			{
				FVector va = quatA.RotateVector(unitCircle[ci / 2]);
				FVector vb = quatB.RotateVector(unitCircle[ci / 2]);

				circleVertices[ci] = va * element_cur.radius + a;
				circleVertices[ci + 1] = vb * element_cur.radius + b;
			}

			const LineElement& evenElement = element_last;
			const LineElement& oddElement = element_cur;

			for (int j = 0; j < m; j += 2)
			{
				quadFactory.pushTriangle(
					circleVertices[(j + 3) % m], circleVertices[j + 1], circleVertices[j + 0],
					uvs[(j + 3) % m], uvs[j + 1], uvs[j + 0],
					oddElement.color, oddElement.color, evenElement.color
				);

				quadFactory.pushTriangle(circleVertices[(j + 2) % m], circleVertices[(j + 3) % m], circleVertices[j + 0],
					uvs[(j + 2) % m], uvs[(j + 3) % m], uvs[j + 0],
					evenElement.color, oddElement.color, evenElement.color
				);
			}
		}
	}
}
