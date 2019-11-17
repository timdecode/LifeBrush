//
//  Utility.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-07-20.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#pragma once

#include "Eigen/Dense"
#include "Algorithm/tcods/Vector.h"
#include "RuntimeMeshComponent.h"

struct FTimStructBox;
struct FGraph;
struct FGraphNode;

#if  WITH_EDITOR
#include "Editor.h"
#include "EditorViewportClient.h"
#include "LevelEditorViewport.h"
#endif

#include "Utility.generated.h"

USTRUCT()
struct LIFEBRUSH_API FDummy_PleaseCompile
{
GENERATED_BODY()

public:

};

namespace std
{
	template<>
	struct hash<FVector>
	{
		std::size_t operator()( const FVector& k ) const
		{
			return FCrc::MemCrc32(&k, sizeof(float) * 3);
		}
	};

	template<>
	struct less<FVector>
	{
		bool operator()(const FVector& a, const FVector& b) const
		{
			return a.X < b.X ? true : (a.X == b.X ? (a.Y < b.Y ? true : (a.Y == b.Y ? a.Z < b.Z : false)) : false);
		}
	};

	template<>
	struct hash<FIntVector>
	{
		std::size_t operator()(const FIntVector& k) const
		{
			return GetTypeHash(k);
		}
	};

	template<>
	struct less<FIntVector>
	{
		bool operator()(const FIntVector& a, const FIntVector& b) const
		{
			return a.X < b.X ? true : (a.X == b.X ? (a.Y < b.Y ? true : (a.Y == b.Y ? a.Z < b.Z : false)) : false);
		}
	};


}

class LIFEBRUSH_API Utility
{
public:
	static FORCEINLINE FVector scale( const FIntVector& intVector, const FVector& floatVector )
	{
		FVector result;
		
		for(int c = 0; c < 3; ++c)
			result[c] = intVector[c] * floatVector[c];

		return result;
	}

	static FORCEINLINE FVector scale(const FVector& a, const FVector& b)
	{
		FVector result;

		for (int c = 0; c < 3; ++c)
			result[c] = a[c] * b[c];

		return result;
	}

	static FString _nextActorLabelInSequence(FString baseName, UWorld * world)
	{
#if WITH_EDITOR
		auto actors = world->PersistentLevel->Actors;

		int32 maxCounter = 0;

		for (auto actor : actors)
		{
			if (!actor->IsValidLowLevel())
				continue;

			FString name = actor->GetActorLabel();

			int32 found = name.Find(baseName, ESearchCase::CaseSensitive, ESearchDir::FromStart);

			if (found == INDEX_NONE)
				continue;

			FString end(*name + baseName.Len());

			if (!end.IsNumeric())
				continue;

			int32 counter = FCString::Atoi(*end);

			if (counter > maxCounter)
				maxCounter = counter;
		}

		return baseName + FString::FromInt(maxCounter + 1);
#else
		return FString(TEXT("EditorOnly"));
#endif

	}

	static URuntimeMeshComponent* duplicateRuntimeMeshComponentToActor(URuntimeMeshComponent * rmc, AActor * actor)
	{
		if (rmc->GetNumSections() == 0)
			return nullptr;

		URuntimeMeshComponent * newRMC = NewObject<URuntimeMeshComponent>(actor);

		newRMC->SetShouldSerializeMeshData(true);

		int32 newRMCSection = 0;

		auto sectionIds = rmc->GetRuntimeMesh()->GetSectionIds();

		auto materials = rmc->GetMaterials();

		auto defaultMaterial = rmc->GetMaterial(0);

		// Duplicate the sections and assign materials
		for (auto& sectionID : sectionIds)
		{
			if (!rmc->DoesSectionExist(sectionID))
				continue;

			auto& material = materials.IsValidIndex(sectionID) ? materials[sectionID] : defaultMaterial;

			auto readonlySection = rmc->GetSectionReadonly(sectionID);

			if (!readonlySection || readonlySection->NumIndices() == 0)
				continue;

			auto builder = MakeRuntimeMeshBuilder(*readonlySection);

			readonlySection->CopyTo(builder, true);

			newRMC->CreateMeshSection(newRMCSection, builder, false, EUpdateFrequency::Infrequent);
			newRMC->SetSectionMaterial(newRMCSection, material);
			newRMC->SetMaterial(newRMCSection, material);

			newRMCSection++;
		}

		newRMC->AttachToComponent(actor->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		actor->AddInstanceComponent(newRMC);

		newRMC->RegisterComponent();

		return newRMC;
	}

	static TArray<FTimStructBox> boxComponents(FGraphNode& node, FGraph& graph);

	static void unboxComponents(TArray<FTimStructBox>& components, FGraphNode& destination, FGraph& graph);


};

struct FQuatAverager
{
protected:
	FVector4 cumulative = FVector4(0.0f, 0.0f, 0.0f, 0.0f);

	// We pick the first quaternion as the hemisphere, other quaternions outside the hemisphere, are flipped back inside.
	FQuat firstRotation = FQuat::Identity; 
	int32 numAdded = 0;

public:
	void average(FQuat newRotation, float weight = 1.0f)
	{
		if (numAdded == 0)
		{
			firstRotation = newRotation;
		}

		newRotation.X *= weight;
		newRotation.Y *= weight;
		newRotation.Z *= weight;
		newRotation.W *= weight;

		numAdded++;

		//Before we add the new rotation to the average (mean), we have to check whether the quaternion has to be inverted. Because
		//q and -q are the same rotation, but cannot be averaged, we have to make sure they are all the same.
		const float DotResult = (newRotation | firstRotation);

		if ( DotResult < 0.0f ) {

			newRotation.X *= -1.0f;
			newRotation.Y *= -1.0f;
			newRotation.Z *= -1.0f;
			newRotation.W *= -1.0f;
		}

		cumulative.W += newRotation.W;
		cumulative.X += newRotation.X;
		cumulative.Y += newRotation.Y;
		cumulative.Z += newRotation.Z;
	}

	FQuat currentUnnormalized() 
	{ 
		const float addDet = 1.0f / float(numAdded);

		return FQuat(cumulative.X * addDet, cumulative.Y * addDet, cumulative.Z * addDet, cumulative.W * addDet);
	}

};

template<typename SubType, typename SourceType>
TArray<SubType*> SubArrayWithClass(TArray<SourceType*> array)
{
    TArray<SubType*> subarray;
    
    for( auto a : array )
    {
        auto b = Cast<SubType>(a);
        if( b )
            subarray.Add(b);
    }
    
    return subarray;
}

static void DrawCircle(ULineBatchComponent * lineBatch, const FVector& origin, const FVector& normal, const FLinearColor& color, const float radius, const float lineWidth, const int numSides = 32)
{
    FVector n = normal;
    n.Normalize();
    
    const FVector x = FVector::CrossProduct(FVector(1.0f, 0.0f, 0.0f), n);
    const FVector y = FVector::CrossProduct(FVector(0.0f, 1.0f, 0.0f), n);
    
    const float	angleDelta = 2.0f * PI / numSides;
    FVector	lastVertex = origin + x * radius;
    
    TArray<FBatchedLine> lines;
    lines.Reserve(numSides);
    
    for( int32 sideIndex = 0; sideIndex < numSides; sideIndex++ )
    {
        const FVector vertex = origin + (x * FMath::Cos(angleDelta * (sideIndex + 1)) + y * FMath::Sin(angleDelta * (sideIndex + 1))) * radius;
        
        lines.Emplace(lastVertex, vertex, color, -1.0f, lineWidth, 255);

        lastVertex = vertex;
    }
    
    lineBatch->DrawLines(lines);
}

template <typename T>
static FORCEINLINE T* LoadObjectFromPath(const FName& Path)
{
    if( Path == NAME_None )
        return nullptr;
    
    return Cast<T>(StaticLoadObject( T::StaticClass(), NULL,*Path.ToString()));
}

// -----------------------------------------------------------------------------
/// Vector conversions
// -----------------------------------------------------------------------------

inline FVector unreal(const Eigen::Vector3f& rhs)
{
    return FVector(rhs.x(), rhs.y(), rhs.z());
}

inline FVector unreal(const DDG::Vector v)
{
    return FVector(v.x, v.y, v.z);
}

inline Eigen::Vector3f eigen(const FVector& rhs)
{
    return Eigen::Vector3f(rhs.X, rhs.Y, rhs.Z);
}

inline Eigen::Vector3f eigen(const DDG::Vector v)
{
    return Eigen::Vector3f(v.x, v.y, v.z);
}

inline DDG::Vector to_tcods(FVector v)
{
    return DDG::Vector(v.X, v.Y, v.Z);
}

inline DDG::Vector to_tcods( Eigen::Vector3f v )
{
    return DDG::Vector( v.x(), v.y(), v.z() );
}


// -----------------------------------------------------------------------------
/// - Quaternion Conversions
// -----------------------------------------------------------------------------

inline Eigen::Quaternionf eigen(const FQuat& rhs)
{
    return Eigen::Quaternionf(rhs.W, rhs.X, rhs.Y, rhs.Z);
}

inline FQuat unreal(const Eigen::Quaternionf&  q)
{
    return FQuat(q.x(), q.y(), q.z(), q.w());
}

// -----------------------------------------------------------------------------
/// - Hyperplane Conversions
// -----------------------------------------------------------------------------

inline Eigen::Hyperplane<float,3> eigen(const FPlane& rhs)
{
    Eigen::Hyperplane<float,3> plane;
    
    
//    plane.coeffs()(0) = rhs[0];
//    plane.coeffs()(1) = rhs[1];
//    plane.coeffs()(2) = rhs[2];
//    plane.coeffs()(3) = rhs[3];
    
    // this is kinda messy, but I didn't see a nice way to do it (coeffs() is hard coded to check for <3 indices)
    plane = Eigen::Hyperplane<float,3>(Eigen::Vector3f(rhs.X, rhs.Y, rhs.Z), rhs.W);
    
    return plane;
}

// -----------------------------------------------------------------------------
/// - Triangles
// -----------------------------------------------------------------------------

// from http://blackpawn.com/texts/pointinpoly/
static bool sameSide( FVector& p1, FVector& p2, FVector& a, FVector& b )
{
	FVector c1 = FVector::CrossProduct( b - a, p1 - a );
	FVector c2 = FVector::CrossProduct( b - a, p2 - a );

	return FVector::DotProduct( c1, c2 ) >= 0;
}

static bool containsPoint( FVector& p, FVector& a, FVector& b, FVector& c )
{
	return sameSide( p, a, b, c ) && sameSide( p, b, a, c ) && sameSide( p, c, a, b );
}

/// <summary>
/// Calculates the intersection line segment between 2 lines (not segments).
/// Returns false if no solution can be found.
/// </summary>
/// <returns></returns>
static bool lineLineIntersection(
	const FVector& u0, const FVector& u1,
	const FVector& v0, const  FVector& v1, 
	FVector& uIntersection, FVector& vIntersection )
{
	// Algorithm is ported from the C# version (by Ronald Holthuizen) of 
	// Paul Bourke's C algorithm at http://paulbourke.net/geometry/pointlineplane/calclineline.cs


	FVector p1 = u0;
	FVector p2 = u1;
	FVector p3 = v0;
	FVector p4 = v1;

	FVector p13 = p1 - p3;
	FVector p43 = p4 - p3;

	if(p43.SizeSquared() < std::numeric_limits<float>::epsilon())
		return false;

	FVector p21 = p2 - p1;
	if(p21.SizeSquared() < std::numeric_limits<float>::epsilon())
		return false;

	double d1343 = p13.X * (double)p43.X + (double)p13.Y * p43.Y + (double)p13.Z * p43.Z;
	double d4321 = p43.X * (double)p21.X + (double)p43.Y * p21.Y + (double)p43.Z * p21.Z;
	double d1321 = p13.X * (double)p21.X + (double)p13.Y * p21.Y + (double)p13.Z * p21.Z;
	double d4343 = p43.X * (double)p43.X + (double)p43.Y * p43.Y + (double)p43.Z * p43.Z;
	double d2121 = p21.X * (double)p21.X + (double)p21.Y * p21.Y + (double)p21.Z * p21.Z;

	double denom = d2121 * d4343 - d4321 * d4321;
	if( FMath::Abs( denom ) < std::numeric_limits<double>::epsilon() ) 
		return false;

	double numer = d1343 * d4321 - d1321 * d4343;

	double mua = numer / denom;
	double mub = (d1343 + d4321 * (mua)) / d4343;

	uIntersection.X = (float)(p1.X + mua * p21.X);
	uIntersection.Y = (float)(p1.Y + mua * p21.Y);
	uIntersection.Z = (float)(p1.Z + mua * p21.Z);

	vIntersection.X = (float)(p3.X + mub * p43.X);
	vIntersection.Y = (float)(p3.Y + mub * p43.Y);
	vIntersection.Z = (float)(p3.Z + mub * p43.Z);

	return true;
}

 enum class IntersectionResult
 {
 	Parallel,
 	Coincident,
 	Intersecting,
 	NotIntersecting	
 };



static 
IntersectionResult intersection2D(
 	const FVector& u0, const FVector& u1,
 	const FVector& v0, const FVector& v1,
 	FVector& intersection )
 {
 	float denom = ((v1.Y - v0.Y)*(u1.X - u0.X)) -
 		((v1.X - v0.X)*(u1.Y - u0.Y));

 	float nume_a = ((v1.X - v0.X)*(u0.Y - v0.Y)) -
 		((v1.Y - v0.Y)*(u0.X - v0.X));

 	float nume_b = ((u1.X - u0.X)*(u0.Y - v0.Y)) -
 		((u1.Y - u0.Y)*(u0.X - v0.X));

 	if(denom == 0.0f)
 	{
 		if(nume_a == 0.0f && nume_b == 0.0f)
 		{
 			return IntersectionResult::Coincident;
 		}
 		return IntersectionResult::Parallel;
 	}

 	float ua = nume_a / denom;
 	float ub = nume_b / denom;

 	if(ua >= 0.0f && ua <= 1.0f && ub >= 0.0f && ub <= 1.0f)
 	{
 		// Get the intersection point.
 		intersection.X = u0.X + ua*(u1.X - u0.X);
 		intersection.Y = u0.Y + ua*(u1.Y - u0.Y);
        intersection.Z = u0.Z + ua*(u1.Z - u0.Z);

 		return IntersectionResult::Intersecting;
 	}

 	return IntersectionResult::NotIntersecting;
 }

static
IntersectionResult intersection2D(
	const FVector& u0, const FVector& u1,
	const FVector& v0, const FVector& v1,
	float &u // the intersection parameter of line segment (u0,u1)
)
{
	float denom = ((v1.Y - v0.Y)*(u1.X - u0.X)) -
		((v1.X - v0.X)*(u1.Y - u0.Y));

	float nume_a = ((v1.X - v0.X)*(u0.Y - v0.Y)) -
		((v1.Y - v0.Y)*(u0.X - v0.X));

	float nume_b = ((u1.X - u0.X)*(u0.Y - v0.Y)) -
		((u1.Y - u0.Y)*(u0.X - v0.X));

	if(denom == 0.0f)
	{
		if(nume_a == 0.0f && nume_b == 0.0f)
		{
			return IntersectionResult::Coincident;
		}
		return IntersectionResult::Parallel;
	}

	u = nume_a / denom;
	float v = nume_b / denom;

	if(u >= 0.0f && u <= 1.0f && v >= 0.0f && v <= 1.0f)
		return IntersectionResult::Intersecting;

	return IntersectionResult::NotIntersecting;
}

// -----------------------------------------------------------------------------
/// - Point Cloud Stuff
// -----------------------------------------------------------------------------
struct LeastSquaresRigidTransform
{
	Eigen::Vector3f translation; // translation component
	Eigen::Quaternionf rotation; // rotation part
};
// Adapted from
// Least-Square Rigid Motion Using SVD - Olga Sorkine
static LeastSquaresRigidTransform leastSquaresRigidTransform( Eigen::Matrix3Xf& P, Eigen::Matrix3Xf& Q, Eigen::VectorXf& w )
{
	using namespace Eigen;

	if(P.size() != Q.size())
		return{ Vector3f::Zero(), Quaternionf::Identity() };


	float wSum = 0.0f;
	for(int i = 0; i < w.size(); ++i)
		wSum += w( i );

	Eigen::Vector3f xCentroid( 0.0f, 0.0f, 0.0f );
	for(int i = 0; i < P.cols(); ++i)
		xCentroid += w( i ) * P.col( i );
	xCentroid /= wSum;

	Eigen::Vector3f yCentroid( 0.0f, 0.0f, 0.0f );
	for(int i = 0; i < P.cols(); ++i)
		yCentroid += w( i ) * Q.col( i );
	yCentroid /= wSum;

	Eigen::Matrix3Xf X = P;
	Eigen::Matrix3Xf Y = Q;

	for(int i = 0; i < X.cols(); ++i)
	{
		X.col( i ) -= xCentroid;
	}

	for(int i = 0; i < Y.cols(); ++i)
	{
		Y.col( i ) -= yCentroid;
	}

	auto wDiagonal = w.asDiagonal();

	MatrixXf S = (X * wDiagonal) * Y.transpose();

	JacobiSVD<MatrixXf> svd( S, ComputeFullU | ComputeFullV );

	MatrixXf U = svd.matrixU();
	MatrixXf V = svd.matrixV();

	// do we contain reflections? we can avoid them with the determinant
	// /see Least-Square Rigid Motion Using SVD - Olga Sorkine
	float det = (V * U.transpose()).determinant();

	Vector3f Det = Vector3f( 1.0f, 1.0f, det );

	auto R = (V * Det.asDiagonal()) * U.transpose();

	Vector3f t = yCentroid - R * xCentroid;

	return{ t, Quaternionf( R.block<3,3>( 0,0 ) ) };
}

// -----------------------------------------------------------------------------
/// - Math
// -----------------------------------------------------------------------------

struct Math
{
	// pi double
	static inline double pid() { return std::acos( -1.0 ); }

	// pi float
	static inline float pif() { return std::acos( -1.0f ); }
};

// -----------------------------------------------------------------------------
/// - Hashing
// -----------------------------------------------------------------------------

namespace rg
{
	// From: http://stackoverflow.com/questions/19195183/how-to-properly-hash-the-custom-struct
	template<typename T>
	void hash_combine (std::size_t & i, const T & v)
	{
		std::hash<T> h;
        i ^= h(v) + 0x9e3779b9 + (i << 6) + (i >> 2);
	}

}

// -----------------------------------------------------------------------------
/// - Scene View
// -----------------------------------------------------------------------------

// Utility
struct SceneViewAndFamily
{
public:
    FSceneView * sceneView = nullptr;
    FSceneViewFamilyContext * context = nullptr;
    
    FMatrix viewProjectionMatrix;
    FMatrix invViewProjectionMatrix;
    FMatrix invViewMatrix;
    FIntRect viewRect;
    FConvexVolume viewFrustum;
    
    bool didInit = false;
    
    ~SceneViewAndFamily()
    {

    }
    
    SceneViewAndFamily(UWorld * world)
    {
        init(world);
    }
    
    void init(UWorld * world)
    {
        if( context )
            delete context, context = nullptr;
        
        bool editor = false;
#if  WITH_EDITOR
        editor = world->WorldType == EWorldType::Editor;
        if( editor )
        {
            
            
            FViewport * viewport = GEditor->GetActiveViewport();
            FEditorViewportClient * viewportClient = (FEditorViewportClient*)viewport->GetClient();
            
            this->context = new FSceneViewFamilyContext(FSceneViewFamily::ConstructionValues(viewport,
                                                                                             viewportClient->GetScene(),
                                                                                             viewportClient->EngineShowFlags)
                                                        .SetRealtimeUpdate(true));
            
            this->sceneView = viewportClient->CalcSceneView(this->context);
        }
#endif
        
        APlayerController * playerController = world->GetFirstPlayerController();
        ULocalPlayer * localPlayer = playerController ? playerController->GetLocalPlayer() : nullptr;
        
        if( !editor && localPlayer && localPlayer->ViewportClient && localPlayer->ViewportClient->Viewport )
        {
            FViewport * viewport = localPlayer->ViewportClient->Viewport;
            
            this->context = new  FSceneViewFamilyContext(FSceneViewFamily::ConstructionValues(localPlayer->ViewportClient->Viewport,
                                                                                              world->Scene,
                                                                                              localPlayer->ViewportClient->EngineShowFlags)
                                                         .SetRealtimeUpdate(true) );
            
            FVector viewLocation;
            FRotator viewRotation;
            this->sceneView = localPlayer->CalcSceneView( this->context, /*out*/ viewLocation, /*out*/ viewRotation, localPlayer->ViewportClient->Viewport );
        }
        
        // be careful, the sceneView may not be active if we exit the render loop and do work on another thread
        // so, we can use this cached matrix to do our projections
        viewProjectionMatrix = sceneView->ViewMatrices.GetViewProjectionMatrix();
        invViewProjectionMatrix = sceneView->ViewMatrices.GetInvProjectionMatrix();
        invViewMatrix = sceneView->ViewMatrices.GetInvViewMatrix();
        
        
        
        viewRect = sceneView->UnscaledViewRect;
        viewFrustum = sceneView->ViewFrustum;
        
        didInit = true;
    }
    
    SceneViewAndFamily()
    {
        sceneView = nullptr;
        context = nullptr;
    }
    
    SceneViewAndFamily(SceneViewAndFamily&& a)
    {
        sceneView = a.sceneView;
        context = a.context;
        
        a.sceneView = nullptr;
        a.context = nullptr;
    }
    
    // worldPoint to screen point in the range [0,1]
    FVector2D worldToNormalizedPoint(const FVector& worldPoint) const
    {
        // this function is adapted from Unreal's
        // SceneView.cpp::ScreenToPixel(const FVector4& ScreenPoint,FVector2D& OutPixelLocation)
        
		// world to screen
		// screen to pixel

        FVector4 screenPoint = viewProjectionMatrix.TransformFVector4(FVector4(worldPoint,1));
        
        float invW = 1.0f / screenPoint.W;

        FVector2D point = FVector2D( 0.5f + screenPoint.X * 0.5f * invW,
									 0.5f - screenPoint.Y * 0.5f * invW);
//
//			OutPixelLocation = FVector2D(
//			UnscaledViewRect.Min.X + (0.5f + ScreenPoint.X * 0.5f * InvW) * UnscaledViewRect.Width(),
//			UnscaledViewRect.Min.Y + (0.5f - Y * 0.5f * InvW) * UnscaledViewRect.Height()
//			);
//        
////        FVector2D screenPoint = FVector2D(viewRect.Min.X + (0.5f + viewPoint.X * 0.5f * invW) * viewRect.Width(),
////                                         viewRect.Min.Y + (0.5f - y * 0.5f * invW) * viewRect.Height());
        
        return point;
    }
public:
    
    static void viewLocationRotation(UWorld * world, FVector& location, FRotator& rotation)
    {
        bool editor = false;

#if WITH_EDITOR
        editor = world->WorldType == EWorldType::Editor;
        if( editor )
        {
            location = GCurrentLevelEditingViewportClient->GetViewLocation();
            rotation = GCurrentLevelEditingViewportClient->GetViewRotation();
            
            return;
        }
#endif
        
        APlayerController * playerController = world->GetFirstPlayerController();
        
        if( !editor && playerController )
        {
            playerController->GetActorEyesViewPoint(location, rotation);
        }
    }
};

