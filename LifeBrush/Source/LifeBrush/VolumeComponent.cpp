// Copyright 2017 Code Monkey Castle, all rights reserved.

#include "LifeBrush.h"
#include "VolumeComponent.h"

#include "ChunkedMarchingCubes.h"

#include "Algorithm/nanoflann.hpp"

static float grad[12][3] = {
	{ 1.0,1.0,0.0 },{ -1.0,1.0,0.0 },{ 1.0,-1.0,0.0 },{ -1.0,-1.0,0.0 },
	{ 1.0,0.0,1.0 },{ -1.0,0.0,1.0 },{ 1.0,0.0,-1.0 },{ -1.0,0.0,-1.0 },
	{ 0.0,1.0,1.0 },{ 0.0,-1.0,1.0 },{ 0.0,1.0,-1.0 },{ 0.0,-1.0,-1.0 }
};

int perm[512] = { 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180, 151, 160, 137, 91, 90, 15, 131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23, 190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33, 88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166, 77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244, 102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196, 135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123, 5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42, 223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9, 129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228, 251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107, 49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254, 138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180 };

static float dot( float x, float y, float z, float* g ) {
	return x*g[0] + y*g[1] + z*g[2];
}

static float noise( float xin, float yin, float zin ) {
	float F3, G3, t, X0, Y0, Z0, x0, y0, z0, s, x1, y1, z1, x2, y2, z2, x3, y3, z3, t0, t1, t2, t3, n0, n1, n2, n3;
	int i, j, k, ii, jj, kk, i1, j1, k1, i2, j2, k2, gi0, gi1, gi2, gi3;

	F3 = 1.0 / 3.0;
	s = (xin + yin + zin)*F3;
	i = xin + s;
	j = yin + s;
	k = zin + s;
	G3 = 1.0 / 6.0;
	t = (i + j + k)*G3;
	X0 = i - t;
	Y0 = j - t;
	Z0 = k - t;
	x0 = xin - X0;
	y0 = yin - Y0;
	z0 = zin - Z0;

	if(x0 >= y0) {
		if(y0 >= z0) {
			i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
		}
		else if(x0 >= z0) {
			i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1;
		}
		else {
			i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1;
		}
	}
	else {
		if(y0 < z0) {
			i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1;
		}
		else if(x0 < z0) {
			i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1;
		}
		else {
			i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
		}
	}

	x1 = x0 - i1 + G3;
	y1 = y0 - j1 + G3;
	z1 = z0 - k1 + G3;
	x2 = x0 - i2 + 2.0*G3;
	y2 = y0 - j2 + 2.0*G3;
	z2 = z0 - k2 + 2.0*G3;
	x3 = x0 - 1.0 + 3.0*G3;
	y3 = y0 - 1.0 + 3.0*G3;
	z3 = z0 - 1.0 + 3.0*G3;

	ii = i & 255;
	jj = j & 255;
	kk = k & 255;

	gi0 = perm[ii + perm[jj + perm[kk]]] % 12;
	gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] % 12;
	gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] % 12;
	gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] % 12;

	t0 = 0.6 - x0*x0 - y0*y0 - z0*z0;
	if(t0<0) {
		n0 = 0.0;
	}
	else {
		t0 *= t0;
		n0 = t0 * t0 * dot( x0, y0, z0, grad[gi0] );
	}

	t1 = 0.6 - x1*x1 - y1*y1 - z1*z1;
	if(t1<0) {
		n1 = 0.0;
	}
	else {
		t1 *= t1;
		n1 = t1 * t1 * dot( x1, y1, z1, grad[gi1] );
	}

	t2 = 0.6 - x2*x2 - y2*y2 - z2*z2;
	if(t2<0) {
		n2 = 0.0;
	}
	else {
		t2 *= t2;
		n2 = t2 * t2 * dot( x2, y2, z2, grad[gi2] );
	}

	t3 = 0.6 - x3*x3 - y3*y3 - z3*z3;
	if(t3<0) {
		n3 = 0.0;
	}
	else {
		t3 *= t3;
		n3 = t3 * t3 * dot( x3, y3, z3, grad[gi3] );
	}

	return 16.0*(n0 + n1 + n2 + n3) + 1.0;
}

static float simplex_noise( int octaves, float x, float y, float z ) {
	float value = 0.0;
	int i;
	for(i = 0; i<octaves; i++) {
		value += noise(
			x*pow( 2, i ),
			y*pow( 2, i ),
			z*pow( 2, i )
		);
	}
	return value;
}

template<typename ElementType>
static void floating_rock( int size, UniformGrid<ElementType>& grid )
{
	for(int x = 1; x < size - 1; ++x)
	{
		float xf = float( x ) / float( size );

		for(int y = 1; y < size - 1; ++y)
		{
			float yf = float( y ) / float( size );

			for(int z = 1; z < size - 1; ++z)
			{
				float zf = float( z ) / float( size );



				float center_falloff = 0.1 / (pow( (xf - 0.5)*1.5, 2 ) + pow( (yf - 0.5)*1.5, 2 ) + pow( (zf - 0.5)*1.5, 2 ));

				float test = simplex_noise( 1, xf * 5, yf * 5, zf * 5 );

				float caves = pow( pow( simplex_noise( 1, xf * 5, yf * 5, zf * 5 ), -3.0 ), 3.0f );
				float density = simplex_noise( 5, xf, yf, zf ) * center_falloff;// *plateau_falloff;

				density -= caves;

				/*			density *= pow(
				noise( (xf + 1)*3.0, (yf + 1)*3.0, (zf + 1)*3.0 ) + 0.4, 1.8
				);*/



				//if(caves<0.5) 
				//	density = 0.0f;

				if(std::isinf( density ) || std::isnan( density ))
					density = 6.0f;



				grid( x, y, z ) = density;
			}
		}
	}
}

template<typename ElementType>
static void floating_rock_original( int size, UniformGrid<ElementType>& grid )
{
	for(int x = 1; x < size - 1; ++x)
	{
		float xf = float( x ) / float( size );

		for(int y = 1; y < size - 1; ++y)
		{
			float yf = float( y ) / float( size );

			for(int z = 1; z < size - 1; ++z)
			{
				float zf = float( z ) / float( size );

				float plateau_falloff = 0.0f;

				if(yf <= 0.8f)
					plateau_falloff = 1.0f;
				else if(0.8f < yf && yf < 0.9f)
					plateau_falloff = 1.0f - (yf - 0.8f) * 10.0f;
				else
					plateau_falloff = 0.0f;


				float center_falloff = 0.1 / (pow( (xf - 0.5)*1.5, 2 ) + pow( (yf - 1.0)*0.8, 2 ) + pow( (zf - 0.5)*1.5, 2 ));

				float caves = pow( simplex_noise( 1, xf * 5, yf * 5, zf * 5 ), 3 );
				float density = (
					simplex_noise( 5, xf, yf*0.5, zf ) *
					center_falloff *
					plateau_falloff
					);

				density *= pow(
					noise( (xf + 1)*3.0, (yf + 1)*3.0, (zf + 1)*3.0 ) + 0.4, 1.8
				);

				if(caves < 0.5)
					density = 0.0f;

				if(std::isinf( density ) || std::isnan( density ))
					density = 100.0f;

				grid( x, y, z ) = density;
			}
		}
	}
}



struct WormSegment
{
	FVector position;

	FVector last;
	FVector next;
};

struct WormCloud
{
	std::vector<WormSegment>  pts;

	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	inline float kdtree_distance( const float *p1, const size_t idx_p2, size_t /*size*/ ) const
	{
		const auto& p2_3f = pts[idx_p2].position;
		const auto& p1_3f = *reinterpret_cast<const FVector*>(p1);

		return (p2_3f - p1_3f).SizeSquared();
	}

	inline float kdtree_get_pt( const size_t idx, int dim ) const
	{
		return pts[idx].position[dim];
	}

	template <class BBOX>
	bool kdtree_get_bbox( BBOX& /*bb*/ ) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<float, WormCloud>,
	WormCloud,
	3 /* dim */
> WormCloudAdaptor;

float wormDistance( WormSegment& segment, FVector point )
{
	float d1 = FMath::PointDistToSegment( point, segment.last, segment.position );
	float d2 = FMath::PointDistToSegment( point, segment.position, segment.next );

	return FMath::Min( d1, d2 );
}

std::vector<FVector> UVolumeComponent::perlin_worm(FIntVector start, const unsigned int nSegments)
{
	FVector p(start);

	std::vector<FVector> segments;
	segments.reserve(nSegments);

	float segmentLength = 1.0f;

	for (unsigned int i = 0; i < nSegments; i++)
	{
		segments.push_back(p);

		FVector offset;

		FVector noisePosition = p * wormNoiseScale;

		float u = simplex_noise(3, noisePosition.X, noisePosition.Y, noisePosition.Z);
		float v = simplex_noise(3, noisePosition.X + 200.0f, noisePosition.Y + 300.0f, noisePosition.Z);

		offset.X = sin(u) * cos(v);
		offset.Y = cos(u) * cos(u);
		offset.Z = sin(v);

		p += offset;
	}

	return segments;
}

void UVolumeComponent::decrepify( UniformGrid<float>& grid, size_t wormSegments, uint32 numWorms_in, float wormRadius, FRandomStream& randomStream )
{
	const FIntVector dimensions = grid.dimensions();

	WormCloud wormCloud;


	for(uint32 wormI = 0; wormI < numWorms_in; wormI++)
	{
		FIntVector start(
			randomStream.RandRange( 0, dimensions.X ),
			randomStream.RandRange( 0, dimensions.Y ),
			randomStream.RandRange( 0, dimensions.Z )
		);

		std::vector<FVector> worm = perlin_worm( start, wormSegments );

		size_t worm_n = worm.size();
		for(int i = 0; i < worm_n; ++i)
		{
			FVector p = worm[i];

			FVector last = i == 0 ? p : worm[i - 1];
			FVector next = (i + 1) < worm_n ? worm[i + 1] : p;

			wormCloud.pts.emplace_back( WormSegment{ p, last, next } );
		}
	}

	WormCloudAdaptor wormIndex( 3, wormCloud, nanoflann::KDTreeSingleIndexAdaptorParams( 2 ) );
	wormIndex.buildIndex();

	if(wormCloud.pts.size() == 0)
		return;

	for(int x = 1; x < dimensions.X - 1; ++x)
	{
		float xf = float( x ) / float( dimensions.X );

		for(int y = 1; y < dimensions.Y - 1; ++y)
		{
			float yf = float( y ) / float( dimensions.Y );

			for(int z = 1; z < dimensions.Z - 1; ++z)
			{
				float zf = float( z ) / float( dimensions.Z );

				float wormD = 0.0f;

				FVector p( x, y, z );

				size_t index;
				float distanceSqrd;

				wormIndex.knnSearch( &p.X, 1, &index, &distanceSqrd );

				WormSegment& segment = wormCloud.pts[index];

				const float radiusSqrd = wormRadius * wormRadius;

				//if(distanceSqrd < radiusSqrd)
				//	grid( x, y, z ) = 0.0f;


				//float distance = FMath::Sqrt( distanceSqrd );

				float distance = wormDistance( segment, p );

				if(distance > wormRadius)
					distance = wormRadius;

				float density = (FMath::Cos( 10.0f * distance / (PI * wormRadius) ) + 1) * (6.0f / 2.0f);



				float& cell = grid( x, y, z );

				cell -= density;

				if(cell < 0.0f)
					cell = 0.0f;

				//float sub = 100.0f * (radius * radius) / distance;


				//float& g = grid( x, y, z );

				//g -= 100.0f * (radius * radius) / distance;

				//if(distance < 0.001f)
				//	g = 0.0f;

				//if(g < 0.0f)
				//	g = 0.0f;

				//grid( x, y, z ) -= distance;
			}
		}
	}

}

template<typename ElementType>
static void emptyVolume(UniformGrid<ElementType>& grid)
{
	FIntVector dim = grid.dimensions();

	for (int z = 0; z < dim.Z; z++)
	{
		for (int y = 0; y < dim.Y; y++)
		{
			for (int x = 0; x < dim.X; x++)
			{
				grid(x, y, z) = 0.0f;
			}
		}
	}
}

template<typename ElementType>
static void sphere( int size, UniformGrid<ElementType>& grid )
{
	FVector centre( 0.5, 0.5, 0.5 );



	for(int x = 1; x < size - 1; ++x)
	{
		float xf = float( x ) / float( size );

		for(int y = 1; y < size - 1; ++y)
		{
			float yf = float( y ) / float( size );

			for(int z = 1; z < size - 1; ++z)
			{
				float zf = float( z ) / float( size );

				FVector p( xf, yf, zf );

				float caves = pow( simplex_noise( 1, xf * 5, yf * 5, zf * 5 ), 3 );

				float falloff = 0.4f - FVector::Dist( p, centre );

				{
					float crust = 0.05f;

					// the falloff function looks like this
					//        ___
					//       /
					//______/
					//
					//      | |
					//     crust thickness
					if(falloff > 0.0f)
						falloff = 1.0f;
					else if(falloff <= 0.0f && falloff > -crust)
						falloff = (falloff + crust) / crust; // linear falloff for the crust thickness
					else
						falloff = 0.0f;
				}


				float density = simplex_noise( 5, xf, yf, zf ) * falloff;

				density += simplex_noise( 3, xf*6.0, yf*6.0, zf*6.0 ) * FMath::Lerp( 0.0f, 3.0f, falloff );

				density -= caves;
				if(caves < 0.5)
					density = 0.0f;

				if(std::isinf( density ) || std::isnan( density ))
					density = 100.0f;

				grid( x, y, z ) = density;
			}
		}
	}
}

// Sets default values for this component's properties
UVolumeComponent::UVolumeComponent()
{
	PrimaryComponentTick.bCanEverTick = false;
}


// Called when the game starts
void UVolumeComponent::BeginPlay()
{
	Super::BeginPlay();

	build();
}

URuntimeMeshComponent* UVolumeComponent::runtimeMesh()
{
	URuntimeMeshComponent * mesh = GetOwner()->FindComponentByClass<URuntimeMeshComponent>();

	return mesh;
}
void UVolumeComponent::build()
{
	FIntVector dimensions( size );
	FVector cellSize( 1.0f );
	FVector extents = Utility::scale( dimensions, cellSize );
	FVector worldMin = extents * -0.5f;

	std::vector<float> samples;

	_grid.init( dimensions, cellSize, worldMin, samples );

	int dSize = dimensions.X - 1;

	switch(method)
	{
	case EBlocksGenerationMethod::FloatingIsland:
		floating_rock( dSize, _grid );
		break;

	case EBlocksGenerationMethod::OriginalFloatingIsland:
		floating_rock_original( dSize, _grid );
		break;

	case EBlocksGenerationMethod::Sphere:
		sphere( dSize, _grid );
		break;

	case EBlocksGenerationMethod::Empty:
		emptyVolume(_grid);

	default:
		break;
	}

	FRandomStream randomStream( randomSeed );
	decrepify( _grid, size, numWorms, wormCellRadius, randomStream );

	FIntVector min = FIntVector::ZeroValue;
	FIntVector max = dimensions;
	auto meshFactory = MarchingCubes::marchingCubes( isoLevel, _grid, min, max, uvScale, _cachedGridIndices, true, flatShading );

	URuntimeMeshComponent * mesh = runtimeMesh();

	if(!mesh)
	{
		mesh = NewObject<URuntimeMeshComponent>( GetOwner(), TEXT( "runtime_voxel_mesh" ) );

		mesh->AttachToComponent( GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );
		mesh->SetMaterial(0, material);

		mesh->RegisterComponent();
	}

	mesh->CreateMeshSection( 0, meshFactory.vertices, meshFactory.indices, meshFactory.normals, meshFactory.uvs /* uvs*/, TArray<FColor>() /* colors */, meshFactory.tangents, false /* create collision */, EUpdateFrequency::Frequent );
	mesh->SetMaterial(0, material);
	mesh->SetMeshSectionCastsShadow( 0, true );
	mesh->SetMeshSectionCollisionEnabled( 0, true );
	mesh->bMultiBodyOverlap = true;
}

UniformGrid<float>& UVolumeComponent::grid()
{
	return _grid;
}

FIntVector UVolumeComponent::componentToIndex( FVector localPoint )
{
	return _grid.sampleIndex( localPoint );
}

FVector UVolumeComponent::indexToComponent( FIntVector index )
{
	return _grid.samplePoint( index );
}

static inline bool overlaps( const FIntVector& p, const FIntVector& min, const FIntVector& max )
{
	return p.X >= min.X && p.X <= max.X &&
		p.Y >= min.Y && p.Y <= max.Y &&
		p.Z >= min.Z && p.Z <= max.Z;
}

FVector _trim(FVector4& v)
{
	return FVector(v.X, v.Y, v.Z);
}

void UVolumeComponent::markDirty( FVector minExtents, FVector maxExtents )
{
	FVector cellSize = _grid.cellSize();

	FIntVector minIndex = _grid.sampleIndex( minExtents );
	FIntVector maxIndex = _grid.sampleIndex( maxExtents );

	if(_marching)
		return;

	_marching = true;

	// we use copy, so that we get a copy of the grid, so that it can still be mutated in the background
	Async<void>( EAsyncExecution::ThreadPool, [=]() mutable 
	{
		std::vector<int32> newIndices;

		auto meshFactory = MarchingCubes::marchingCubes( isoLevel, _grid, minIndex, maxIndex, uvScale, newIndices, true, flatShading );

		AsyncTask( ENamedThreads::GameThread, [=]() mutable
		{
			URuntimeMeshComponent * mesh = runtimeMesh();

			if(!mesh)
			{
				mesh = NewObject<URuntimeMeshComponent>( GetOwner(), TEXT( "runtime_voxel_mesh" ) );
				mesh->AttachToComponent( GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform );

				mesh->RegisterComponent();
			}

			if (mesh->GetNumSections())
			{
				auto section = mesh->GetSectionReadonly(0);

				int32 n = section->NumIndices();
				int32 i = 0;
				while (i < n)
				{
					int32 i0 = section->GetIndex(i++);
					int32 i1 = section->GetIndex(i++);
					int32 i2 = section->GetIndex(i++);

					FVector v0 = section->GetPosition(i0);
					FVector4 n0 = section->GetNormal(i0);
					FVector2D u0 = section->GetUV(i0, 0);

					FVector v1 = section->GetPosition(i1);
					FVector4 n1 = section->GetNormal(i1);
					FVector2D u1 = section->GetUV(i1, 0);

					FVector v2 = section->GetPosition(i2);
					FVector4 n2 = section->GetNormal(i2);
					FVector2D u2 = section->GetUV(i2, 0);

					FIntVector vi0 = _grid.inflate(_cachedGridIndices[i0 / 3]);

					// check if the grid index is contained in our regenerated area
					if (overlaps(vi0, minIndex, maxIndex))
						continue;

					meshFactory.pushTriangle(v0, v1, v2, _trim(n0), _trim(n1), _trim(n2), u0, u1, u2);

					newIndices.push_back(_cachedGridIndices[i0 / 3]);
				}
			}			

			mesh->CreateMeshSection( 0, meshFactory.vertices, meshFactory.indices, meshFactory.normals, meshFactory.uvs /* uvs*/, TArray<FColor>() /* colors */, meshFactory.tangents, false /* create collision */, EUpdateFrequency::Frequent );

			_cachedGridIndices = newIndices;

			_marching = false;
		} );
	} );
}




































std::vector<FVector> UChunkedVolumeComponent::perlin_worm(FIntVector start, const unsigned int nSegments)
{
	FVector p(start);

	std::vector<FVector> segments;
	segments.reserve(nSegments);

	float segmentLength = 1.0f;

	for (unsigned int i = 0; i < nSegments; i++)
	{
		segments.push_back(p);

		FVector offset;

		FVector noisePosition = p * wormNoiseScale;

		float u = simplex_noise(3, noisePosition.X, noisePosition.Y, noisePosition.Z);
		float v = simplex_noise(3, noisePosition.X + 200.0f, noisePosition.Y + 300.0f, noisePosition.Z);

		offset.X = sin(u) * cos(v);
		offset.Y = cos(u) * cos(u);
		offset.Z = sin(v);

		p += offset;
	}

	return segments;
}

void UChunkedVolumeComponent::decrepify(UniformGrid<float>& grid, size_t wormSegments, uint32 numWorms_in, float wormRadius, FRandomStream& randomStream)
{
	const FIntVector dimensions = grid.dimensions();

	WormCloud wormCloud;


	for (uint32 wormI = 0; wormI < numWorms_in; wormI++)
	{
		FIntVector start(
			randomStream.RandRange(0, dimensions.X),
			randomStream.RandRange(0, dimensions.Y),
			randomStream.RandRange(0, dimensions.Z)
		);

		std::vector<FVector> worm = perlin_worm(start, wormSegments);

		size_t worm_n = worm.size();
		for (int i = 0; i < worm_n; ++i)
		{
			FVector p = worm[i];

			FVector last = i == 0 ? p : worm[i - 1];
			FVector next = (i + 1) < worm_n ? worm[i + 1] : p;

			wormCloud.pts.emplace_back(WormSegment{ p, last, next });
		}
	}

	WormCloudAdaptor wormIndex(3, wormCloud, nanoflann::KDTreeSingleIndexAdaptorParams(2));
	wormIndex.buildIndex();

	if (wormCloud.pts.size() == 0)
		return;

	for (int x = 1; x < dimensions.X - 1; ++x)
	{
		float xf = float(x) / float(dimensions.X);

		for (int y = 1; y < dimensions.Y - 1; ++y)
		{
			float yf = float(y) / float(dimensions.Y);

			for (int z = 1; z < dimensions.Z - 1; ++z)
			{
				float zf = float(z) / float(dimensions.Z);

				float wormD = 0.0f;

				FVector p(x, y, z);

				size_t index;
				float distanceSqrd;

				wormIndex.knnSearch(&p.X, 1, &index, &distanceSqrd);

				WormSegment& segment = wormCloud.pts[index];

				const float radiusSqrd = wormRadius * wormRadius;

				//if(distanceSqrd < radiusSqrd)
				//	grid( x, y, z ) = 0.0f;


				//float distance = FMath::Sqrt( distanceSqrd );

				float distance = wormDistance(segment, p);

				if (distance > wormRadius)
					distance = wormRadius;

				float density = (FMath::Cos(10.0f * distance / (PI * wormRadius)) + 1) * (6.0f / 2.0f);



				float& cell = grid(x, y, z);

				cell -= density;

				if (cell < 0.0f)
					cell = 0.0f;

				//float sub = 100.0f * (radius * radius) / distance;


				//float& g = grid( x, y, z );

				//g -= 100.0f * (radius * radius) / distance;

				//if(distance < 0.001f)
				//	g = 0.0f;

				//if(g < 0.0f)
				//	g = 0.0f;

				//grid( x, y, z ) -= distance;
			}
		}
	}

}




// Sets default values for this component's properties
UChunkedVolumeComponent::UChunkedVolumeComponent()
{
	PrimaryComponentTick.bCanEverTick = true;
	bIsActive = true;
}


// Called when the game starts
void UChunkedVolumeComponent::BeginPlay()
{
	Super::BeginPlay();

	init();
}

void UChunkedVolumeComponent::TickComponent(float DeltaTime, enum ELevelTick TickType, FActorComponentTickFunction *ThisTickFunction)
{
	Super::TickComponent(DeltaTime, TickType, ThisTickFunction);

	if (_marching)
		return;

	_marching = true;

	auto gridAccess = std::move(_gridAccess);
	auto dirtyAcces = std::move(_dirtyAccess);

	_gridAccess = nullptr;
	_dirtyAccess = nullptr;
	
	Async<void>(EAsyncExecution::ThreadPool, [this, gridAccess, dirtyAcces]() mutable
	{
		if (gridAccess) gridAccess(_grid);
		if (dirtyAcces) dirtyAcces(_grid);

		_marching = false;
	});
}

URuntimeMeshComponent* UChunkedVolumeComponent::runtimeMesh()
{
	URuntimeMeshComponent * mesh = GetOwner()->FindComponentByClass<URuntimeMeshComponent>();

	if (!mesh)
	{
		mesh = NewObject<URuntimeMeshComponent>(GetOwner(), TEXT("runtime_voxel_mesh"));
		mesh->AttachToComponent(GetOwner()->GetRootComponent(), FAttachmentTransformRules::KeepRelativeTransform);

		mesh->RegisterComponent();
	}

	mesh->SetMeshSectionCastsShadow(0, true);
	mesh->SetMeshSectionCollisionEnabled(0, false);
	mesh->SetMaterial(0, material);
	mesh->bMultiBodyOverlap = true;

	return mesh;
}

void UChunkedVolumeComponent::init()
{
	FVector cellSize(1.0f);
	FIntVector dimensions(size);

	_cellSize = cellSize;
	for (int c = 0; c < 3; ++c)
		_invCellSize[c] = 1.0f / cellSize[c];

	_grid.init(dimensions, cellSize);
}

void UChunkedVolumeComponent::initializeTestVolume()
{
	init();

	const int radius = 5;
	const float maxValue = 10.0f;

	FIntVector start(-128, -radius, -radius);
	FIntVector end(0, radius, radius);

	FIntVector index;

	for (index.Z = start.Z; index.Z <= end.Z; index.Z++)
	{
		for (index.Y = start.Y; index.Y <= end.Y; index.Y++)
		{
			for (index.X = start.X; index.X <= end.X; index.X++)
			{
				float r = std::sqrt(std::pow(index.Z, 2.0f) + std::pow(index.Y, 2.0f));

				float d = (radius - r) / float(radius) * maxValue;

				_grid.set(index,d);
			}
		}
	}

	for (index.Z = start.Z; index.Z <= end.Z; index.Z++)
	{
		for (index.Y = start.Y; index.Y <= end.Y; index.Y++)
		{
			for (index.X = start.X; index.X <= end.X; index.X++)
			{
				float r = std::sqrt(std::pow(index.Z, 2.0f) + std::pow(index.Y, 2.0f));

				float d = (radius - r) / float(radius) * maxValue;

				_grid.set(index, d);
			}
		}
	}


	PaddedUniformGrid<float> chunk0 = _grid.chunkAtChunkIndex(FIntVector(0, 0, 0));
	PaddedUniformGrid<float> chunk1 = _grid.chunkAtChunkIndex(FIntVector(0, 0, -1));

	float a0 = chunk0(0, 0, -1);
	float a1 = chunk1(0, 0, 63);

	float b0 = chunk0(0, 0, 0);
	float b1 = chunk1(0, 0, 64);

	if (a0 != a1)
		UE_LOG(LogTemp, Warning, TEXT("a0 != a1"));

	if (b0 != b1)
		UE_LOG(LogTemp, Warning, TEXT("b0 != b1"));
}

void UChunkedVolumeComponent::markTestVolumeDirty()
{
	FIntVector min = _grid.gridIndexFromChunkIndex(_grid.minChunkndex());
	FIntVector max = _grid.gridIndexFromChunkIndex(_grid.maxChunkIndex()) + _grid.chunkDimensions() - FIntVector(1);

	markDirtyIndex(min, max);
}

FIntVector UChunkedVolumeComponent::componentToIndex(FVector localPoint)
{
	return FIntVector(localPoint * _invCellSize);
}

FVector UChunkedVolumeComponent::indexToComponent(FIntVector index)
{ 
	return Utility::scale(index, _cellSize);
}

void UChunkedVolumeComponent::markDirtyIndex(FIntVector minGridIndex, FIntVector maxGridIndex)
{
	_dirtyAccess = [=](ChunkGrid<float>& grid)  
	{
		auto chunkGeometries = ChunkedMarchingCubes::marchingCubes_byChunk(isoLevel, grid, minGridIndex, maxGridIndex, uvScale, true, flatShading);

		AsyncTask(ENamedThreads::GameThread, [=]() mutable
		{
			URuntimeMeshComponent * mesh = runtimeMesh();

			for (auto& pair : chunkGeometries)
			{
				FIntVector chunkIndex = pair.first;

				ChunkedMarchingCubes::ChunkGeometry& chunkGeometry = pair.second;

				SmoothMeshFactory& meshFactory = chunkGeometry.factory;
				std::vector<FIntVector> chunkIndices = chunkGeometry.cellIndices;


				auto chunkIter = _chunkInfo.find(chunkIndex);

				bool isNewChunk = chunkIter == _chunkInfo.end();

				ChunkInfo& chunkInfo = isNewChunk ? _chunkInfo[chunkIndex] : chunkIter->second;

				if (isNewChunk)
				{
					chunkInfo.sectionIndex = _chunkInfo.size() - 1;
				}
				// we are updating a chunk, gather the old, non-overlapping mesh information for the section
				else if (!isNewChunk && mesh->DoesSectionExist(chunkInfo.sectionIndex))
				{
					auto section = mesh->GetSectionReadonly(chunkInfo.sectionIndex);

					int32 n = section->NumIndices();
					int32 i = 0;
					while (i < n)
					{
						int32 i0 = section->GetIndex(i++);
						int32 i1 = section->GetIndex(i++);
						int32 i2 = section->GetIndex(i++);

						FVector v0 = section->GetPosition(i0);
						FVector4 n0 = section->GetNormal(i0);
						FVector2D u0 = section->GetUV(i0, 0);

						FVector v1 = section->GetPosition(i1);
						FVector4 n1 = section->GetNormal(i1);
						FVector2D u1 = section->GetUV(i1, 0);

						FVector v2 = section->GetPosition(i2);
						FVector4 n2 = section->GetNormal(i2);
						FVector2D u2 = section->GetUV(i2, 0);

						// check if the grid index is contained in our regenerated area
						FIntVector cachcedCellIndex = chunkInfo.chunkCellIndexForVertex[i0 / 3];

						FIntVector gridCachcedCellIndex = _grid.gridIndexFromChunkIndex(chunkIndex) + cachcedCellIndex;

						if (overlaps(gridCachcedCellIndex, minGridIndex, maxGridIndex))
							continue;

						meshFactory.pushTriangle(v0, v1, v2, _trim(n0), _trim(n1), _trim(n2), u0, u1, u2);

						chunkIndices.push_back(cachcedCellIndex);
					}
				}

				chunkInfo.chunkCellIndexForVertex = chunkIndices;

				mesh->CreateMeshSection(
					chunkInfo.sectionIndex,
					meshFactory.vertices,
					meshFactory.indices,
					meshFactory.normals,
					meshFactory.uvs /* uvs*/,
					TArray<FColor>() /* colors */,
					meshFactory.tangents,
					false /* create collision */,
					EUpdateFrequency::Frequent
				);

				mesh->SetMaterial(chunkInfo.sectionIndex, material);

				mesh->SetMeshSectionCollisionEnabled(chunkInfo.sectionIndex, false);
			}
		});
	};
}

void UChunkedVolumeComponent::markDirty(FVector minExtents, FVector maxExtents)
{
	FIntVector minGridIndex = _grid.sampleIndex(minExtents);
	FIntVector maxGridIndex = _grid.sampleIndex(maxExtents);

	markDirtyIndex(minGridIndex, maxGridIndex);
}
