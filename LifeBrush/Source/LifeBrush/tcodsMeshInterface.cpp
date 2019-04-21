//
//  tcodsMeshInterface.cpp
//  RegionGrowing
//
//  Created by Timothy Davison on 2015-08-21.
//  Copyright (c) 2015 Timothy Davison. All rights reserved.
//

#include "LifeBrush.h"
#include "tcodsMeshInterface.h"

#include "Algorithm/tcods/HalfEdge.h"
#include "Algorithm/tcods/MeshIO.h"
#include "Algorithm/SurfaceIndex.h"
#include "RuntimeMeshComponent.h"
#include "RuntimeMeshLibrary.h"
#include "MeshFactory.h"
#include "ChunkedMarchingCubes.h"

#include <vcg/complex/complex.h>

#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clustering.h>





#include "StaticMeshResources.h"
#include "Utility.h"


#include <vector>
#include <unordered_set>
#include <map>
#include <algorithm>




auto tcodsMeshInterfaceBase::indexOfVertex(const FVector& uStaticMeshVertex, int32& index_out) -> bool
{
	auto found = _vertexLookup.find( uStaticMeshVertex );

    if(found == _vertexLookup.end() )
        return false;
    
	index_out = found->second.vertexIndex;
    
    return true;
}

auto tcodsMeshInterfaceBase::indexOfVertex(const FVector& uStaticMeshVertex, VertexIndex& index_out) -> bool
{
	auto found = _vertexLookup.find(uStaticMeshVertex);

	if (found == _vertexLookup.end())
		return false;

	index_out = found->second;

	return true;
}

void tcodsMeshInterfaceBase::rebuildBvh()
{
	_initBVH();
}

auto tcodsMeshInterfaceBase::_countTriangles() -> size_t
{
	size_t count = 0;

	for( auto& pair : _sections )
	{
		count += pair.second->faces.size();
	}

	return count;
}

void tcodsMeshInterfaceBase::_initBVH()
{
    std::vector<bvh::Triangle> triangles;
    triangles.reserve(_countTriangles());

	for (auto& pair : _sections)
	{
		tcods::Mesh& mesh = *pair.second.get();

		FSurfaceIndex surfaceIndex;

		surfaceIndex.sectionIndex = pair.first;
		surfaceIndex.faceIndex = 0;

		for (const tcods::Face& face : mesh.faces)
		{
			auto a = face.he->from->position;
			auto b = face.he->next->from->position;
			auto c = face.he->next->next->from->position;

			triangles.emplace_back();
			bvh::Triangle& t = triangles.back();
			t.vertices = { { unreal(a), unreal(b), unreal(c) } };
			t.surfaceIndex = surfaceIndex;

			++surfaceIndex.faceIndex;
		}
	}
    
    bvh.build(triangles);
}




void tcodsMeshInterfaceBase::_initNearestKDTree(FBox limits)
{
	using namespace vcg;

	// count the number of triangles
	size_t n = 0;

	for (auto& pair : _sections)
		n += pair.second->faces.size();

	// build a list of triangle samples
	std::vector<TriangleSample> samples;
	samples.reserve(n);

	FBox aabb(EForceInit::ForceInitToZero);

	for (auto& pair : _sections)
	{
		auto sectionIndex = pair.first;
		auto& mesh = *pair.second;

		// grab all the triangles from this mesh
		FSurfaceIndex index;
		index.sectionIndex = sectionIndex;
		index.faceIndex = 0;

		for (auto& face : mesh.faces)
		{
			samples.emplace_back();
			TriangleSample& sample = samples.back();

			auto he = face.he;

			int32 ti = 0;
			do 
			{
				assert(ti < 3);

				FVector p = unreal(he->from->position);

				aabb += p;

				sample.triangle[ti] = p;

				ti++;
				he = he->next;
			} while (he != face.he);

			bool inLimits = sample.hackIntersectsBox(limits);

			inLimits = true;

			if (inLimits)
			{
				sample.surfaceIndex = index;

				index.faceIndex++;
			}
			else // discard
				samples.pop_back();
		}
	}

	_faceTree = std::make_unique<vcg::KdTreeFace>(samples, aabb);
}

auto tcodsMeshInterfaceBase::_clear()
{
	_sections.clear();

	_vertexLookup.clear();
	_halfEdgeVertex_to_runtimeMeshVertices.clear();
}


auto tcodsMeshInterfaceBase::nearestPointOnMesh(const FVector& point) -> SurfacePoint
{
	SurfacePoint result;

	if (!_faceTree)
		return result;

	vcg::KdTreeFace::ObjectMarker marker;
	float distance = 0.0f;
	FVector nearestPoint;

	auto triangleSample = _faceTree->doQueryClosest(point, nearestPoint, distance, marker);

	if( !triangleSample )
		return result;

	result.point = nearestPoint;
	result.surfaceIndex = triangleSample->surfaceIndex;

	// we should really check the point against the limits, and project the point onto the limits
	// we could do this with edge segments

	return result;
}


void tcodsMeshInterfaceBase::_floodVisit(
	const int32 vertex_in, 
	const std::vector< std::vector<int32> >& verticesConnectedToVertex, 
	std::set<int>& unvisitedVertices, 
	std::unordered_set<int32>& connectedVertices
)
{
	using namespace std;
	using namespace tcods;

	using Index = tcods::MeshIO::Index;

	std::vector<int32> vertexStack;

	vertexStack.push_back(vertex_in);

	do
	{
		int32 vertex = vertexStack.back();
		vertexStack.pop_back();

		auto foundUnvisited = unvisitedVertices.find(vertex);

		// we have visited this vertex, nothing to do
		if (foundUnvisited == unvisitedVertices.end())
			continue;

		// we are visiting this vertex, remove it from the set
		unvisitedVertices.erase(foundUnvisited);

		connectedVertices.insert(vertex);

		// visit the connected
		auto& connected = verticesConnectedToVertex[vertex];
		
		vertexStack.insert(vertexStack.end(), connected.begin(), connected.end());

	} while (vertexStack.size());
}

auto tcodsMeshInterfaceBase::_extractSections(tcods::MeshIO::MeshData& data) -> std::vector<tcods::MeshIO::MeshData>
{
	using namespace std;
	using namespace tcods;
	using namespace DDG;

	using Index = tcods::MeshIO::Index;

	vector<MeshIO::MeshData> result;

	// first, associate faces with each vertex
	vector< vector<int32> > verticesConnectedToVertex(data.positions.size());

	for (int fi = 0; fi < data.indices.size(); ++fi)
	{
		auto& face = data.indices[fi];

		for (Index& index : face)
		{
			auto& vertexFaces = verticesConnectedToVertex[index.position];

			for (Index& index : face)
			{
				vertexFaces.push_back(index.position);
			}
		}
	}

	// track unvisited faces
	set<int> unvisitedVertices;
	for (int i = 0; i < data.positions.size(); ++i) unvisitedVertices.insert(i);

	while (unvisitedVertices.size())
	{
		// grab the first unvisited vertex
		auto begin = unvisitedVertices.begin();

		int unvisitedI = *begin;
		auto& unvisitedFace = data.indices[unvisitedI];
		
		// idea, visit all the faces in unvisited and check if they reference something in connected
		unordered_set<int> connectedVertices;

		_floodVisit(unvisitedI, verticesConnectedToVertex, unvisitedVertices, connectedVertices);

		// convert them to a mesh-data
		result.emplace_back();

		MeshIO::MeshData& sectionData = result.back();

		auto& indices = sectionData.indices;
		auto& positions = sectionData.positions;
		auto& normals = sectionData.normals;
		auto& texcoords = sectionData.texcoords;

		std::unordered_map<int, int> vertexToIndex;

		for (auto& face : data.indices)
		{
			bool shouldAdd = [&]() {
				for (auto& index : face)
				{
					if (connectedVertices.find(index.position) != connectedVertices.end())
					{
						return true;
					}
				}

				return false;
			}();

			if( !shouldAdd )
				continue;

			indices.emplace_back();
			vector<Index>& newFace = indices.back();
			newFace.reserve(face.size());

			for (auto& index : face)
			{
				auto found = vertexToIndex.find(index.position);

				if (found == vertexToIndex.end())
				{
					int32 newp = positions.size();

					vertexToIndex[index.position] = newp;

					positions.push_back(data.positions[index.position]);
					if (data.normals.size() > index.position)
						normals.push_back(data.normals[index.position]);
					else
						normals.push_back(Vector(0.0f, 0.0f, 1.0f));

					if( data.texcoords.size() > index.position ) 
						texcoords.push_back(data.texcoords[index.position]);
					else
						normals.push_back(Vector(0.0f, 0.0f, 0.0f));

					int32 newn = normals.size() == 0 ? -1 : normals.size() - 1;
					int32 newt = texcoords.size() == 0 ? -1 : texcoords.size() - 1;

					newFace.push_back(Index(newp, newn, newt));
				}
				else
				{
					int32 i = found->second;

					newFace.push_back(Index(i, i, i));
				}
			}


		}

		// nuke it, nothing here
		if (indices.size() == 0)
		{
			result.pop_back();
		}
	}

	return result;
}

void tcodsMeshInterfaceBase::_buildRuntimeMeshSection(URuntimeMeshComponent * runtimeMesh, tcods::Mesh& mesh, uint32_t section, FTransform transform, bool clipToLimits, FBox limits)
{
	// create a runtime static mesh directly from the tcods half-edge data structure so that our face indices are in sync
	TArray<FVector> runtimePositions;
	TArray<FVector2D> uvs;
	TArray<FVector> normals;
	TArray<int32> runtimeTriangles;

	FQuat rotateNormal = transform.GetRotation();

	std::array<FVector2D,3> uvs_base = { FVector2D(0.0,0.0), FVector2D(0.0, 1.0), FVector2D(1.0, 1.0) };

	// faces in tcods are ordered already by face.index
	int32 fi = 0;
	for (auto& face : mesh.faces)
	{
		auto halfEdge = face.he;



		// clip positions to limits
		// include triangles that have at least one vertex inside the limits
		int inLimitCount = 0;
		int vertexCount = 0;

		do
		{
			FVector position = unreal(halfEdge->from->position);

			bool inside = limits.IsInside(position);

			if (inside) inLimitCount++;

			halfEdge = halfEdge->next;

			vertexCount++;
		} while (halfEdge != face.he);

		// skip the triangle if the whole thing is outside the limits
		// if one vertex is in the limits, project outside vertices to the limits and include it
		if ((inLimitCount > 0 || !clipToLimits) && vertexCount == 3)
		{
			int ti = 0;
			do
			{
				FVector position = unreal(halfEdge->from->position);

				if (clipToLimits)
					position = limits.GetClosestPointTo(position);

				position = transform.TransformPosition(position);



				_halfEdgeVertex_to_runtimeMeshVertices[halfEdge->from->index].push_back(runtimePositions.Num());

				normals.Add(rotateNormal.RotateVector(unreal(halfEdge->normal)));
				runtimePositions.Add(position);
				uvs.Add(uvs_base[ti]);

				FVector2D texCoord(halfEdge->texcoord.x, halfEdge->texcoord.y);



				halfEdge = halfEdge->next;
				ti++;

			} while (halfEdge != face.he);

			check(ti <= 3);

			runtimeTriangles.Add(fi + 0);
			runtimeTriangles.Add(fi + 1);
			runtimeTriangles.Add(fi + 2);

			fi += 3;


			// we need flipped normals for the other way
			int32 base = runtimePositions.Num() - 1;

			for (int c = 0; c < 3; ++c)
			{
				auto n = -normals[base - c];
				normals.Add(n);

				auto p = runtimePositions[base - c];
				runtimePositions.Add(p);

				auto u = uvs[base - c];
				uvs.Add(u);
			}

			runtimeTriangles.Add(fi + 0);
			runtimeTriangles.Add(fi + 1);
			runtimeTriangles.Add(fi + 2);

			fi += 3;
		}
	}


	runtimeMesh->CreateMeshSection(section, runtimePositions, runtimeTriangles, normals, uvs, TArray<FColor>() /* colors */, TArray<FRuntimeMeshTangent>(), false, EUpdateFrequency::Frequent);
}

void tcodsMeshInterfaceBase::_buildRuntimeMesh(URuntimeMeshComponent * runtimeMesh, bool clipToLimits, FBox limits)
{
	using namespace tcods;

	runtimeMesh->ClearAllMeshSections();

	const FTransform transform = runtimeMesh->GetComponentTransform().Inverse();

	for( auto& pair : _sections )
	{
		Mesh& mesh = *pair.second.get();
		uint32_t section = pair.first;

		_buildRuntimeMeshSection(runtimeMesh, mesh, section, transform, clipToLimits, limits);
	}

}

// ---------------------------------------------------------
// tcodsUStaticMeshInterface
// ---------------------------------------------------------

auto tcodsUStaticMeshInterface::_toMeshData(UStaticMesh& uMesh, const FTransform toWorld) -> MeshDataAndVertexLookup
{
	using namespace tcods;

	MeshDataAndVertexLookup result;

	MeshIO::MeshData& data = result.meshData;
	std::unordered_map<FVector, VertexIndex>& vertexLookup = result.vertexLookup;


	if (uMesh.RenderData == nullptr) return result;
	if (uMesh.RenderData->LODResources.Num() == 0) return result;


	// build the TCODS half-edge data structure
	
	auto& positions = data.positions;
	auto& indices = data.indices;
	auto& normals = data.normals;
	auto& texcoords = data.texcoords;

	FStaticMeshLODResources& resource = uMesh.RenderData->LODResources[0 /* meshComponent->PreviousLODLevel */];

	FPositionVertexBuffer& positionBuffer = resource.PositionVertexBuffer;
	FRawStaticIndexBuffer& indexBuffer = resource.IndexBuffer;
	FStaticMeshVertexBuffer& vertexBuffer = resource.VertexBuffer;

	// now we need to process it with tcods
	using namespace tcods;
	using namespace DDG;

	uint32 vCount = positionBuffer.GetNumVertices();

	// Unreal stores duplicate vertices, we need to resolve this with a map
	std::vector<uint32> vertexIndices;
	vertexIndices.resize(vCount);

	// but each face has its own set of normals for each vertex in the face
	bool hasUvs = vertexBuffer.GetNumTexCoords() == vertexBuffer.GetNumVertices();

	for (uint32 i = 0; i < vCount; ++i)
	{
		const FVector& v = positionBuffer.VertexPosition(i);
		const auto& normal = vertexBuffer.VertexTangentZ(i);
		
		FVector2D uv = hasUvs ? vertexBuffer.GetVertexUV(i, 0) : FVector2D(0.0f, 0.0f);

		normals.push_back(to_tcods(normal));
		texcoords.push_back(Vector(uv.X, uv.Y, 0.0f));

		auto found = vertexLookup.find(v);
		if (found == vertexLookup.end())
		{
			vertexLookup[v].sectionIndex = 0;
			vertexLookup[v].vertexIndex = positions.size();

			vertexIndices[i] = positions.size();

			positions.push_back(to_tcods(toWorld.TransformPosition(v)));
		}
		else
		{
			vertexIndices[i] = found->second.vertexIndex;
		}
	}

	FIndexArrayView indexArrayView = indexBuffer.GetArrayView();
	int32 iCount = indexArrayView.Num();
	for (int32 i = 0; i + 2 < iCount; i += 3)
	{
		// let's assume triangles
		std::vector<MeshIO::Index> face(3);

		// need to consider the winding order wrt to the face normal
		const int32 indexArrayView0 = indexArrayView[i + 0];
		const int32 indexArrayView1 = indexArrayView[i + 1];
		const int32 indexArrayView2 = indexArrayView[i + 2];

		Vector n0 = normals[indexArrayView0];
		Vector n1 = normals[indexArrayView1];
		Vector n2 = normals[indexArrayView2];

		Vector n = (n0 + n1 + n2) / 3.0;

		int32 i0 = vertexIndices[indexArrayView0];
		int32 i1 = vertexIndices[indexArrayView1];
		int32 i2 = vertexIndices[indexArrayView2];

		face[0] = MeshIO::Index(i0, indexArrayView[i + 0], indexArrayView[i + 0]);
		face[1] = MeshIO::Index(i1, indexArrayView[i + 1], indexArrayView[i + 1]);
		face[2] = MeshIO::Index(i2, indexArrayView[i + 2], indexArrayView[i + 2]);

		indices.push_back(face);
	}

	return result;
}


auto tcodsUStaticMeshInterface::buildMesh(
	UStaticMesh& uMesh,
	const FTransform& toWorld,
	FBox worldLimits,
	URuntimeMeshComponent * runtimeMesh) -> bool
{
	using namespace tcods;

	_clear();

	_worldLimits = worldLimits;

	if(runtimeMesh) runtimeMesh->ClearAllMeshSections();

	MeshDataAndVertexLookup wholeMeshDataAndLookup = _toMeshData(uMesh, toWorld);

	MeshIO::MeshData& wholeData = wholeMeshDataAndLookup.meshData;
	_vertexLookup = wholeMeshDataAndLookup.vertexLookup;


	// process into sections
	std::vector<MeshIO::MeshData> extractedSectionData = _extractSections(wholeData);

	int32 sectionIndex = 0;

	for( MeshIO::MeshData& sectionData : extractedSectionData)
	{
		_sections[sectionIndex] = std::make_unique<Mesh>();
		Mesh& mesh = *_sections[sectionIndex].get();

		MeshIO::buildMesh(sectionData, mesh);

		mesh.topologyChange();
		mesh.geometryChange();


		sectionIndex++;
	}

	if( runtimeMesh ) _buildRuntimeMesh(runtimeMesh, false, worldLimits);

	_initBVH();
	_initNearestKDTree(worldLimits);

	return true;
}



// ---------------------------------------------------------
// tcodsChunkedGridMeshInterface
// ---------------------------------------------------------


void tcodsChunkedGridMeshInterface::buildMesh(ChunkGrid<float>& chunkGrid, const FTransform& toWorld, float isoLevel, FBox worldLimits, URuntimeMeshComponent * runtimeMesh)
{
	_clear();

	_worldLimits = worldLimits;

	auto mergedData = _toMergedMeshData(chunkGrid, isoLevel, toWorld);

	_buildSectionsFromMeshData(mergedData);

	if (runtimeMesh)
	{
		_buildRuntimeMesh(runtimeMesh, true, worldLimits);
	}

	_initBVH();
	_initNearestKDTree(worldLimits);
}

tcods::MeshIO::MeshData tcodsChunkedGridMeshInterface::_toTcods(SmoothMeshFactory& factory)
{
	using namespace tcods;
	using namespace DDG;

	// merge all of the mesh chunks into a single mesh data
	tcods::MeshIO::MeshData tcodsData;

	std::vector<Vector>& positions = tcodsData.positions;
	std::vector<Vector>& texcoords = tcodsData.texcoords;
	std::vector<Vector>& normals = tcodsData.normals;
	std::vector< std::vector< MeshIO::Index > >& indices = tcodsData.indices;

	// reserve
	positions.reserve(factory.vertices.Num());
	texcoords.reserve(factory.uvs.Num());
	normals.reserve(factory.normals.Num());
	indices.reserve(factory.indices.Num() / 3);

	// copy
	for (auto& v : factory.vertices)
		positions.push_back(to_tcods(v));

	for (auto& uv : factory.uvs)
		texcoords.emplace_back(uv.X, uv.Y, 0.0f);

	for (auto& n : factory.normals)
		normals.push_back(to_tcods(n));

	for (int32 i = 0; i < factory.indices.Num(); i += 3)
	{
		int32* triangle = &factory.indices[i];

		indices.emplace_back(3);
		auto& face = indices.back();

		face[0] = MeshIO::Index(triangle[0], triangle[0], triangle[0]);
		face[1] = MeshIO::Index(triangle[1], triangle[1], triangle[1]);
		face[2] = MeshIO::Index(triangle[2], triangle[2], triangle[2]);
	}

	return tcodsData;
}


auto tcodsChunkedGridMeshInterface::_toMergedMeshData(ChunkGrid<float> &chunkGrid, float isoLevel, FTransform toWorld) -> tcods::MeshIO::MeshData
{
	using namespace tcods;
	using namespace DDG;

	FVector2D uvScale(1.0f, 1.0f);

	FIntVector minChunkIndex = chunkGrid.minChunkndex();
	FIntVector maxChunkIndex = chunkGrid.maxChunkIndex();

	FIntVector minGridIndex = chunkGrid.gridIndexFromChunkIndex(minChunkIndex);
	FIntVector maxGridIndex = chunkGrid.gridIndexFromChunkIndex(maxChunkIndex) + chunkGrid.chunkDimensions() - FIntVector(1,1,1);

	const auto chunkGeometries = ChunkedMarchingCubes::marchingCubes_byChunk(isoLevel, chunkGrid, minGridIndex, maxGridIndex, uvScale, true, true);

	SmoothMeshFactory mergedFactory;

	// merge all of the mesh chunks into a single mesh data
	for (auto pair : chunkGeometries)
	{
		FIntVector chunkIndex = pair.first;
		ChunkedMarchingCubes::ChunkGeometry& chunkGeometry = pair.second;

		SmoothMeshFactory& meshFactory = chunkGeometry.factory;

		meshFactory.transform(toWorld);

		mergedFactory.append(meshFactory);
	}

	// remove duplicates
	mergedFactory = _clean(mergedFactory);
	//mergedFactory = mergedFactory.createMesh_mergeNearbyVertices(0.1f);
	//mergedFactory.removeDuplicateFaces();

	tcods::MeshIO::MeshData mergedData = _toTcods(mergedFactory);

	return mergedData;
}

void tcodsChunkedGridMeshInterface::_buildSectionsFromMeshData(tcods::MeshIO::MeshData& data)
{
	using namespace tcods;
	using namespace DDG;


	std::vector<MeshIO::MeshData> dataSections = _extractSections(data);

	// now build a mesh for each section
	uint32_t sectionIndex = 0;
	for (MeshIO::MeshData& sectionData : dataSections)
	{
		_sections[sectionIndex] = std::make_unique<Mesh>();

		Mesh& sectionMesh = *_sections[sectionIndex].get();

		int32 errorCode = MeshIO::buildMesh(sectionData, sectionMesh);
		if (errorCode > 0)
		{
			UE_LOG(LogTemp, Warning, TEXT("_buildMeshFromMeshData bad mesh. Section %d code %d"), sectionIndex, errorCode);

			bool hasIsolated = MeshIO::hasIsolatedVertices(sectionMesh);
			bool hasNonManifold = MeshIO::hasNonManifoldVertices(sectionMesh);

			UE_LOG(LogTemp, Warning, TEXT("hasIsolatedVertices %d hasNonManifoldVertices %d"), hasIsolated, hasNonManifold);
		}

		sectionMesh.topologyChange();
		sectionMesh.geometryChange();

		sectionIndex++;
	}
}

// The class prototypes.
using namespace vcg;
using namespace tri;

// forward declarations
class CFace;
class CVertex;
class CHEdge;
class CEdge;
class MyPolyVertex;

struct CUsedTypes : public vcg::UsedTypes< vcg::Use<CVertex>::AsVertexType, vcg::Use<CFace>::AsFaceType > {};

// Mesh of triangles
class CVertex : public Vertex<
	CUsedTypes,
	vertex::BitFlags,
	vertex::Coord3f,
	vertex::Normal3f,
	vertex::TexCoord2f,
	vertex::VFAdj
> {};

class CFace : public Face<
	CUsedTypes,
	face::VertexRef,
	face::Normal3f,
	face::BitFlags,
	face::FFAdj,
	face::VFAdj
> {};

class MyMesh : public vcg::tri::TriMesh< std::vector<CVertex>, std::vector<CFace> > {};


std::unique_ptr<MyMesh> _toMyMesh(SmoothMeshFactory& meshFactory)
{
	std::unique_ptr<MyMesh> result = std::make_unique<MyMesh>();

	MyMesh& myMesh = *result;

	int32 nv = meshFactory.vertices.Num();
	int32 ni = meshFactory.indices.Num();

	MyMesh::VertexIterator vi = vcg::tri::Allocator<MyMesh>::AddVertices(myMesh, nv);
	MyMesh::FaceIterator fi = vcg::tri::Allocator<MyMesh>::AddFaces(myMesh, ni / 3);

	std::vector<MyMesh::VertexPointer> vertexPointers;
	vertexPointers.reserve(meshFactory.vertices.Num());

	for (int i = 0; i < ni; ++i)
	{
		auto& vertex = meshFactory.vertices[i];
		auto& normal = meshFactory.normals.Num() > i ? meshFactory.normals[i] : FVector::ZeroVector;
		auto& texcoord = meshFactory.uvs.Num() > i ? meshFactory.uvs[i] : FVector2D::ZeroVector;

		vertexPointers.push_back(&*vi);

		vi->P() = MyMesh::CoordType(&vertex.X);
		vi->T() = vcg::TexCoord2f(texcoord.X, texcoord.Y);
		vi->N() = MyMesh::CoordType(&normal.X);

		++vi;
	}

	for (int i = 0; i < ni; i += 3)
	{
		int32* triangle = &meshFactory.indices[i];

		fi->V(0) = vertexPointers[triangle[0]];
		fi->V(1) = vertexPointers[triangle[1]];
		fi->V(2) = vertexPointers[triangle[2]];

		++fi;
	}

	return result;
}

FVector unreal(const vcg::Point3f& p)
{
	return FVector(p.X(), p.Y(), p.Z());
}

FVector2D unreal(const vcg::Point2f& p)
{
	return FVector2D(p.X(), p.Y());
}

FVector2D unreal(const vcg::TexCoord2f& p)
{
	return FVector2D(p.U(), p.V());
}

SmoothMeshFactory _toMeshFactory(MyMesh& mesh)
{
	SmoothMeshFactory factory;

	auto& vertices = factory.vertices;
	auto& normals = factory.normals;
	auto& uvs = factory.uvs;
	auto& indices = factory.indices;

	// reserve space
	vertices.Reserve(mesh.VN());
	normals.Reserve(mesh.VN());
	uvs.Reserve(mesh.VN());
	indices.Reserve(mesh.FN() * 3);

	// map past deleted vertices
	std::vector<int> toExportedVertexIndices(mesh.vert.size());

	// export vertices
	int vi = 0;
	for (auto& vertex : mesh.vert)
	{
		if (vertex.IsD())
		{
			vi++;
			continue;
		}

		toExportedVertexIndices[vi] = vertices.Num();

		vertices.Add(unreal(vertex.P()));
		normals.Add(unreal(vertex.N()));
		uvs.Add(unreal(vertex.T()));

		vi++;
	}

	// export face indices
	for (auto& fi : mesh.face)
	{
		if(fi.IsD())
			continue;

		auto vi0 = vcg::tri::Index(mesh, fi.V(0));
		auto vi1 = vcg::tri::Index(mesh, fi.V(1));
		auto vi2 = vcg::tri::Index(mesh, fi.V(2));

		indices.Add(toExportedVertexIndices[vi0]);
		indices.Add(toExportedVertexIndices[vi1]);
		indices.Add(toExportedVertexIndices[vi2]);
	}

	return factory;
}

SmoothMeshFactory tcodsChunkedGridMeshInterface::_clean(SmoothMeshFactory& mergedFactory)
{
 	typedef vcg::tri::Clean<MyMesh> Cleaner  ;

	auto myMesh = _toMyMesh(mergedFactory);

	int removedDegenerateFaces = Cleaner::RemoveDuplicateFace(*myMesh);
	int removedDuplicateFaces = Cleaner::RemoveDuplicateFace(*myMesh);
	int removedDuplicateVertices = Cleaner::RemoveDuplicateVertex(*myMesh);
	int removedUnreferencedVertices = Cleaner::RemoveUnreferencedVertex(*myMesh);
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(*myMesh);

	int removedNonManifoldFaces = Cleaner::RemoveNonManifoldFace(*myMesh);
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(*myMesh);

	int removedNonManifoldVertices = Cleaner::RemoveNonManifoldVertex(*myMesh);
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(*myMesh);

	removedUnreferencedVertices += Cleaner::RemoveUnreferencedVertex(*myMesh);
	vcg::tri::UpdateTopology<MyMesh>::FaceFace(*myMesh);

	vcg::tri::Allocator<MyMesh>::CompactEveryVector(*myMesh);

	return _toMeshFactory(*myMesh);
}





































// Keeping this code here, because it has some useful ideas
//void tcodsMeshInterfaceBase::_buildSamplePoints(FBox limits)
//{
//	_nearestPointSamples.pts.clear();
//
//	FBox bounds = FBox(EForceInit::ForceInitToZero);
//
//	for (auto& meshPair : _sections)
//	{
//		for (const tcods::Vertex& vertex : meshPair.second->vertices)
//		{
//			FVector fPosition = unreal(vertex.position);
//
//			bounds += fPosition;
//		}
//	}
//
//	FVector size = bounds.GetSize();
//	float max = size.GetMax();
//
//	float stepSize = max / sampleSubdivisions;
//
//	bvh::Ray ray;
//
//	std::vector<SamplePoint>& hits = _nearestPointSamples.pts;
//
//	for (unsigned int dim = 0; dim < 3; ++dim)
//	{
//		unsigned int dimI = (dim + 1) % 3;
//		unsigned int dimJ = (dim + 2) % 3;
//
//		ray.direction = FVector::ZeroVector;
//		ray.direction[dim] = 1.0f;
//
//		ray.origin[dim] = bounds.Min[dim];
//
//		for (float i = bounds.Min[dimI] - 2.0f * stepSize; i <= bounds.Max[dimI] + stepSize; i += stepSize)
//		{
//			ray.origin[dimI] = i;
//
//			for (float j = bounds.Min[dimJ] - 2.0f * stepSize; j < bounds.Max[dimJ] + stepSize; j += stepSize)
//			{
//				ray.origin[dimJ] = j;
//
//				auto intersections = bvh.getIntersections(ray, false);
//
//				for (auto& intersection : intersections)
//				{
//					const FVector p = intersection.t * ray.direction + ray.origin;
//
//					if (limits.IsValid && limits.IsInside(p))
//						hits.push_back({ p.X, p.Y, p.Z, intersection.triangle.surfaceIndex });
//				}
//			}
//		}
//	}
//}