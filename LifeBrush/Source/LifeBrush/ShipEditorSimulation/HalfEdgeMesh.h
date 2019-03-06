// Copyright 2016 Timothy Davison, all rights reserved.

#pragma once

#include <vector>
#include <unordered_map>

#include "Utility.h"
#include "MeshFactory.h"

#include "HalfEdgeMesh.generated.h"

struct FGraph;

struct FHalfEdgeVertex;
struct FHalfEdgeHullFace;
struct FHullHalfEdge;
struct HalfEdgeMesh;

typedef int32 VertexIndex;
typedef int32 FaceIndex;
typedef int32 HalfEdgeIndex;

// a reference to an element in a vector that is invariant to changes in the vectors size
// but that is variant to reordering
template<typename TElement>
struct ElementReference
{
public:
    int32 i = -1;
    std::vector<TElement>* v = nullptr;

public:
    ElementReference() {}
    ElementReference( ElementReference&& other )
    {
        i = other.i;
        v = other.v;
    }

    ElementReference(const ElementReference& other )
    {
        i = other.i;
        v = other.v;
    }

    ElementReference( int32 i_in, std::vector<TElement>* v_in ) : i( i_in ), v( v_in ) {}

    ElementReference<TElement>& operator=( const ElementReference<TElement>& other )
    {
        if(&other == this)
            return *this;

        i = other.i;
        v = other.v;

        return *this;
    }

    TElement* operator ->()
    {
        return &(*v)[i];
    }

    inline operator bool() const
    {
        return i >= 0;
    }

    inline bool operator==(const ElementReference<TElement>& other) const
    {
        return (i == other.i) && (v == other.v);
    }

    inline bool operator!=( const ElementReference<TElement>& other ) const
    {
        return !(*this == other);
    }
};

typedef ElementReference<FHalfEdgeVertex> VertexReference;
typedef ElementReference<FHullHalfEdge> HalfEdgeReference;
typedef ElementReference<FHalfEdgeHullFace> FaceReference;

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FHalfEdgeVertex
{
    GENERATED_USTRUCT_BODY()

    FHalfEdgeVertex() {} 

    FHalfEdgeVertex( FVector position_in ) : position( position_in ) {};

    void invalidate()
    {
        _out = -1;
        node = -1;
    }

    bool isValid() { return _out >= 0; }

    bool _deleted = false;

    void setDeleted( bool deleted ) { _deleted = deleted; }
    bool isDeleted() { return _deleted; }

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    FVector position;

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 node = -1; // the FGraphNode index that this vertex is attached to

    HalfEdgeReference out( HalfEdgeMesh& mesh );
    void setOut( HalfEdgeMesh& mesh, HalfEdgeReference& edge );

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _out = -1;

    template<typename TLambda>
    void visit( HalfEdgeMesh& mesh, TLambda&& visitor );
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FHullHalfEdge
{
    GENERATED_USTRUCT_BODY()

    FHullHalfEdge() {}

    FHullHalfEdge( FVector normal_in ) : normal( normal_in ) {};

    // disconnects the half-edge from everything in the mesh
    void invalidate()
    {
        _from = -1;
        _next = -1;
        _flip = -1;
        _face = -1;
        graphEdge = -1;
    }

    bool isValid() { return _next >= 0; }

    bool _deleted = false;

    void setDeleted( bool deleted ) { _deleted = deleted;  }
    bool isDeleted() { return _deleted; }
    bool isBoundary() { return _face < 0; }

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    FVector normal; // vertex normal for the from vertex
        
    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _from = -1; // the from vertex

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _next = -1; // the next edge in the face

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _flip = -1; // the sibling of the half-edge

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _face = -1; // the incident face

    VertexReference from( HalfEdgeMesh& mesh );
    void setFrom( HalfEdgeMesh& mesh, VertexReference& vertex );

    HalfEdgeReference next( HalfEdgeMesh& mesh );
    void setNext( HalfEdgeMesh& mesh, HalfEdgeReference& edge );

    HalfEdgeReference previous( HalfEdgeMesh& mesh );

    HalfEdgeReference flip( HalfEdgeMesh& mesh );
    void setFlip( HalfEdgeMesh& mesh, HalfEdgeReference& edge );

    FaceReference face( HalfEdgeMesh& mesh ); 
    void setFace( HalfEdgeMesh& mesh, FaceReference& face_in );

    int32 graphEdge = -1; // the FGraphEdge index that this edge is attached to
};

USTRUCT( BlueprintType )
struct LIFEBRUSH_API FHalfEdgeHullFace
{
    GENERATED_USTRUCT_BODY()

    FHalfEdgeHullFace() {}
    FHalfEdgeHullFace(int32 type_) : type(type_) {}

    void invalidate()
    {
        _halfEdge = -1;
    }

    bool isValid() { return _halfEdge >= 0; }

    bool _deleted = false;

    void setDeleted( bool deleted ) { _deleted = deleted; }
    bool isDeleted() { return _deleted; }

    bool isBoundary( HalfEdgeMesh& mesh );
    HalfEdgeReference boundary( HalfEdgeMesh& mesh );

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 _halfEdge = -1; // one of the half-edges in this face

    HalfEdgeReference halfEdge( HalfEdgeMesh& mesh );
    void setHalfEdge( HalfEdgeMesh& mesh, HalfEdgeReference& edge );

    UPROPERTY( EditAnywhere, BlueprintReadWrite, Category = "ShipEditor" )
    int32 type; // the hull type
};



struct LIFEBRUSH_API HalfEdgeMesh
{
public:
    std::vector<FHalfEdgeVertex> vertices;
    std::vector<FHullHalfEdge> halfEdges;
    std::vector<FHalfEdgeHullFace> faces; 
         
    template<class... Args>
    ElementReference<FHalfEdgeVertex> emplace_back_vertex( Args&&... args )
    {
        vertices.emplace_back( std::forward<Args>( args )... );
        return{ int32(vertices.size() - 1), &vertices };
    }

    template<class... Args>
    ElementReference<FHullHalfEdge> emplace_back_halfEdge( Args&&... args )
    {
        halfEdges.emplace_back( std::forward<Args>( args )...);
        return{ int32(halfEdges.size() - 1), &halfEdges };
    }

    template<class... Args>
    ElementReference<FHalfEdgeHullFace> emplace_back_face( VertexReference v0, VertexReference v1, VertexReference v2, Args&&... args );

    void removeVertex( int32 vi, FGraph& graph );
    void removeFace( int32 fi, FGraph& graph );

    HalfEdgeReference edgeBetween( VertexReference& a, VertexReference& b );

    std::unordered_map<int32, SmoothMeshFactory> buildMesh( FGraph& graph );        

    // Returns the clockwise-most half edge around the vertex that doesn't have a flip.
    // The half-edge will originate from the vertex (it's HalfEdge::from will be the vertex).
    HalfEdgeReference cwMostEdge( VertexReference vertex );

    // Returns the counter-clockwise-most half edge around the vertex that doesn't have a flip.
    // The half-edge will  point at the vertex (it's HalfEdge::from will not be the vertex).
    HalfEdgeReference ccwMostEdge( VertexReference& vertex );

    // Visit all of the half-edges whose from reference is the passed vertex.
    // The order of visitation is immune to modification of the mesh (it's calculated before being visited).
    template<typename TLambda>
    void visitOutgoing( VertexReference& vertex, TLambda&& visitor );

    template<typename TLambda>
    void visitFaceEdges( VertexReference& vertex, TLambda&& visitor );

    // Visit all the half-edges of the passed face.
    // The order of visitation is immune to modification of the mesh (it's calculated before being visited).
    template<typename TLambda>
    void visit( FaceReference& face, TLambda&& visitor );

    // Visit all of the half-edges in the loop starting from start.
    // the order of visitation is immune to modification of the mesh (it's calculated before being visited).
    template<typename TLambda>
    void visit( HalfEdgeReference& start, TLambda&& visitor );
private:
    void _eraseVertex( size_t vi, FGraph& graph );
    void _eraseVertices( std::vector<size_t> indices, FGraph& graph );

    void _eraseEdges( std::vector<size_t> indices, FGraph& graph );
    void _eraseFaces( std::vector<size_t> indices, FGraph& graph );

    void _popBackVertexTo( size_t vi, FGraph& graph );
    void _popBackEdgeTo( size_t ei, FGraph& graph );
    void _popBackFaceTo( size_t fi, FGraph& graph );

};


// Vertex Defintion
inline HalfEdgeReference FHalfEdgeVertex::out( HalfEdgeMesh& mesh ) { return{ _out, &mesh.halfEdges }; }
inline void FHalfEdgeVertex::setOut( HalfEdgeMesh& mesh, HalfEdgeReference& edge ) { _out = edge.i; }

template<class... Args>
ElementReference<FHalfEdgeHullFace>
HalfEdgeMesh::emplace_back_face( VertexReference v0, VertexReference v1, VertexReference v2, Args&&... args )
{
    faces.emplace_back( std::forward<Args>( args )... );

    FaceReference face( faces.size() - 1, &faces );

    HalfEdgeReference e[3];
    VertexReference v[3] = { v0, v1, v2 }; 

    for(int i = 0; i < 3; ++i)
        e[i] = edgeBetween( v[i], v[(i + 2) % 3] ); // (i - 1) % 3 

    for(int i = 0; i < 3; ++i)
    {
        HalfEdgeReference& edge = e[i];

        if(!edge)
        {
            edge = emplace_back_halfEdge();
            auto flip = emplace_back_halfEdge();
             
            edge->setFrom( *this, v[i] );
            flip->setFrom( *this, v[(i + 2) % 3] ); // (i - 1) % 3

            edge->setFlip( *this, flip );
            flip->setFlip( *this, edge );
        }

        if( v[i]->_out < 0 )
            v[i]->setOut( *this, edge );
    }

    for(int i = 0; i < 3; ++i)
    {
        e[i]->setNext( *this, e[(i + 2) % 3]); // (i - 1) % 3
        e[i]->setFace( *this, face ); 
    }

    for(int i = 0; i < 3; ++i)
    {
        auto& vertex = v[i];

        // loop ccw to find the boundary edge
        const auto outStart = vertex->out( *this );
        auto out = outStart;

        do
        {
            if(out->isBoundary() )
                break;

            out = out->next( *this )->next( *this )->flip( *this );
        } while( outStart != out );

        // loop cw to find the incoming boundary edges
        // connect them
        const auto inStart = vertex->out( *this )->flip(*this);
        auto in = inStart;

        do 
        {
            if(in->isBoundary())
                break;

            in = in->next( *this )->flip( *this );
        } while( in != inStart );

        in->setNext( *this, out );
    }

    face->setHalfEdge( *this, e[0] ); 

    return face;
}

template<typename TLambda>
inline void HalfEdgeMesh::visitFaceEdges( VertexReference& vertex, TLambda && visitor )
{
    const auto vertexOut = vertex->out( *this );
    auto edge = vertexOut;

    if(!vertexOut)
        return;

    // do the work in an array so that we can allow the visitor to mutate the graph
    std::vector<HalfEdgeReference> toVisit;

    // loop cw
    do 
    {
        if(!edge->isBoundary())
            toVisit.push_back( edge );

        edge = edge->flip( *this )->next( *this );
    } while(vertexOut != edge );

    // visit cw
    for(auto& edgeToVisit : toVisit)
        visitor( edgeToVisit );
}

template<typename TLambda>
inline void HalfEdgeMesh::visitOutgoing( VertexReference& vertex, TLambda && visitor )
{
    const auto vertexOut = vertex->out( *this );
    auto edge = vertexOut;

    if(!vertexOut)
        return;

    // do the work in an array so that we can allow the visitor to mutate the graph
    std::vector<HalfEdgeReference> toVisit;

    // loop ccw
    do
    {
        toVisit.push_back( edge );

        edge = edge->flip( *this )->next( *this );
    } while(vertexOut != edge);

    // visit ccw
    for(auto& edgeToVisit : toVisit)
        visitor( edgeToVisit );
}

template<typename TLambda>
void HalfEdgeMesh::visit( FaceReference& face, TLambda&& visitor )
{
    // do the work in an array so that we can allow the visitor to mutate the graph
    std::vector<HalfEdgeReference> toVisit;

    auto start = face->halfEdge( *this );
    auto edge = start;

    if(!start)
        return;

    do 
    {
        toVisit.push_back( edge );

        edge = edge->next( *this );
    } while (edge != start);

    size_t n = toVisit.size();
    for(int i = 0; i < n; ++i)
    {
        visitor( toVisit[i], toVisit[(n + i - 1) % n] );
    }
}



template<typename TLambda>
void HalfEdgeMesh::visit( HalfEdgeReference& start, TLambda&& visitor )
{
    if(!start) 
        return;

    // do the work in an array so that we can allow the visitor to mutate the graph
    std::vector<HalfEdgeReference> toVisit;

    auto edge = start;
    do 
    {
        toVisit.push_back( edge );

        edge = edge->next( *this );
    } while (edge != start);

    for(auto edgeToVisit : toVisit)
        visitor( edgeToVisit );
}

// HalfEdge Definition
inline VertexReference FHullHalfEdge::from( HalfEdgeMesh& mesh ) { return{ _from, &mesh.vertices }; }
inline void FHullHalfEdge::setFrom( HalfEdgeMesh& mesh, VertexReference& vertex ) { _from = vertex.i; }

inline HalfEdgeReference FHullHalfEdge::next( HalfEdgeMesh& mesh ) { return{ _next, &mesh.halfEdges }; }
inline void FHullHalfEdge::setNext( HalfEdgeMesh& mesh, HalfEdgeReference& edge_in ) { _next = edge_in.i; }

inline HalfEdgeReference FHullHalfEdge::previous( HalfEdgeMesh& mesh ) 
{
    HalfEdgeReference start = next( mesh );
    HalfEdgeReference prev;
    HalfEdgeReference edge;
    
    do 
    {
        auto nextEdge = edge->next( mesh );

        if(nextEdge == start)
            return prev;
        
        prev = edge;
        edge = edge->next( mesh );
    } while (edge != start);

    return HalfEdgeReference();
}

inline HalfEdgeReference FHullHalfEdge::flip( HalfEdgeMesh& mesh ) { return{ _flip, &mesh.halfEdges }; }
inline void FHullHalfEdge::setFlip( HalfEdgeMesh& mesh, HalfEdgeReference& edge ) { _flip = edge.i; }

inline FaceReference FHullHalfEdge::face( HalfEdgeMesh& mesh ) { return{ _face, &mesh.faces }; }
inline void FHullHalfEdge::setFace( HalfEdgeMesh& mesh, FaceReference& face_in ) { _face = face_in.i; }

// Face Definition
inline HalfEdgeReference FHalfEdgeHullFace::halfEdge( HalfEdgeMesh& mesh ) { return{ _halfEdge, &mesh.halfEdges }; }
inline void FHalfEdgeHullFace::setHalfEdge( HalfEdgeMesh& mesh, HalfEdgeReference& edge ) { _halfEdge = edge.i; }

struct HullMeshArchiveProxy
{
    static FArchive& archive( FArchive& archive, HalfEdgeMesh& mesh )
    {
        if(archive.IsSaving()) 
        {
            size_t count = mesh.vertices.size();
            archive << count;

            UScriptStruct * scriptStruct = FHalfEdgeVertex::StaticStruct();
            for(auto& v : mesh.vertices)
            {
                scriptStruct->SerializeItem( archive, &v, nullptr );
            }

            count = mesh.halfEdges.size();
            archive << count;

            scriptStruct = FHullHalfEdge::StaticStruct();
            for(auto& he : mesh.halfEdges)
            {
                scriptStruct->SerializeItem( archive, &he, nullptr );
            }

            count = mesh.faces.size();
            archive << count;

            scriptStruct = FHalfEdgeHullFace::StaticStruct();
            for(auto& face : mesh.faces)
            {
                scriptStruct->SerializeItem( archive, &face, nullptr );
            }
        }
        else
        {
            size_t count;
            archive << count;
            mesh.vertices.resize( count );

            UScriptStruct * scriptStruct = FHalfEdgeVertex::StaticStruct();
            for(auto& v : mesh.vertices)
            {
                scriptStruct->SerializeItem( archive, &v, nullptr );
            }

            archive << count;
            mesh.halfEdges.resize( count );

            scriptStruct = FHullHalfEdge::StaticStruct();
            for(auto& he : mesh.halfEdges)
            {
                scriptStruct->SerializeItem( archive, &he, nullptr );
            }

            archive << count;
            mesh.faces.resize( count );

            scriptStruct = FHalfEdgeHullFace::StaticStruct();
            for(auto& face : mesh.faces)
            {
                scriptStruct->SerializeItem( archive, &face, nullptr );
            }
        }

        return archive;
    }
};

FORCEINLINE FArchive& operator<<( FArchive& archive, HalfEdgeMesh& mesh )
{
    return HullMeshArchiveProxy::archive( archive, mesh );   
}