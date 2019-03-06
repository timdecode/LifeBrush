// Copyright 2016 Timothy Davison, all rights reserved.

#include "LifeBrush.h"

#include "HalfEdgeMesh.h"
#include "Graph.h"

void HalfEdgeMesh::removeVertex( int32 vi_, FGraph& graph )
{
    size_t vi = size_t( vi_ ); 

    VertexReference vertex( vi_, &vertices );

    if(!vertex->isDeleted())
        return;


}

void HalfEdgeMesh::removeFace( int32 faceIndex, FGraph& graph )
{
    FaceReference face( faceIndex, &faces ); 

    std::vector<size_t> edgesToRemove;

    visit( face, [&]( HalfEdgeReference& edge, HalfEdgeReference& previous ) {
        // disconnect any vertices from the loop
        VertexReference from = edge->from( *this );

        auto flip = edge->flip( *this );

        if(from->out( *this ) == previous)
        {
            if(flip)
                from->setOut( *this, flip );
            else
            {
                flip = previous->flip( *this );

                if(flip)
                {
                    auto flipPrevious = flip->next( *this )->next( *this );
                    from->setOut( *this, flipPrevious );
                }
                // this vertex is dead
                else
                {
                    _eraseVertex( from.i, graph );
                }
            }
        }

        // break the flip
        if(flip)
            flip->_flip = -1;

        edgesToRemove.push_back( edge.i );
    } );
    _eraseEdges( edgesToRemove, graph );
}

HalfEdgeReference HalfEdgeMesh::edgeBetween( VertexReference& a, VertexReference& b )
{
    HalfEdgeReference found;

    visitOutgoing( a, [&]( HalfEdgeReference& edge )
    {
        if(edge->next( *this )->from( *this ) == b)
            found = edge;
    } );

    return found;
}

void HalfEdgeMesh::_eraseVertex( size_t vi, FGraph& graph )
{
    std::vector<size_t> toErase;  
    toErase.push_back( vi );

    _eraseVertices( toErase, graph );
}


void HalfEdgeMesh::_eraseVertices(std::vector<size_t> indices, FGraph& graph )
{
    for(auto vi : indices)
    {
        auto& vertex = vertices[vi];
     
        if( !vertex.isValid() ) 
            continue;

        FGraphNode& node = graph.node( vertex.node );

        vertex.invalidate();
    }
}

void HalfEdgeMesh::_eraseEdges( std::vector<size_t> indices, FGraph& graph )
{
    for(auto i : indices)
    {
        auto& edge = halfEdges[i];

        if(!edge.isValid())
            continue;

        edge.invalidate();
    }
}

void HalfEdgeMesh::_eraseFaces( std::vector<size_t> indices, FGraph& graph )
{
    for(auto i : indices)
    {
        auto& face = faces[i];

        if(!face.isValid())
            continue;

        face.invalidate();
    }
}

std::unordered_map<int32, SmoothMeshFactory> HalfEdgeMesh::buildMesh( FGraph& graph )
{
    std::unordered_map<int32, SmoothMeshFactory> factories;

    for(auto& face : faces)
    {
        if( !face.isValid() )
            continue;

        auto& factory = factories[face.type];

        auto e0 = face.halfEdge( *this );
        auto v0 = e0->from( *this );

        auto e1 = e0->next( *this );
        do
        {
            auto e2 = e1->next( *this );

            auto v1 = e1->from( *this );
            auto v2 = e2->from( *this );

            factory.pushTriangle(
                v0->position, v1->position, v2->position,
                e0->normal, e1->normal, e2->normal
            ); 

            FGraphNode& n0 = graph.node( v0->node );
            FGraphNode& n1 = graph.node( v1->node ); 
            FGraphNode& n2 = graph.node( v2->node );

            factory.pushTriangle(
                n2.position, n1.position, n0.position,
                -e2->normal, -e1->normal, -e0->normal 
            );

            e1 = e1->next( *this );
        } while(e1 != face.halfEdge( *this ));
    }

    return factories;
}


// Returns the counter-clockwise-most half edge around the vertex that doesn't have a flip.
// The half-edge will  point at the vertex (it's HalfEdge::from will not be the vertex).
HalfEdgeReference HalfEdgeMesh::ccwMostEdge( VertexReference& vertex )
{

    HalfEdgeReference start = vertex->out( *this );
    HalfEdgeReference edge = start;

    do
    {
        HalfEdgeReference nextNext = edge->next( *this )->next( *this );
        HalfEdgeReference flip = nextNext->flip( *this );

        if(!flip)
            return nextNext;

        edge = flip;

    } while(start != edge);

    // there isn't a ccw-most half-edge (this is a closed loop)
    return HalfEdgeReference();
}

// Returns the clockwise-most half edge around the vertex that doesn't have a flip.
// The half-edge will originate from the vertex (it's HalfEdge::from will be the vertex).
HalfEdgeReference HalfEdgeMesh::cwMostEdge( VertexReference vertex )
{

    HalfEdgeReference start = vertex->out( *this );
    HalfEdgeReference edge = start;

    do
    {
        HalfEdgeReference flip = edge->flip( *this );

        if(!flip)
            return edge;

        edge = flip->next( *this );

    } while(start != edge);

    // there isn't a cw-most half-edge (this is a closed loop)
    return HalfEdgeReference();
}

bool FHalfEdgeHullFace::isBoundary( HalfEdgeMesh& mesh )
{
    return mesh.halfEdges[_halfEdge].flip( mesh )->isBoundary();
}

HalfEdgeReference FHalfEdgeHullFace::boundary( HalfEdgeMesh& mesh )
{
    return mesh.halfEdges[_halfEdge].flip( mesh );
}
