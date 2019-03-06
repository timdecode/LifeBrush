//============================================================
// HalfEdge.cpp
// Keenan Crane
//

#include "LifeBrush.h"

#include <map>
#include <unordered_set>
#include <cmath>
#include <queue>
#include <float.h>
#include <assert.h>
#include "HalfEdge.h"
#include "Connection.h"
#include "MeshIO.h"

#include <string>
#include <fstream>
#include <sstream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;
using namespace DDG;

namespace tcods
{
    int Vertex :: valence( void ) const
    {
        HalfEdgeIter he = out;
        int n = 0;
        
        do
        {
            n++;
            he = he->flip->next;
        }
        while( he != out );
        
        return n;
    }
    
    bool Vertex :: onBoundary( void ) const
    {
		// visit each edge out of this vertex
        HalfEdgeCIter he = out;
        do
        {
            if( he->onBoundary )
            {
                return true;
            }
            he = he->flip->next;
        }
        while( he != out );
        
        return false;
    }
    
    double Vertex :: defect( void ) const
    {
        double sum = 0.;
        
        // iterate over incident triangles
        HalfEdgeIter he = out;
        do
        {
            // grab vertices
            Vector p1 = he->from->position;
            Vector p2 = he->next->from->position;
            Vector p3 = he->next->next->from->position;
            
            // subtract incident angle from sum
            Vector u1 = ( p2 - p1 );
            Vector u2 = ( p3 - p1 );
            sum += atan2( cross(u1,u2).norm(), dot(u1,u2) );
            
            he = he->flip->next;
        }
        while( he != out );
        
        return 2.*M_PI - sum;
    }

	vector<HalfEdge> isolated; // all isolated vertices point to isolated.begin()

	bool Vertex::isIsolated(void) const
	{
		return out == isolated.begin();
	}


    void Edge :: updateStar( void )
    {
        double sum = 0.;
        
        // compute the dual/primal length ratio by adding up the
        // cotangents of the two angles opposing the current edge
        HalfEdgeCIter h = he;
        do
        {
            Vector& a( h->from->position );
            Vector& b( h->next->from->position );
            Vector& c( h->next->next->from->position );
            
            // compute the cotangent of the angle of the current
            // triangle that opposes the current edge
            Vector u = a-c;
            Vector v = b-c;
            double cotTheta = dot( u,v )/cross( u,v ).norm();
            
            sum += .5*cotTheta;
            
            h = h->flip;
        }
        while( h != he );
        
        star = max( sum, 0. );
    }
    
    bool Edge :: onBoundary( void ) const
    {
        return he->onBoundary || he->flip->onBoundary;
    }
    
    double Face :: area( void ) const
    {
        const Vector& a( he->from->position );
        const Vector& b( he->next->from->position );
        const Vector& c( he->next->next->from->position );
        
        return .5 * cross((b-a),(c-a)).norm();
    }
    
    Vector Face :: normal( void ) const
    {
        const Vector& a( he->from->position );
        const Vector& b( he->next->from->position );
        const Vector& c( he->next->next->from->position );
        
        return cross((b-a),(c-a)).unit();
    }

	Vector Face::smoothNormal( Vector p )
	{
		const Vector& a( he->from->position );
		const Vector& b( he->next->from->position );
		const Vector& c( he->next->next->from->position );

		const Vector v0 = b - a;
		const Vector v1 = c - a;
		const Vector v2 = p - a;

		float d00 = dot(v0 , v0);
		float d01 = dot(v0 , v1);
		float d11 = dot(v1 , v1);
		float d20 = dot(v2 , v0);
		float d21 = dot(v2 , v1);

		float denom = d00 * d11 - d01 * d01;

		float alpha = (d11 * d20 - d01 * d21) / denom;
		float beta = (d00 * d21 - d01 * d20) / denom;
		float u = 1.0f - alpha - beta;

		const Vector& na = he->normal;
		const Vector& nb = he->next->normal;
		const Vector& nc = he->next->next->normal;

		const Vector n = alpha * na + beta * nb + u * nc;

		return n;
	}

	// \see https://blogs.msdn.microsoft.com/rezanour/2011/08/07/barycentric-coordinates-and-point-in-triangle-tests/
	bool Face::contains( Vector p ) const
	{
		const Vector& a( he->from->position );
		const Vector& b( he->next->from->position );
		const Vector& c( he->next->next->from->position );

		const Vector u = b - a;
		const Vector v = c - a;
		const Vector w = p - a;

		const Vector v_x_w = cross(v,w);
		const Vector v_x_u = cross(v,u);

		if(dot(v_x_w , v_x_u) < 0.0)
			return false;

		const Vector u_x_w = cross(u,w);
		const Vector u_x_v = cross(u,v);

		if(dot(u_x_w , u_x_v) < 0.0)
			return false;

		double denominator = u_x_v.norm();

		float r = v_x_w.norm() / denominator;
		float t = u_x_w.norm() / denominator;

		return (r + t <= 1.0);
	}
    
    double Face :: circumradius( void ) const
    {
        const Vector& a( he->from->position );
        const Vector& b( he->next->from->position );
        const Vector& c( he->next->next->from->position );
        double u = (a-b).norm();
        double v = (b-c).norm();
        double w = (c-a).norm();
        
        return (u*v*w)/sqrt((u+v-w)*(u-v+w)*(-u+v+w)*(u+v+w));
    }
    
    Vector Face :: barycenter( void ) const
    {
        const Vector& a( he->from->position );
        const Vector& b( he->next->from->position );
        const Vector& c( he->next->next->from->position );
        
        return (a+b+c)/3.;
    }
    
    void Face :: frame( Vector& e1, Vector& e2 ) const
    {
        const Vector& a( he->from->position );
        const Vector& b( he->next->from->position );
        const Vector& c( he->next->next->from->position );
        
        e1 = ( b - a ).unit();
        e2 = c - a;
        e2 = ( e2 - dot(e2,e1)*e1 ).unit();
    }
    
    Vector Face :: toLocal( Vector p ) const
    {
        Vector origin = he->from->position;
        Vector e1, e2; frame( e1, e2 );
        
        p -= origin;
        
        return Vector( dot(e1,p), dot(e2,p), 0. );
    }
    
    Vector Face :: toGlobal( Vector q ) const
    {
        Vector origin = he->from->position;
        Vector e1, e2; frame( e1, e2 );
        
        return origin + q.x*e1 + q.y*e2;
    }
    
    bool Face :: isBoundary( void ) const
    {
        return he->onBoundary;
    }
    
    Mesh :: Mesh( const Mesh& mesh )
    {
        *this = mesh;
    }
    
    Mesh :: ~Mesh( void )
    {
        if( connection != nullptr )
        {
            delete connection;
        }
    }
    
    // iterator comparison operators -- needed to build STL maps of iterators in Mesh::operator=
    class HalfEdgeCIterCompare { public: bool operator()( const HalfEdgeCIter& i, const HalfEdgeCIter& j ) const { return &*i < &*j; } };
    class   VertexCIterCompare { public: bool operator()( const   VertexCIter& i, const   VertexCIter& j ) const { return &*i < &*j; } };
    class     FaceCIterCompare { public: bool operator()( const     FaceCIter& i, const     FaceCIter& j ) const { return &*i < &*j; } };
    class     EdgeCIterCompare { public: bool operator()( const     EdgeCIter& i, const     EdgeCIter& j ) const { return &*i < &*j; } };
    
    const Mesh& Mesh :: operator=( const Mesh& mesh )
    {
        map< HalfEdgeCIter, HalfEdgeIter, HalfEdgeCIterCompare > halfedgeOldToNew;
        map<   VertexCIter,   VertexIter,   VertexCIterCompare >   vertexOldToNew;
        map<     EdgeCIter,     EdgeIter,     EdgeCIterCompare >     edgeOldToNew;
        map<     FaceCIter,     FaceIter,     FaceCIterCompare >     faceOldToNew;
        
        // copy geometry from the original mesh and create a
        // map from pointers in the original mesh to
        // those in the new mesh
        halfedges.clear(); for( HalfEdgeCIter i = mesh.halfedges.begin(); i != mesh.halfedges.end(); i++ ) halfedgeOldToNew[ i ] = halfedges.insert( halfedges.end(), *i );
        vertices.clear(); for(   VertexCIter i =  mesh.vertices.begin(); i !=  mesh.vertices.end(); i++ )   vertexOldToNew[ i ] =  vertices.insert(  vertices.end(), *i );
        edges.clear(); for(     EdgeCIter i =     mesh.edges.begin(); i !=     mesh.edges.end(); i++ )     edgeOldToNew[ i ] =     edges.insert(     edges.end(), *i );
        faces.clear(); for(     FaceCIter i =     mesh.faces.begin(); i !=     mesh.faces.end(); i++ )     faceOldToNew[ i ] =     faces.insert(     faces.end(), *i );
        
        // ``search and replace'' old pointers with new ones
        for( HalfEdgeIter i = halfedges.begin(); i != halfedges.end(); i++ )
        {
            i->next = halfedgeOldToNew[ i->next ];
            i->flip = halfedgeOldToNew[ i->flip ];
            i->from =   vertexOldToNew[ i->from ];
            i->edge =     edgeOldToNew[ i->edge ];
            i->face =     faceOldToNew[ i->face ];
        }
        
        for( VertexIter i = vertices.begin(); i != vertices.end(); i++ )
        {
            i->out = halfedgeOldToNew[ i->out ];
        }
        
        for( EdgeIter i = edges.begin(); i != edges.end(); i++ ) i->he = halfedgeOldToNew[ i->he ];
        for( FaceIter i = faces.begin(); i != faces.end(); i++ ) i->he = halfedgeOldToNew[ i->he ];
        
        return *this;
    }
    
    void Mesh :: indexElements( void )
    {
        int i;
        
        index2he.resize( halfedges.size());
        index2vertex.resize( vertices.size());
        index2edge.resize( edges.size());
        index2face.resize( faces.size());
        
        i = 0; for( HalfEdgeIter he = halfedges.begin(); he != halfedges.end(); he++ ) {     index2he[i]=he; he->index=i++; }
        i = 0; for(   VertexIter  v =  vertices.begin();  v !=  vertices.end();  v++ ) { index2vertex[i]=v;   v->index=i++; }
        i = 0; for(     EdgeIter  e =     edges.begin();  e !=     edges.end();  e++ ) {   index2edge[i]=e;   e->index=i++; }
        i = 0; for(     FaceIter  f =     faces.begin();  f !=     faces.end();  f++ ) {   index2face[i]=f;   f->index=i++; }
    }
    
    HalfEdgeIter Mesh :: halfedge( int index )
    {
        return index2he[ index ];
    }
    
    VertexIter Mesh :: vertex( int index )
    {
        return index2vertex[ index ];
    }
    
    EdgeIter Mesh :: edge( int index )
    {
        return index2edge[ index ];
    }
    
    FaceIter Mesh :: face( int index )
    {
        return index2face[ index ];
    }
    
    HalfEdgeCIter Mesh :: halfedge( int index ) const
    {
        return index2he[ index ];
    }
    
    VertexCIter Mesh :: vertex( int index ) const
    {
        return index2vertex[ index ];
    }
    
    EdgeCIter Mesh :: edge( int index ) const
    {
        return index2edge[ index ];
    }
    
    FaceCIter Mesh :: face( int index ) const
    {
        return index2face[ index ];
    }
    
    int Mesh :: nFaces( void ) const
    {
        int n = 0;
        
        for( FaceCIter f = faces.begin(); f != faces.end(); f++ )
        {
            if( !f->he->onBoundary )
            {
                n++;
            }
        }
        
        return n;
    }
    
    double Mesh :: parallelTransport( double phi, HalfEdgeCIter he )
    // given an angle phi relative to the canonical reference frame
    // of he->face, returns the angle parallel transported across he
    // using the Levi-Civita connection, expressed relative to the
    // canonical frame of he->flip->face
    {
        // get (oriented) direction along shared edge
        VertexIter u = he->from;
        VertexIter v = he->flip->from;
        Vector e = v->position - u->position;
        if( u->index > v->index ) e = -e;
        
        // compute angle adjustments between canonical frames
        Vector e1, e2; he->face->frame( e1, e2 );
        Vector f1, f2; he->flip->face->frame( f1, f2 );
        double deltaIJ = atan2( dot(e,e2), dot(e,e1) );
        double deltaJI = atan2( dot(e,f2), dot(e,f1) );
        
        // transport phi
        return ( phi - deltaIJ ) + deltaJI;
    }

	double Mesh::parallelTransport_2( double phi, HalfEdgeCIter he )
		// given an angle phi relative to the canonical reference frame
		// of he->face, returns the angle parallel transported across he
		// using the Levi-Civita connection, expressed relative to the
		// canonical frame of he->flip->face
		// Tim: this version doesn't flip e based on vertex index
	{
		// get (oriented) direction along shared edge
		VertexIter u = he->from;
		VertexIter v = he->flip->from;
		Vector e = v->position - u->position;

		// compute angle adjustments between canonical frames
		Vector e1, e2; he->face->frame( e1, e2 );
		Vector f1, f2; he->flip->face->frame( f1, f2 );
		double deltaIJ = atan2( dot(e,e2), dot(e,e1) );
		double deltaJI = atan2( dot(e,f2), dot(e,f1) );

		// transport phi
		return (phi - deltaIJ) + deltaJI;
	}
    
    double Mesh :: defect( const Cycle& c )
    {
        double theta = 0.;
        
        for( Cycle::const_iterator he = c.begin(); he != c.end(); he++ )
        {
            theta = parallelTransport( theta, *he );
        }
        
        while( theta >=  M_PI ) theta -= 2.*M_PI;
        while( theta <  -M_PI ) theta += 2.*M_PI;
        
        return -theta;
    }
    
    void Mesh :: integralCurve( FaceIter initialFace,
                               const Vector& initialPoint,
                               double initialAngle,
                               std::vector<Vector>& curve,
                               std::vector<Vector>& normals,
                               int maxPts )
    {
        static int curveIndex = 0; // unique ID for each curve generated
        FaceIter f = initialFace; // current face
        Vector x = initialPoint; // current point
        double alpha = f->alpha + initialAngle;
        Vector u( cos(alpha), sin(alpha), 0. ); // current direction
        
        curve.clear();
        curve.push_back( f->toGlobal( x ));
        normals.push_back( f->normal() );
        
        for( int i = 0; i < maxPts; i++ )
        {
            // stop if we enter a virtual boundary face
            if( f->isBoundary() ) break;
            
            // stop if we've already visited this face
            if( f->curveIndex == curveIndex ) break;
            
            HalfEdgeCIter h[3];
            Vector q[3];
            
            h[0] = f->he;
            h[1] = f->he->next;
            h[2] = f->he->next->next;
            for( int j = 0; j < 3; j++ )
            {
                q[j] = f->toLocal( h[j]->from->position );
            }
            
            // intersect ray with each edge
            const double eps = 1e-6; // tolerance for intersecting the edge we're coming from
            double tMin = DBL_MAX; // minimum distance to any edge
            HalfEdgeCIter hMin = h[0]; // closest edge
            for( int j = 0; j < 3; j++ )
            {
                int k = (j+1)%3;
                Vector c = q[k]-q[j];
                c = Vector( -c.y, c.x, 0. ).unit();
                double d = dot(-c,q[j]);
                double t = -(dot(c,x)+d)/dot(c,u);
                
                if( t > eps && t < tMin )
                {
                    tMin = t;
                    hMin = h[j];
                }
            }
            
            // stop if intersection test yields dubious results
            if( tMin == DBL_MAX || tMin > 2.*f->circumradius() )
            {
                break;
            }
            
            // move current point to intersection point
            x += tMin*u;
            
            curve.push_back( f->toGlobal( x ));
            normals.push_back( f->normal() );
            
            // get pointer to next triangle
            FaceIter g = hMin->flip->face;
            
            // rewrite x in local coords of next triangle
            x = g->toLocal( f->toGlobal( x ));
            
            // rewrite u in local coords of next triangle
            double thetaIJ = hMin->edge->theta;
            if( hMin->from->index > hMin->flip->from->index ) thetaIJ = -thetaIJ;
            double beta = parallelTransport( atan2(u.y,u.x), hMin ) - thetaIJ;
            u = Vector( cos(beta), sin(beta), 0. );
            
            // indicate which curve last passed through f
            f->curveIndex = curveIndex;
            
            // move to next triangle
            f = g;
        }
        
        curveIndex++; // keep track of which curve we're integrating
    }

	std::pair<FaceIter, Vector>
	Mesh::integralWalk
	( 
		FaceIter initialFace,
		const Vector& initialPoint,
		double initialAngle,
		int maxPts,
		double distanceToTravel 
	)
	{
        // tim: we don't use face->curveIndex because it's not thread-safe
        // therefore, we use an unordered map to track visited faces, still pretty fast
		//static int curveIndex = 0; // unique ID for each curve generated
        std::unordered_set<int> visited;

		FaceIter f = initialFace; // current face
		Vector x = initialPoint; // current point
		double alpha = f->alpha + initialAngle;
		Vector u( cos( alpha ), sin( alpha ), 0. ); // current direction
		double distance = distanceToTravel;

		for(int i = 0; i < maxPts; i++)
		{
			// stop if we enter a virtual boundary face
			if(f->isBoundary()) break;

			// stop if we've already visited this face
            if(visited.find( f->index ) != visited.end())
                break;

			HalfEdgeCIter h[3];
			Vector q[3];

			h[0] = f->he;
			h[1] = f->he->next;
			h[2] = f->he->next->next;
			for(int j = 0; j < 3; j++)
			{
				q[j] = f->toLocal( h[j]->from->position );
			}

			// intersect ray with each edge
			const double eps = 1e-6; // tolerance for intersecting the edge we're coming from
			double tMin = DBL_MAX; // minimum distance to any edge
			HalfEdgeCIter hMin = h[0]; // closest edge
			for(int j = 0; j < 3; j++)
			{
				int k = (j + 1) % 3;
				Vector c = q[k] - q[j];
				c = Vector( -c.y, c.x, 0. ).unit();
				double d = -dot(c,q[j]);
				double t = -(dot(c,x) + d) / (dot(c,u));

				if(t > eps && t < tMin)
				{
					tMin = t;
					hMin = h[j];
				}
			}

			// stop if intersection test yields dubious results
			if(tMin == DBL_MAX || tMin > 2.*f->circumradius())
			{
				break;
			}

			// are we tired yet? (have we finished walking)
			if(distance - tMin < 0.0)
			{
				x += distance*u;
                
				return { f, x };
			}

			// move current point to intersection point
			x += tMin*u;
            distance -= tMin;

			Vector xGlobal = f->toGlobal( x );

			// get pointer to next triangle
			FaceIter g = hMin->flip->face;

			// rewrite x in local coords of next triangle
			x = g->toLocal( xGlobal );

			// rewrite u in local coords of next triangle
			double thetaIJ = hMin->edge->theta;
			//if(hMin->from->index > hMin->flip->from->index) thetaIJ = -thetaIJ;

			double beta = parallelTransport_2( atan2( u.y, u.x ), hMin ) - thetaIJ;

			u = Vector( cos( beta ), sin( beta ), 0. );

			// indicate which curve last passed through f
            visited.insert( f->index );

			// move to next triangle
			f = g;
		}

		return{ faces.end(), Vector() };
	}
    
    void Mesh :: topologyChange( void )
    {
        indexElements();
        
        buildTreeCotreeDecomposition();
        generatorIndices.resize( nGenerators(), 0. );
        
        if( connection != nullptr )
        {
            delete connection;
            connection = nullptr;
        }
    }
    
    void Mesh :: geometryChange( void )
    {
        for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
        {
            e->updateStar();
        }
        
        if( connection != nullptr )
        {
            delete connection;
            connection = nullptr;
        }
    }
    
    void Mesh :: computeFrameAngles( double initialAngle )
    {
        if( transportRoot->constraintAngle >= 0. )
            transportRoot->alpha = transportRoot->constraintAngle;
        else
            transportRoot->alpha = initialAngle;
        
        for( vector<TransportData>::const_iterator td  = transportOrder.begin();
            td != transportOrder.end();
            td ++ )
        {
            // alphaJ = alphaI + delta - sign*omega
            *(td->alphaJ) = *(td->alphaI) + td->delta - td->sign*(*(td->omega));
            
            // here we're transporting the angle "alphaI" in triangle i across
            // a shared edge to get the angle "alphaJ" in neighboring triangle j;
            // "delta" compensates for the difference between the local
            // reference frames, "omega" is the angle of the connection, and
            // "sign" accounts for the fact that the direction of transport
            // might not be consistent with the orientation of the shared
            // dual edge.
        }
    }
    
    bool Mesh :: hasBoundary( void ) const
    {
        for( HalfEdgeCIter he = halfedges.begin();
            he != halfedges.end();
            he++ )
        {
            if( he->onBoundary )
            {
                return true;
            }
        }
        return false;
    }
    
    bool inPrimalSpanningTree( HalfEdgeCIter he )
    {
        VertexCIter v = he->from;
        VertexCIter w = he->flip->from;
        
        return v->parent == w || w->parent == v;
    }
    
    bool inDualSpanningTree( HalfEdgeCIter he )
    {
        FaceCIter f = he->face;
        FaceCIter g = he->flip->face;
        
        return f->parent == g || g->parent == f;
    }
    
    void Mesh :: buildPrimalSpanningTree( void )
    {
        VertexIter root = vertices.begin();
        while( root->onBoundary() ) root++;
        
        for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
        {
            v->parent = v;
        }
        
        queue<VertexIter> Q;
        Q.push( root );
        while( !Q.empty() )
        {
            VertexIter v = Q.front(); Q.pop();
            
            HalfEdgeIter he = v->out;
            do
            {
                VertexIter w = he->flip->from;
                
                if( w->parent == w &&
                   w != root &&
                   !w->onBoundary() )
                {
                    w->parent = v;
                    Q.push( w );
                }
                
                he = he->flip->next;
            }
            while( he != v->out );
        }
    }
    
    void Mesh :: buildDualSpanningCoTree( void )
    {
        FaceIter root = faces.begin();
        while( root->isBoundary() ) root++;
        
        for( FaceIter f = faces.begin(); f != faces.end(); f++ )
        {
            f->parent = f;
        }
        
        queue<FaceIter> Q;
        Q.push( root );
        while( !Q.empty() )
        {
            FaceIter f = Q.front(); Q.pop();
            
            HalfEdgeIter he = f->he;
            do
            {
                FaceIter g = he->flip->face;
                
                if( g->parent == g &&
                   g != root &&
                   !inPrimalSpanningTree( he ) &&
                   !g->isBoundary() )
                {
                    g->parent = f;
                    Q.push( g );
                }
                
                he = he->next;
            }
            while( he != f->he );
        }
    }
    
    void Mesh :: buildTreeCotreeDecomposition( void )
    {
        buildPrimalSpanningTree();
        buildDualSpanningCoTree();
    }
    
    HalfEdgeIter sharedHalfEdge( VertexIter& v, VertexIter& w )
    {
        HalfEdgeIter he = v->out;
        do
        {
            if( he->flip->from == w )
            {
                return he;
            }
            
            he = he->flip->next;
        }
        while( he != v->out );
        
        assert( 0 );
        
        return he;
    }
    
    HalfEdgeIter sharedHalfEdge( FaceIter& f, FaceIter& g )
    {
        HalfEdgeIter he = f->he;
        do
        {
            if( he->flip->face == g )
            {
                return he;
            }
            
            he = he->next;
        }
        while( he != f->he );
        
        assert( 0 );
        
        return he;
    }
    
    int Mesh :: nGenerators( void ) const
    {
        int n = 0;
        
        for( EdgeCIter e = edges.begin(); e != edges.end(); e++ )
        {
            if( e->onBoundary() ) continue;
            
            if( !inPrimalSpanningTree( e->he ) &&
               !inDualSpanningTree( e->he ))
            {
                n++;
            }
        }
        
        return n;
    }
    
    bool Mesh :: isDualBoundaryLoop( const Cycle& cycle )
    {
        if( cycle.size() == 0 ) return false;
        
        return cycle[0]->from->onBoundary() ||
        cycle[0]->flip->from->onBoundary();
    }
    
    void Mesh :: appendDualGenerators( vector<Cycle>& cycles )
    {
        for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
        {
            if( e->onBoundary() ) continue;
            
            if( !inPrimalSpanningTree( e->he ) &&
               !inDualSpanningTree( e->he ))
            {
                Cycle g, c1, c2;
                FaceIter f;
                
                g.push_back( e->he );
                
                f = e->he->flip->face;
                while( f != f->parent )
                {
                    c1.push_back( sharedHalfEdge( f, f->parent ));
                    f = f->parent;
                }
                
                f = e->he->face;
                while( f != f->parent )
                {
                    c2.push_back( sharedHalfEdge( f, f->parent ));
                    f = f->parent;
                }
                
                int m = c1.size()-1;
                int n = c2.size()-1;
                while( c1[m] == c2[n] ) { m--; n--; }
                for( int i = 0; i <= m; i++ ) g.push_back( c1[i] );
                for( int i = n; i >= 0; i-- ) g.push_back( c2[i]->flip );
                
                // make sure that boundary loops wind around the boundary in a consistent direction
                if( isDualBoundaryLoop( g ))
                {
                    if( g[0]->next->from->onBoundary() )
                    {
                        unsigned int n = g.size();
                        for( unsigned int i = 0; i < n; i++ )
                        {
                            g[i] = g[i]->flip;
                        }
                        
                        for( unsigned int i = 0; i < n/2; i++ )
                        {
                            swap( g[i], g[n-1-i] );
                        }
                    }
                }
                
                cycles.push_back( g );
            }
        }
    }
    
    void Mesh :: appendDirectionalConstraints( vector<Cycle>& cycles, vector<double>& holonomies )
    {
        // first point all faces to themselves to indicate that they have not yet
        // been added to the tree; meanwhile look for a constrained face to serve
        // as the root for our constraint tree (if there aren't any constrained
        // faces, an arbitrary face will work just fine)
        transportRoot = faces.begin();
        for( FaceIter f = faces.begin(); f != faces.end(); f++ )
        {
            f->cParent = f;
            if( f->constraintAngle >= 0. )
            {
                transportRoot = f;
            }
        }
        
        // do a breadth first search on the dual edges and cache the traversal
        // order and angles (we'll need this tree later to construct a global
        // frame starting with a known direction)
        queue<FaceIter> Q;
        Q.push( transportRoot );
        while( !Q.empty() )
        {
            FaceIter f = Q.front(); Q.pop();
            
            // visit neighboring faces
            HalfEdgeIter he = f->he;
            do
            {
                FaceIter g = he->flip->face;
                if( g->cParent == g &&
                   g != transportRoot &&
                   !g->isBoundary() )
                {
                    // point the current neighbor to its parent in the
                    // traversal and enqueue it
                    g->cParent = f;
                    Q.push( g );
                    
                    // also, cache transport via Levi-Civita across this
                    // edge (can then transport via the trivial connection
                    // by simply subtracting connection angles once we
                    // compute them)
                    TransportData td;
                    td.delta = parallelTransport( 0., he );
                    td.sign = he->from->index > he->flip->from->index ? -1. : 1.;
                    td.omega = &(he->edge->theta);
                    td.alphaI = &(he->face->alpha);
                    td.alphaJ = &(he->flip->face->alpha);
                    transportOrder.push_back( td );
                    
                    // check if this face is constrained; if so, we need to add a constraint
                    if( g->constraintAngle >= 0. )
                    {
                        Cycle constraint;
                        double alpha = g->constraintAngle;
                        
                        // follow this face back up the tree to the most
                        // recent constrained ancestor, adding all the halfedges
                        // in between to a new constraint; also compute the difference
                        // between the two constrained directions relative to transport
                        // via the Levi-Civita connection
                        FaceIter h = g;
                        do
                        {
                            HalfEdgeIter he = sharedHalfEdge( h, h->cParent );
                            alpha = parallelTransport( alpha, he );
                            constraint.push_back( he );
                            h = h->cParent;
                        } while( h->constraintAngle < 0. );
                        
                        // add this new constraint to the set of all constraints
                        cycles.push_back( constraint );
                        double gamma = h->constraintAngle;
                        Vector u1( cos(gamma), sin(gamma), 0. );
                        Vector u2( cos(alpha), sin(alpha), 0. );
                        holonomies.push_back( acos(dot(u1,u2)) );
                    }
                }
                he = he->next;
            } while( he != f->he );
        }
    }
    
    double tipAngle( const Vector& x, const Vector& a, const Vector& b )
    // returns the angle between (a-x) and (b-x)
    {
        Vector u = ( a - x ).unit();
        Vector v = ( b - x ).unit();
        
        return atan2( cross(u,v).norm(), dot(u,v) );
    }
    
    double Mesh :: boundaryLoopCurvature( const Cycle& cycle )
    {
        double totalK = 0.;
        
        // get a halfedge of the "virtual" face bounded by the current cycle
        VertexCIter v0 = cycle[0]->flip->next->from;
		HalfEdgeIter he0 = v0->out;
		{
			// try the 'from' vertex of the cycle
			do
			{
				he0 = he0->flip->next;
			} while (!he0->onBoundary && he0 != v0->out);

			// try the 'next' vertex of the cycle
			if (!he0->onBoundary)
			{
				v0 = cycle[0]->flip->from;
				he0 = v0->out;
				do
				{
					he0 = he0->flip->next;
				} while (!he0->onBoundary && he0 != v0->out);
			}
		}
        
        
        // compute a "virtual" vertex in the middle of this loop
        Vector c( 0., 0., 0. );
        HalfEdgeCIter he = he0;
        int boundaryLength = 0;
        do
        {
            c += he->from->position;
            boundaryLength++;
            he = he->next;
        }
        while( he != he0 );
        c /= (double) boundaryLength;
        
        // compute the curvature around the center vertex
        double K = 2.*M_PI;
        he = he0;
        do
        {
            Vector a = he->from->position;
            Vector b = he->next->from->position;
            K -= tipAngle( c, a, b );
            he = he->next;
        }
        while( he != he0 );
        totalK += K;
        
        // add the curvature around each of the boundary vertices, using
        // the following labels:
        //    c - virtual center vertex of boundary loop (computed above)
        //    d - current boundary vertex (we walk around the 1-ring of this vertex)
        //    a,b - consecutive interior vertices in 1-ring of d
        //    e,f - boundary vertices adjacent to d
        he = he0;
        do
        {
            VertexCIter v = he->from;
            Vector d = v->position;
            
            K = 2.*M_PI;
            
            HalfEdgeCIter he2 = v->out;
            do
            {
                if( he2->onBoundary )
                {
                    Vector f = he2->next->from->position;
                    K -= tipAngle( d, f, c );
                }
                else
                {
                    Vector a = he2->next->from->position;
                    Vector b = he2->next->next->from->position;
                    K -= tipAngle( d, a, b );
                    
                    if( he2->flip->onBoundary )
                    {
                        Vector e = he2->flip->from->position;
                        K -= tipAngle( d, c, e );
                    }
                }
                
                he2 = he2->flip->next;
            }
            while( he2 != v->out );
            
            totalK += K;
            
            he = he->next;
        }
        while( he != he0 );
        
        return totalK;
    }
    
    void Mesh :: computeTrivialConnection( void )
    {
        if( connection == nullptr )
        {
            connection = new Connection( *this );
        }
        
        connection->update();
        computeFrameAngles( fieldAngle );
    }
    
    string getFilenameExtension( const string& filename )
    {
        int lastDot = filename.find_last_of( '.' );
        int extLength = filename.length() - lastDot - 1;
        
        string extension = filename.substr( lastDot+1, extLength );
        
        transform( extension.begin(), extension.end(), extension.begin(), ::tolower );
        
        return extension;
    }
    
    void Mesh :: read( const string& filename )
    {
        // make sure we can open the specified file
        std::ifstream in( filename.c_str() );
        if( !in.is_open() )
        {
            cerr << "Error: couldn't open file " << filename << " for input." << endl;
            exit( 1 );
        }
        
        // select format according to filename extension
        string extension = getFilenameExtension( filename );
        if( extension == "obj" )
        {
            MeshIO::readOBJ( in, *this );
        }
        else
        {
            cerr << "Error: input file format " << extension << " not supported!" << endl;
            exit( 1 );
        }
        
        // apply any updates needed after a change to topology or geometry
        topologyChange();
        geometryChange();
    }
    
    void Mesh :: write( const string& filename )
    {
        // make sure we can open the specified file
        ofstream out( filename.c_str() );
        if( !out.is_open() )
        {
            cerr << "Error: couldn't open the file " << filename << " for output." << endl;
            exit( 1 );
        }
        
        // select format according to filename extension
        string extension = getFilenameExtension( filename );
        if( extension == "obj" )
        {
            MeshIO::writeOBJ( out, *this );
        }
        else if( extension == "eobj" )
        {
            MeshIO::writeEOBJ( out, *this );
        }
        else if( extension == "objx" )
        {
            MeshIO::writeOBJX( out, *this );
        }
        else if( extension == "jvx" )
        {
            vector<int> singularityIndices;
            for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
            {
                if( v->k != 0. )
                {
                    singularityIndices.push_back( v->index );
                }
            }
            
            MeshIO::writeJVX( out, *this, singularityIndices );
        }
        else
        {
            cerr << "Error: output file format " << extension << " not supported!" << endl;
            exit( 1 );
        }
    }
    
    int Mesh :: eulerCharacteristic( void ) const
    {
        return vertices.size() - edges.size() + faces.size();
    }
}

