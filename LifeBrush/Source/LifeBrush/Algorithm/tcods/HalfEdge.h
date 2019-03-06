//============================================================
// HalfEdge.h
// 
// Mesh represents a polygonal mesh via a half edge data structure.
//
// Texture coordinates and vertex normals are stored on half edges
// to allow different values at each triangle corner.
//
// Boundaries are handled by creating "virtual" faces for each boundary
// cycle and marking all HalfEdges incident on these faces.
//

#ifndef TCODS_HALFEDGE_H
#define TCODS_HALFEDGE_H

#include "Vector.h"
#include <iosfwd>
#include <vector>
#include <set>
#include <string>

namespace tcods
{
   // forward declarations
   class HalfEdge;
   class Vertex;
   class Edge;
   class Face;
   class Mesh;
   class MeshIO;
   class Connection;


   // element iterators -- *CIter and *Iter are for const- and non-const access, respectively
   // most routines will want to work with element iterators rather than element references
   typedef std::vector<HalfEdge> :: iterator HalfEdgeIter; typedef std::vector<HalfEdge> :: const_iterator HalfEdgeCIter;
   typedef std::vector<Vertex>   :: iterator   VertexIter; typedef std::vector<Vertex>   :: const_iterator   VertexCIter;
   typedef std::vector<Edge>     :: iterator     EdgeIter; typedef std::vector<Edge>     :: const_iterator     EdgeCIter;
   typedef std::vector<Face>     :: iterator     FaceIter; typedef std::vector<Face>     :: const_iterator     FaceCIter;

   class HalfEdge
   {
      public:
         HalfEdgeIter next; // next halfedge in face
         HalfEdgeIter flip; // opposite halfedge across edge
         VertexIter from;   // originating vertex
         EdgeIter edge;     // incident edge
         FaceIter face;     // incident face

         int index;       // unique ID for this HalfEdge (0-based)
         bool onBoundary; // true if this halfedge is contained in a boundary face
         DDG::Vector normal;   // vertex normal for "from" vertex
		 DDG::Vector texcoord; // texture coordinates for "from" vertex
   };

   class Vertex
   {
      public:
         int valence( void ) const;     // returns number of incident edges
         double defect( void ) const;   // returns discrete Gaussian curvature given by 2 pi minus sum of tip angle defects
         bool onBoundary( void ) const; // returns true iff this vertex is contained in the boundary
		 bool isIsolated(void) const;
         HalfEdgeIter out; // one of the outgoing halfedges

         int index;         // unique ID for this Vertex (0-based)
		 DDG::Vector position;   // coordinates in R^3
         double k;          // singularity index (equal to target holonomy divided by 2 pi)
         VertexIter parent; // parent in tree-cotree decomposition
   };

   class Edge
   {
      public:
         bool onBoundary( void ) const; // returns true iff this edge is contained in the boundary
         void updateStar( void );       // recomputes value of Hodge star for this edge

         HalfEdgeIter he; // one of the halfedges of this edge

         int index;    // unique ID for this Edge (0-based)
         double theta; // connection between incident faces
         double star;  // Hodge star on primal 1-forms
   };

   class Face
   {
      public:
         double area( void ) const;                  // triangle area
		 DDG::Vector normal( void ) const;                // returns unit normal with orientation relative to halfedge circulation
         double circumradius( void ) const;          // returns circumradius
		 DDG::Vector barycenter( void ) const;            // returns mean of vertex positions
         void frame(DDG::Vector& e1, DDG::Vector& e2 ) const; // computes a canoncal basis for the tangent plane
		 DDG::Vector toLocal(DDG::Vector p ) const;           // maps point in R^3 to canonical basis
		 DDG::Vector toGlobal(DDG::Vector q ) const;          // maps point in canonical basis to R^3
         bool isBoundary( void ) const;              // returns true iff this face is a "virtual" face
		 DDG::Vector smoothNormal(DDG::Vector p );
		 bool contains(DDG::Vector p ) const;
         HalfEdgeIter he; // one of the halfedges of this face

         int index;              // unique ID for this Face (0-based)
         double alpha;           // smooth frame angle
         FaceIter parent;        // parent in tree-cotree decomposition
         int curveIndex;         // index of last curve to pass through this face
         FaceIter cParent;       // parent in constraint tree
         double constraintAngle = -1.; // negative implies unconstrained
   };

   class Mesh
   {
      friend class MeshIO;

      public:
         typedef std::vector<HalfEdgeIter> Cycle;
         Mesh( void ) {}
         Mesh( const Mesh& mesh );                                // copy constructor
         virtual ~Mesh( void );                                   // destructor
         const Mesh& operator=( const Mesh& mesh );               // assignment operator
         void topologyChange( void );                             // must be called by any method that modifies mesh connectivity
         void geometryChange( void );                             // must be called by any method that modifies mesh geometry
		 void indexElements(void);								  // assigns a unique ID to each mesh element
		 int nFaces( void ) const;                                // number of faces, excluding virtual faces
         int nGenerators( void ) const;                           // returns the number of independent noncontractible loops
         bool hasBoundary( void ) const;                          // returns true iff the mesh has boundary
         void computeFrameAngles( double initialAngle );          // computes a global direction field using the trivial connection
         void appendDualGenerators( std::vector<Cycle>& cycles ); // appends generators on the dual graph to "cycles"
         static bool isDualBoundaryLoop( const Cycle& cycle );    // returns true only if cycle is along boundary
         double boundaryLoopCurvature( const Cycle& cycle );      // computes the Riemannian holonomy around a boundary loop
         void computeTrivialConnection( void );                   // compute trivial connection with specified singularities
         void  read( const std::string& filename );               // write to disk; filetype is inferred from path
         void write( const std::string& filename );               // read from disk; filetype is inferred from path
         int eulerCharacteristic( void ) const;                   // returns chi = V-E+F (including virtual faces)
         void appendDirectionalConstraints( std::vector<Cycle>& cycles, std::vector<double>& holonomies );

         static double parallelTransport( double phi, HalfEdgeCIter he ); // transport a direction across an edge using Levi-Civita
		 double parallelTransport_2( double phi, HalfEdgeCIter he ); // Tim: version of parallelTransport that doesn't flip vectors based on index
         static double defect( const Cycle& c ); // computes angle defect resulting from parallel transport around a dual cycle c
         
         void integralCurve( FaceIter initialFace,
                             const DDG::Vector& initialPoint,
                             double initialAngle,
                             std::vector<DDG::Vector>& curve,
                             std::vector<DDG::Vector>& normals,
                             int maxPts = 3000 );


         std::pair<FaceIter, DDG::Vector>
         integralWalk( FaceIter initialFace,
                       const DDG::Vector& initialPoint,
                       double initialAngle,
                       int maxPts,
                       double distanceToTravel );

         // mesh element storage
         std::vector<HalfEdge> halfedges;
         std::vector<Vertex>   vertices;
         std::vector<Edge>     edges;
         std::vector<Face>     faces;

         // get mesh element iterators by index
         // (std::vector::operator[] would instead return a reference)
         HalfEdgeIter  halfedge( int index );
         VertexIter      vertex( int index );
         EdgeIter          edge( int index );
         FaceIter          face( int index );
         HalfEdgeCIter halfedge( int index ) const;
         VertexCIter     vertex( int index ) const;
         EdgeCIter         edge( int index ) const;
         FaceCIter         face( int index ) const;

         std::vector<double> generatorIndices; // target holonomy around generators divided by 2pi
         double fieldAngle;                    // initial angle for global direction field

      protected:
         void buildTreeCotreeDecomposition( void ); // updates the tree-cotree decomposition
         void buildPrimalSpanningTree( void );      // builds a spanning tree of primal edges
         void buildDualSpanningCoTree( void );      // builds a spanning tree of dual edges that do not cross the primal tree

         // stores the element iterator associated with each element ID
         std::vector<HalfEdgeIter> index2he;
         std::vector<  VertexIter> index2vertex;
         std::vector<    EdgeIter> index2edge;
         std::vector<    FaceIter> index2face;

         Connection* connection = nullptr; // represents a trivial connection on the mesh

         // transportOrder caches data needed to update direction field from connection
         struct TransportData
         {
            double delta;
            double sign;
            double* omega;
            double* alphaI;
            double* alphaJ;
         };
         std::vector<TransportData> transportOrder;
         FaceIter transportRoot;
   };
}

#endif

