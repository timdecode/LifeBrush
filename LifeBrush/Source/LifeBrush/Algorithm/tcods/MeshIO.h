//============================================================
// MeshIO.h
//
// MeshIO contains a collection of static members for polygon
// mesh input/output via the Mesh class.
//

#ifndef TCODS_MESHIO_H
#define TCODS_MESHIO_H

#include <iosfwd>
#include <string>
#include <sstream>
#include <string>
#include <vector>

#include "Vector.h"

namespace tcods
{
    class Mesh;

    class MeshIO
    {
    public:
        
        class Index
        {
        public:
            Index( void )
            {}
            
            Index( int p, int t, int n )
            : position( p ), texcoord( t ), normal( n )
            {}
            
            int position;
            int texcoord;
            int normal;
        };
        
        class MeshData
        {
        public:
            std::vector<DDG::Vector> positions;
            std::vector<DDG::Vector> texcoords;
            std::vector<DDG::Vector> normals;
            std::vector< std::vector< Index > > indices;
        };

        static int buildMesh( const MeshData& data, Mesh& mesh );
		static void buildMesh2(const MeshData& data, Mesh& mesh);

		static bool hasIsolatedVertices(const Mesh& mesh);
		static bool hasNonManifoldVertices(const Mesh& mesh);

		// routines for reading and writing various mesh formats
        // assumes valid, open streams that point to the target file
        
        // read
        static void readOBJ( std::istream& in, Mesh& mesh ); // reads WavefrontOBJ
        
        // write
        static void  writeOBJ( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ
        static void writeOBJX( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ plus tangent vector per vertex
        static void writeEOBJ( std::ostream& out, const Mesh& mesh ); // writes WavefrontOBJ plus tangent vector per face
        
        static void  writeJVX( std::ostream& out, const Mesh& mesh, const std::vector<int>& singularityIndices );
        // writes XML containing mesh data, plus pair of orthogonal tangent vector fields, singularity
        // locations, and "matchings" between neighboring triangles which define a quad-covering of the surface
        // (can be used as input for QuadCover)
        
        
    protected:
        static void readMeshData( std::istream& in, MeshData& data );
        static void readPosition( std::stringstream& ss, MeshData& data );
        static void readTexCoord( std::stringstream& ss, MeshData& data );
        static void readNormal  ( std::stringstream& ss, MeshData& data );
        static void readFace    ( std::stringstream& ss, MeshData& data );
        static Index parseFaceIndex( const std::string& token );
        static void preallocateMeshElements( const MeshData& data, Mesh& mesh );
		static void sizeToFit(Mesh& mesh, Mesh& newMesh);
	};
}

#endif

