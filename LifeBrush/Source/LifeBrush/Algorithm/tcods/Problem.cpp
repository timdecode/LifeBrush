#include "LifeBrush.h"

#include "Problem.h"
#include "HalfEdge.h"
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

namespace tcods
{
   Problem :: Problem( void )
   {}
   
   Problem :: Problem( istream& in )
   {
      read( in );
   }
   
   void Problem :: solve( void ) const
   {
      Mesh mesh;
   
      mesh.read( inputPath.c_str() );
   
      int V = mesh.vertices.size();
      for( vector<Singularity>::const_iterator v  = singularities.begin();
                                               v != singularities.end();
                                               v ++ )
      {
         int i = v->first;
         double k = v->second;
   
         if( i < 0 || i >= V )
         {
            cerr << "Warning: singularity requested at vertex " << i+1 << " (mesh has only " << V << " vertices!)" << endl;
         }
         else
         {
            mesh.vertex( i )->k = k;
         }
      }
   
      int G = mesh.nGenerators();
      for( vector<Generator>::const_iterator g  = generators.begin();
                                             g != generators.end();
                                             g ++ )
      {
         int i = g->first;
         double k = g->second;
   
         if( i < 0 || i >= G )
         {
            cerr << "Warning: additional holonomy requested around generator " << i+1 << " (mesh has only " << G << " generators!)" << endl;
         }
         else
         {
            mesh.generatorIndices[ i ] = k;
         }
      }
   
      mesh.fieldAngle = fieldAngle;
   
      mesh.computeTrivialConnection();
   
      mesh.write( outputPath.c_str() );
   }
   
   void Problem :: read( istream& in )
   {
      singularities.clear();
      generators.clear();
   
      string s;
      while( getline( in, s ))
      {
         stringstream line( s );
         string token;
   
         line >> token;
         transform( token.begin(), token.end(), token.begin(), ::tolower );
   
         if( token == "in" )
         {
            line >> inputPath;
         }
         else if( token == "out" )
         {
            line >> outputPath;
         }
         else if( token == "vertex" )
         {
            Singularity s;
            line >> s.first >> s.second;
            singularities.push_back( s );
         }
         else if( token == "generator" )
         {
            Generator g;
            line >> g.first >> g.second;
            generators.push_back( g );
         }
         else if( token == "angle" )
         {
            line >> fieldAngle;
         }
      }
   }
}

