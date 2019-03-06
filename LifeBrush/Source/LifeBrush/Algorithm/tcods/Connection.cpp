#include "LifeBrush.h"

#include "LinearContext.h"
#include "Real.h"

#include "Connection.h"

// MS doesn't define M_PI without this... ummm okay
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <iostream>
#include "SparseMatrix.h"

using namespace std;

namespace DDG
{
	extern LinearContext context;
}

namespace tcods
{
	using namespace DDG;

   Connection :: Connection( Mesh& _mesh )
   : mesh( _mesh ),
     K(),
     b(),
     QR( nullptr )
   {
      build();
   }
   
   Connection :: ~Connection( void )
   {
	   if (QR != nullptr)
	   {
		   SuiteSparseQR_free<double>(&QR, DDG::context);
		   QR = nullptr;
	   }
   }
   
   double sgn( double x )
   {
      if( x > 0. ) return  1.;
      if( x < 0. ) return -1.;
      return 0.;
   }
   
   void Connection :: build( void )
   // build the system of constraint equations
   {
      // construct a list of basis cycles, represented by a vector
      // of oriented edges (really: halfedges)
      std::vector<Mesh::Cycle> basisCycles;
   
      // append contractible basis cycles
      append1RingBases( basisCycles );
   
      unsigned int nContractibleCycles = basisCycles.size();
   
      // append noncontractible basis cycles
      mesh.appendDualGenerators( basisCycles );
   
      nBasisCycles = basisCycles.size();
   
      // append derivative constraints used to specify directional constraints in faces
      vector<double> directionalHolonomies;
      mesh.appendDirectionalConstraints( basisCycles, directionalHolonomies );
   
      // build constraint matrix
	  DDG::SparseMatrix<Real> A( mesh.edges.size(), basisCycles.size() );
      buildCycleMatrix( A, basisCycles );
      applyCotanWeights( A );
   
      // prefactor [Q,R,E] = qr( A' )
      if( QR != nullptr) { SuiteSparseQR_free<double>( &QR, DDG::context ); QR = nullptr; }
      QR = SuiteSparseQR_factorize <double> (7, -2., A.to_cholmod(), DDG::context);
      
      // compute Riemannian holonomy of basis cycles
      K = DDG::DenseMatrix<Real>( basisCycles.size(), 1 );
      b = DDG::DenseMatrix<Real>( basisCycles.size(), 1 );
      generatorOnBoundary.resize( mesh.nGenerators() );
      for( unsigned int i = 0; i < nContractibleCycles; i++ )
      {
         K( i ) = -mesh.vertex( i )->defect();
      }

	  size_t nContratibleCycles_plus_nGenerators = nContractibleCycles + mesh.nGenerators();
      for( unsigned int i = nContractibleCycles;
                        i < nContratibleCycles_plus_nGenerators;
                        i++ )
      {
         if( Mesh::isDualBoundaryLoop( basisCycles[i] ))
         {
            K( i ) = -mesh.boundaryLoopCurvature( basisCycles[i] );
            generatorOnBoundary[ i-nContractibleCycles ] = true;
         }
         else
         {
            K( i ) = -mesh.defect( basisCycles[i] );
            generatorOnBoundary[ i-nContractibleCycles ] = false;
         }
      }
      
      // specify change in angle for each directional constraint
      for( unsigned int i = 0; i < directionalHolonomies.size(); i++ )
      {
         K( i+nBasisCycles ) = -directionalHolonomies[i];
      }
      
      // setup the right hand side using the Riemannian holonomy, ignoring
      // singularities for now; also make sure rhsChanged() is initially true
      for( int i = 0; i < b.length(); i++ )
      {
         b( i ) = K( i );
      }
   }
   
   void Connection :: setupRHS( void )
   // add 2*pi*k to the right hand side, where k is the vector of singularity/generator indices
   {
      double indexSum = 0;
   
      // iterate over vertices
      for( VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         int i = vertex2row[ v->index ]; // get the row index of the current vertex
   
		 // Tim: v is on a boundary \see append1RingBases
		 if (i < 0)
			 continue;

         b(i) = K(i) + 2*M_PI*v->k;
   
         indexSum += v->k;
      }
   
      // iterate over generators
      int nContractibleCycles = nBasisCycles - mesh.generatorIndices.size();
      for( int i = 0; i < (int)mesh.generatorIndices.size(); i++ )
      {
         int j = nContractibleCycles + i;
         b(j) = K(j) + 2.*M_PI*mesh.generatorIndices[i];
   
         if( generatorOnBoundary[i] )
         {
            indexSum += mesh.generatorIndices[i];
         }
      }
   
      // display a warning if the sum of singular indices does not
      // add up to the Euler characteristic
      if( abs( indexSum-mesh.eulerCharacteristic() ) > 1e-7 )
      {
         cerr << endl;
         cerr << "  *************************************************************" << endl;
         cerr << "    Warning: indices do not add up to Euler characteristic!" << endl;
         cerr << "             (solution may have unwanted singularities)" << endl;
         cerr << endl;
         cerr << "             Euler characteristic: " << mesh.eulerCharacteristic() << endl;
         cerr << "                   sum of indices: " << indexSum << endl;
         cerr << "  *************************************************************" << endl;
         cerr << endl;
      }
   }
   
   void Connection :: resetRHS( void )
   {
      // make a copy of the right hand side used in the most recent solve
      for( int i = 0; i < b.length(); i++ )
      {
         b(i) = K(i);
      }
   }
   
   bool Connection :: update( void )
   {
	   DenseMatrix<Real> x;
	   DenseMatrix<Real> y;
   
      // specify singular values in right hand side b
      setupRHS();
   
      // solve y = R'\(E'*b)
      y = SuiteSparseQR_solve (SPQR_RTX_EQUALS_ETB, QR, b.to_cholmod(), DDG::context) ;
   
      // compute w = Q*y
      x = SuiteSparseQR_qmult (SPQR_QX, QR, y.to_cholmod(), DDG::context) ;
   
      applyCotanWeights( x ); 
   
      for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++ )
      {
         e->theta = x( e->index );
      }
   
      // restore original right hand side
      resetRHS();
   
      return true;
   }
   
   void Connection :: applyCotanWeights( DDG::SparseMatrix<Real>& A )
   {

	   for (auto& e : A)
	   {
         int row = e.first.second;

         Real& val = e.second;
         double s = mesh.edge( row )->star;
   
         val *= sqrt( s );
      }
   }
   
   void Connection :: applyCotanWeights(DenseMatrix<Real>& x )
   {
      for( EdgeIter e = mesh.edges.begin(); e != mesh.edges.end(); e++ )
      {
         x( e->index ) *= sqrt( e->star );
      }
   }
   
   void Connection :: buildCycleMatrix(DDG::SparseMatrix<Real>& A, vector<Mesh::Cycle>& cycles ) const
   {
      for( unsigned int l = 0; l < cycles.size(); l++ )
      {
         for( Mesh::Cycle::iterator h  = cycles[l].begin();
                                    h != cycles[l].end();
                                    h ++ )
         {
            int k = (*h)->edge->index;
            int i = (*h)->from->index;
            int j = (*h)->flip->from->index;
   
            if( i > j ) A( k, l ) = -1.;
            else        A( k, l ) =  1.;
         }
      }
   }
   
   void Connection :: append1RingBases( vector<Mesh::Cycle>& cycles )
   {
      // contractible bases
      vertex2row.resize( mesh.vertices.size() );
      for( VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++ )
      {
         if( v->onBoundary() )
         {
            vertex2row[ v->index ] = -1;
            continue;
         }
   
         Mesh::Cycle c;
         HalfEdgeIter he = v->out;
         do
         {
            c.push_back( he );
            he = he->flip->next;
         }
         while( he != v->out );
   
         vertex2row[ v->index ] = cycles.size();
         cycles.push_back( c );
      }
   }
}

