#include "LifeBrush.h"

#include "CMWrapper.h"
#include <algorithm>
#include <cassert>

using namespace std;

namespace cm
{
   Common :: Common( void )
   {
      cholmod_l_start( &common );
   }

   Common :: ~Common( void )
   {
      cholmod_l_finish( &common );
   }

   Common :: operator cholmod_common*( void )
   {
      return &common;
   }

   Sparse :: Sparse( Common& _common, int _m, int _n, int _xtype )
   : common( _common ),
     m( _m ),
     n( _n ),
     xtype( _xtype ),
     A( nullptr )
   {}

   Sparse::Sparse(const Sparse& B) : A(nullptr), common(B.common)
   {
	   *this = B;
   }

   const Sparse& Sparse::operator=(const Sparse& other)
	   // copies B
   {
	   if (A)
	   {
		   cholmod_l_free_sparse(&A, common);
		   A = nullptr;
	   }

	   m = other.m;
	   n = other.n;
	   data = other.data;

	   return *this;
   }

   Sparse :: ~Sparse( void )
   {
	   if (A)
	   {
		   cholmod_l_free_sparse(&A, common);
		   A = nullptr;
	   }
   }

   int Sparse :: size( int dim ) const
   {
      if( dim == 1 ) return m;
      if( dim == 2 ) return n;
      return 0;
   }

   int Sparse :: length( void ) const
   {
      return max( m, n );
   }

   void Sparse :: zero( double rVal, double iVal )
   {
      EntryValue val( rVal, iVal );

      for( EntryMap::iterator i = data.begin(); i != data.end(); i++ )
      {
         i->second = val;
      }
   }

   void Sparse :: transpose( void )
   {
      EntryMap transposed;

      for( const_iterator e = data.begin(); e != data.end(); e++ )
      {
         int i = e->first.first;
         int j = e->first.second;
         transposed[ EntryIndex( j, i ) ] = e->second;
      }

      data = transposed;
      swap( m, n );
   }

   void Sparse :: horzcat( const Sparse& A, const Sparse& B )
   {
      assert( A.m == B.m );
      assert( A.xtype == B.xtype );

      data.clear();
      m = A.m;
      n = A.n + B.n;
      xtype = A.xtype;

      for( const_iterator e = A.begin(); e != A.end(); e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         data[ EntryIndex( j, i ) ] = e->second;
      }

      for( const_iterator e = B.begin(); e != B.end(); e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         data[ EntryIndex( j+A.n, i ) ] = e->second;
      }
   }

   void Sparse :: vertcat( const Sparse& A, const Sparse& B )
   {
      assert( A.n == B.n );
      assert( A.xtype == B.xtype );

      data.clear();
      m = A.m + B.m;
      n = A.n;
      xtype = A.xtype;

      for( const_iterator e = A.begin(); e != A.end(); e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         data[ EntryIndex( j, i ) ] = e->second;
      }

      for( const_iterator e = B.begin(); e != B.end(); e++ )
      {
         int i = e->first.second;
         int j = e->first.first;
         data[ EntryIndex( j, i+A.m ) ] = e->second;
      }
   }

   cholmod_sparse* Sparse :: operator*( void )
   {
      if( A ) { cholmod_l_free_sparse( &A, common ); A = nullptr; }

      int nzmax = data.size();
      int sorted = true;
      int packed = true;
      int stype = 0;
      A = cholmod_l_allocate_sparse( m, n, nzmax, sorted, packed, stype, xtype, common );
      //A->i = (void*) new SuiteSparse_long[ nzmax ];
      //A->x = (void*) new double[ nzmax ];
      //if( xtype == CHOLMOD_ZOMPLEX )
      //{
      //   A->z = (void*) new double[ nzmax ];
      //}

      // build compressed matrix (note that EntryMap stores entries in column-major order)
      double* pr = (double*) A->x;
      double* pi = (double*) A->z;
      SuiteSparse_long* ir = (SuiteSparse_long*) A->i;
      SuiteSparse_long* jc = (SuiteSparse_long*) A->p;
      int i = 0;
      int j = -1;
      for( EntryMap::const_iterator e = data.begin(); e != data.end(); e++ )
      {
         int c = e->first.first;
         if( c != j )
         {
            for( int k = j+1; k <= c; k++ )
            {
               jc[k] = i;
            }
            j = c;
         }

         ir[i] = e->first.second;
         pr[i] = e->second.first;
         if( xtype == CHOLMOD_ZOMPLEX )
         {
            pi[i] = e->second.second;
         }
         i++;
      }
      for( int k = j+1; k <= n; k++ )
      {
         jc[k] = i;
      }

      return A;
   }

   double& Sparse :: operator()( int row, int col )       { return retrieveEntry( row, col, CHOLMOD_REAL    ); }
   double  Sparse :: operator()( int row, int col ) const { return retrieveEntry( row, col, CHOLMOD_REAL    ); }
   double& Sparse ::          r( int row, int col )       { return retrieveEntry( row, col, CHOLMOD_REAL    ); }
   double  Sparse ::          r( int row, int col ) const { return retrieveEntry( row, col, CHOLMOD_REAL    ); }
   double& Sparse ::          i( int row, int col )       { return retrieveEntry( row, col, CHOLMOD_ZOMPLEX ); }
   double  Sparse ::          i( int row, int col ) const { return retrieveEntry( row, col, CHOLMOD_ZOMPLEX ); }

   Sparse::iterator Sparse :: begin( void )
   {
      return data.begin();
   }

   Sparse::const_iterator Sparse :: begin( void ) const
   {
      return data.begin();
   }

   Sparse::iterator Sparse :: end( void )
   {
      return data.end();
   }

   Sparse::const_iterator Sparse :: end( void ) const
   {
      return data.end();
   }

   double& Sparse :: retrieveEntry( int row, int col, int c )
   {
      EntryIndex index( col, row );

      EntryMap::const_iterator entry = data.find( index );
      if( entry == data.end())
      {
         data[ index ] = EntryValue( 0., 0. );
      }

      if( c == CHOLMOD_REAL ) return data[ index ].first;
      else             return data[ index ].second;
   }

   double Sparse :: retrieveEntry( int row, int col, int c ) const
   {
      EntryIndex index( col, row );

      EntryMap::const_iterator entry = data.find( index );

      if( entry == data.end())
      {
         return 0;
      }

      if( c == CHOLMOD_REAL ) return entry->second.first;
      else                    return entry->second.second;
   }

   Dense :: Dense( Common& _common, int _m, int _n, int _xtype )
   : common( _common ),
     m( _m ),
     n( _n ),
     xtype( _xtype ),
     data( nullptr ),
     rData( nullptr ),
     iData( nullptr )
   {
      int d = m; // leading dimension
      data = cholmod_l_allocate_dense( m, n, d, xtype, common );
      rData = (double*) data->x;
      if( xtype == CHOLMOD_ZOMPLEX )
      {
         iData = (double*) data->z;
      }
   }

   Dense :: Dense( const Dense& A )
   : common( A.common )
   {
      *this = A;
   }

   Dense :: ~Dense( void )
   {
	   if (data)
	   {
		   cholmod_l_free_dense(&data, common);
		   data = nullptr;
	   }
   }

   cholmod_dense* Dense :: operator*( void )
   {
      return data;
   }

   void Dense :: initializeFromCopy( void )
   {
      assert( data->dtype == CHOLMOD_DOUBLE );

      xtype = data->xtype;
      m = data->nrow;
      n = data->ncol;
      rData = (double*) data->x;
      if( data->xtype == CHOLMOD_ZOMPLEX )
      {
         iData = (double*) data->z;
      }
   }

   const Dense& Dense :: operator=( const Dense& A )
   {
      if( data ) { cholmod_l_free_dense( &data, common ); data = nullptr; }

      data = cholmod_l_copy_dense( A.data, common );

      initializeFromCopy();

      return *this;
   }

   const Dense& Dense :: operator=( cholmod_dense* A )
   {
      if( data ) { cholmod_l_free_dense( &data, common ); data = nullptr; }

      data = A;

      initializeFromCopy();

      return *this;
   }

   int Dense :: size( int dim ) const
   {
      if( dim == 1 ) return m;
      if( dim == 2 ) return n;
      return 0;
   }

   int Dense :: length( void ) const
   {
      return max( m, n );
   }

   void Dense :: zero( double rVal, double iVal )
   {
      if( rData )
      for( int i = 0; i < m*n; i++ )
      {
         rData[i] = rVal;
      }

      if( iData )
      for( int i = 0; i < m*n; i++ )
      {
         iData[i] = iVal;
      }
   }

   void Dense :: horzcat( const Dense& A, const Dense& B )
   {
      assert( A.m == B.m );
      assert( A.xtype == B.xtype );

      m = A.m;
      n = A.n + B.n;
      xtype = A.xtype;

      if( data ) { cholmod_l_free_dense( &data, common ); data = nullptr; }
      data = cholmod_l_allocate_dense( m, n, m, xtype, common );
      rData = (double*) data->x;
      if( xtype == CHOLMOD_ZOMPLEX )
      {
         iData = (double*) data->z;
      }

      for( int i = 0; i < A.m; i++ )
      for( int j = 0; j < A.n; j++ )
      {
         rData[i+m*j] = A.rData[i+A.m*j];
         if( xtype == CHOLMOD_ZOMPLEX )
         {
            iData[i+m*j] = A.iData[i+A.m*j];
         }
      }

      for( int i = 0; i < B.m; i++ )
      for( int j = 0; j < B.n; j++ )
      {
         rData[i+m*(j+A.n)] = B.rData[i+B.m*j];
         if( xtype == CHOLMOD_ZOMPLEX )
         {
            iData[i+m*(j+A.n)] = B.iData[i+B.m*j];
         }
      }
   }

   void Dense :: vertcat( const Dense& A, const Dense& B )
   {
      assert( A.n == B.n );
      assert( A.xtype == B.xtype );

      m = A.m + B.m;
      n = A.n;
      xtype = A.xtype;

      if( data ) { cholmod_l_free_dense( &data, common ); data = nullptr; }
      data = cholmod_l_allocate_dense( m, n, m, xtype, common );
      rData = (double*) data->x;
      if( xtype == CHOLMOD_ZOMPLEX )
      {
         iData = (double*) data->z;
      }

      for( int i = 0; i < A.m; i++ )
      for( int j = 0; j < A.n; j++ )
      {
         rData[i+m*j] = A.rData[i+A.m*j];
         if( xtype == CHOLMOD_ZOMPLEX )
         {
            iData[i+m*j] = A.iData[i+A.m*j];
         }
      }

      for( int i = 0; i < B.m; i++ )
      for( int j = 0; j < B.n; j++ )
      {
         rData[(i+A.m)+m*j] = B.rData[i+B.m*j];
         if( xtype == CHOLMOD_ZOMPLEX )
         {
            iData[(i+A.m)+m*j] = B.iData[i+B.m*j];
         }
      }
   }

   double& Dense :: operator()( int row, int col )       { return rData[ row + m*col ]; }
   double  Dense :: operator()( int row, int col ) const { return rData[ row + m*col ]; }
   double& Dense ::          r( int row, int col )       { return rData[ row + m*col ]; }
   double  Dense ::          r( int row, int col ) const { return rData[ row + m*col ]; }
   double& Dense ::          i( int row, int col )       { return iData[ row + m*col ]; }
   double  Dense ::          i( int row, int col ) const { return iData[ row + m*col ]; }
   double& Dense :: operator()( int index )              { return rData[ index ]; }
   double  Dense :: operator()( int index ) const        { return rData[ index ]; }
   double& Dense ::          r( int index )              { return rData[ index ]; }
   double  Dense ::          r( int index ) const        { return rData[ index ]; }
   double& Dense ::          i( int index )              { return iData[ index ]; }
   double  Dense ::          i( int index ) const        { return iData[ index ]; }
}

