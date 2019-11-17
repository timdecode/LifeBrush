//============================================================
// CMWrapper.h
// 
// Wrapper classes for CHOLMOD sparse matrices and CHOLMOD
// common objects.  Many features of CHOLMOD are not
// supported -- this wrapper has mainly been used for
// solving real linear systems via SPQR.
//

#ifndef CMWRAPPER_H
#define CMWRAPPER_H

#include <cholmod.h>
#include <map>

namespace cm
{
   class Common
   {
      public:
         Common( void );
         ~Common( void );

         operator cholmod_common*( void );
         // allows cm::Common to be treated as a cholmod_common*

      protected:
         cholmod_common common;
   };

   class Dense
   {
      public:
         Dense( Common& common, int m = 0, int n = 0, int xtype = CHOLMOD_REAL );
         // initialize an mxn matrix of doubles
         // xtype is either CHOLMOD_REAL or CHOLMOD_ZOMPLEX
         
         Dense( const Dense& A );
         // copy constructor
         
         ~Dense( void );
         // destructor

         cholmod_dense* operator*( void );
         // dereference operator gets pointer to underlying cholmod_dense data structure

         const Dense& operator=( const Dense& A );
         // copies A

         const Dense& operator=( cholmod_dense* A );
         // gets pointer to A; will deallocate A upon destruction

         int size( int dim ) const;
         // returns the size of the dimension specified by scalar dim

         int length( void ) const;
         // returns the size of the largest dimension

         void zero( double rVal = 0., double iVal = 0. );
         // sets all elements to rVal+iVal*i

         void horzcat( const Dense& A, const Dense& B );
         // replaces the current matrix with [ A, B ]

         void vertcat( const Dense& A, const Dense& B );
         // replaces the current matrix with [ A; B ]

         double& operator()( int row, int col );
         double  operator()( int row, int col ) const;
         double& r( int row, int col );
         double  r( int row, int col ) const;
         // access real part of element (row,col)
         // note: uses 0-based indexing

         double& i( int row, int col );
         double  i( int row, int col ) const;
         // access imaginary part of element (row,col)
         // note: uses 0-based indexing
         
         double& operator()( int index );
         double  operator()( int index ) const;
         // access real part of element ind of a vector
         // note: uses 0-based indexing
         
         double& r( int index );
         double  r( int index ) const;
         // access real part of element ind of a vector
         // note: uses 0-based indexing
         
         double& i( int index );
         double  i( int index ) const;
         // access imaginary part of element ind of a vector
         // note: uses 0-based indexing

      protected:
         void initializeFromCopy( void );

         Common& common;
         int m, n;
         int xtype;
         cholmod_dense* data;
         double* rData;
         double* iData;
   };

   class Sparse
   {
      public:
         Sparse( Common& common, int m = 0, int n = 0, int xtype = CHOLMOD_REAL );
         // initialize an mxn matrix of doubles
         // xtype is either CHOLMOD_REAL or CHOLMOD_ZOMPLEX

		 Sparse(const Sparse& B);
			 // copy constructor

         ~Sparse( void );
         
         cholmod_sparse* operator*( void );
         // dereference operator gets pointer to underlying cholmod_sparse data structure

         int size( int dim ) const;
         // returns the size of the dimension specified by scalar dim

         int length( void ) const;
         // returns the size of the largest dimension

         void zero( double rVal = 0., double iVal = 0. );
         // sets all nonzero elements to rVal+iVal*i
         
         void transpose( void );
         // replaces this matrix with its transpose

         void horzcat( const Sparse& A, const Sparse& B );
         // replaces the current matrix with [ A, B ]

         void vertcat( const Sparse& A, const Sparse& B );
         // replaces the current matrix with [ A; B ]

         double& operator()( int row, int col );
         double  operator()( int row, int col ) const;
         double& r( int row, int col );
         double  r( int row, int col ) const;
         // access real part of element (row,col)
         // note: uses 0-based indexing

         double& i( int row, int col );
         double  i( int row, int col ) const;
         // access imaginary part of element (row,col)
         // note: uses 0-based indexing
         
         typedef std::pair<int,int> EntryIndex; // NOTE: column THEN row! (makes it easier to build compressed format)
         typedef std::pair<double,double> EntryValue;
         typedef std::map<EntryIndex,EntryValue> EntryMap;
         typedef EntryMap::iterator       iterator;
         typedef EntryMap::const_iterator const_iterator;

               iterator begin( void );
         const_iterator begin( void ) const;
               iterator   end( void );
         const_iterator   end( void ) const;

		 const Sparse& operator=(const Sparse& other); // copies B;

      protected:
         Common& common;
         int m, n;
         int xtype;
         cholmod_sparse* A = nullptr;
         EntryMap data;

         double& retrieveEntry( int row, int col, int c );
         double  retrieveEntry( int row, int col, int c ) const;
   };
}

#endif
