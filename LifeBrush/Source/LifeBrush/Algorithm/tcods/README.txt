===============================================================================

                               TCODS v0.1

                              Keenan Crane
                           September 18, 2010

===============================================================================


--------------
0. CONTENTS
--------------

  0. CONTENTS
  1. ABOUT
  2. DEPENDENCIES
  3. BUILDING
  4. TESTING
  5. USAGE
  6. FILE FORMATS
  7. SOURCE CODE
  8. LICENSE


-----------
1. ABOUT
-----------

This archive contains a C++ implementation of an algorithm for computing
direction fields with prescribed singularities on polygon meshes.  It is
based on the paper

  K. Crane, M. Desbrun, and P. Schršder, "Trivial Connections on Discrete
  Surfaces," Computer Graphics Forum (SGP) 2010.

The command-line utility included with this implementation takes a manifold,
triangulated mesh (with or without boundary) plus a specification of
singular points as input, and produces the same mesh plus a single tangent
vector in each face.  The direction field is guaranteed to have exactly the
requested singularities so long as the sum of the requested singular
indices adds up to the Euler characteristic X = |V|-|E|+|F| of the input
surface.

PLEASE NOTE that this is research code and has not been tested extensively
(in fact, not at all!) on other systems.  Luckily, it was written by a
friendly researcher who is happy to help you out via email in your times
of trouble: keenan@cs.caltech.edu


------------------
2. DEPENDENCIES
------------------

TCODS depends on Tim Davis' SuiteSparseQR sparse QR factorization library:

http://www.cise.ufl.edu/research/sparse/SPQR/

which in turn depends on SuiteSparse and METIS:

http://www.cise.ufl.edu/research/sparse/SuiteSparse/
http://glaros.dtc.umn.edu/gkhome/views/metis

as well as some (hopefully optimized!) BLAS/LAPACK implementation.  On
UNIX-like systems you will probably end up needing the libraries

  bamd.a
  libcamd.a
  libcolamd.a
  libccolamd.a
  libcholmod.a
  libspqr.a

from SuiteSparse and

  libmetis.a

from METIS.  If you want to avoid compiling all of SuiteSparse, you can simply
type "make" in each of the appropriate Lib directories (e.g., AMD/Lib) after
setting up UFConfig and copying the resulting library (.a) files to
/usr/local/bin or some other appropriate place.  Further instructions on
building SuiteSparse and its dependencies can be found on the SuiteSparse home
page.

On Mac OS X, the easiest way to link to an efficient BLAS
implementation is by adding the framework

  -framework Accelerate

On other platforms, Kazushige Goto's GotoBLAS library is a popular choice:

  http://www.tacc.utexas.edu/tacc-projects/gotoblas2/


--------------
3. BUILDING
--------------

On UNIX-like systems, you should simply need to type

  make

at the command line in the root install directory.  This should build the
command-line utility "tcods."  On Windows it may be simplest to install Cygwin:

  http://www.cygwin.com/

which includes an implementation of the GNUMake system which can be used to
execute the Makefiles for TCODS, SuiteSparse, and METIS.  If you need to build
code that does not depend on Cygwin DLLs, MINGW is an option:

  http://www.mingw.org/

Finally, it should not be too difficult to setup a VisualStudio project for
TCODS itself, though I have no experience building SuiteSparse, METIS, or BLAS
with VisualStudio.


-------------
4. TESTING
-------------

To test that TCODS built correctly, type

  ./tcods test/problem.txt
  diff test/solution.eobj test/reference_solution.eobj

at the command line.  Depending on output formatting and the libraries used
for BLAS, linear solver, etc., your result may not match the reference
solution exactly, but it should be pretty close.

For additional testing on Mac OS X, you can compare with the results
computed by Comb:

  http://www.cs.caltech.edu/~keenan/project_tcods.html
  

-----------
5. USAGE
-----------

The tcods executable takes as input a file describing a "problem," i.e., an
input mesh file and a specification of singularities.  Input is given by a
standard ASCII text file containing lines of the following format in any
order:

  in [path to input mesh file]
  
  out [path to output data]
  
  vertex [0-based vertex ID] [target holonomy]
  ...
  
  generator [0-based generator ID] [target holonomy]
  ...
  
  angle [initial field angle]

Terms in square brackets [] are specified by the user; vertex and generator IDs
must be within a valid range (namely, nonnegative integers less than the total
number of vertices and non-contractible cycles in the mesh,
respectively).  An example input file can be found in test/problem.txt

Once a problem file has been created, the problem can be solved by typing

  ./tcods problem.txt

at the command line.  On Windows, it may be convenient to create a batch (.bat)
file that you can run from Explorer (an example is provided).


------------------
6. FILE FORMATS
------------------

The current version of TCODS supports several simple, human-readable mesh
formats.  These formats are described in the file FORMATS.txt.  Additional
formats can be added by modifying the MeshIO class and the Mesh::read() method
in HalfEdge.cpp.


-----------------
7. SOURCE CODE
-----------------

As an alternative to the command-line utility, you can of course interface with
the code directly.  I've done my best to clean up and comment, but (as already
noted) the code evolved along with the research, so be prepared for some
funkiness!

Basic usage can be seen in the method Problem::solve() from Problem.cpp.  For
instance, here's how you might put a few arbitrary singularities on a bunny
mesh:

  Mesh mesh;
  mesh.read( "bunny.obj" );          // read a polygon mesh
  mesh.vertex( 0 )->k = 1.;          // specify some singularities
  mesh.vertex( 100 )->k = .5;
  mesh.vertex( 200 )->k = .5;
  mesh.fieldAngle = 1.23;            // specify the field direction
  mesh.computeTrivialConnection();   // compute the solution
  mesh.write( "hairybunny.eobj" );   // write the mesh plus direction field

The code also includes basic support for generating integral curves of the
direction field in Mesh::integralCurve() (see HalfEdge.h).  For instance,
the following code will compute a single integral curve starting at the
center of the first face of the mesh in an arbitrary direction:

  vector<Vector> curve, normals;
  FaceIter f = mesh.faces.begin();
  double initialAngle = 0.;
  Vector initialPoint = f->toLocal( f->barycenter() );
  mesh.integralCurve( initialFace,
                      initialPoint,
                      initialAngle,
                      curve,
                      normals );

The resulting vectors "curve" and "normals" contain lists of positions and
surface normals along the curve.  This routine uses piecewise linear
interpolation of the direction field and may not produce nice results in highly
singular regions.


-------------
8. LICENSE
-------------

TCODS is covered by the following free software license based on the 2-clause
BSD license.  This license is compatible with most other free software
licenses, including the GNU General Public License.

*
* Copyright 2010 Keenan Crane. All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without modification,
* are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, this
*    list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
* WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
* SHALL THE FREEBSD PROJECT OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
* OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
* ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* The views and conclusions contained in the software and documentation are those
* of the author and should not be interpreted as representing official policies,
* either expressed or implied, of any other person or institution.
*

