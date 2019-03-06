# tcods v0.1 Makefile
# -------------------

#
# You must specify paths to METIS, SuiteSparse, and BLAS if they
# are not kept in the default paths (e.g., /usr/local/lib, etc.).
# You must also specify a BLAS library (on MacOS X you can simply
# use -framework Accelerate)
#

# library and include paths for METIS
METIS_LIB =
METIS_INC =

# library and include paths for SuiteSparse
SUITESPARSE_LIB =
SUITESPARSE_INC =

# library paths and flags for BLAS and LAPACK
BLAS_INC =
BLAS_LIB = -framework Accelerate
LAPACK_INC =
LAPACK_LIB =

################################################################################

CPP = g++
LD = g++
CPPFLAGS = -O3 -Wall -Werror
LDFLAGS = -O3 -Wall -Werror
LIBS = $(METIS_LIB) $(SUITESPARSE_LIB) $(BLAS_LIB) $(LAPACK_LIB) -lm -lamd -lcamd -lcolamd -lccolamd -lcholmod -lspqr -lmetis
INCLUDE = $(METIS_INC) $(SUITESPARSE_INC) $(BLAS_INC) $(LAPACK_INC)
OBJS = main.o HalfEdge.o Vector.o MeshIO.o Connection.o Problem.o CMWrapper.o
TARGET = tcods

all: $(OBJS)
	$(LD) $(LDFLAGS) $(LIBS) *.o -o $(TARGET)

main.o: main.cpp HalfEdge.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) main.cpp

HalfEdge.o: HalfEdge.cpp HalfEdge.h Vector.h MeshIO.h Connection.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) HalfEdge.cpp

Vector.o: Vector.cpp Vector.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) Vector.cpp

MeshIO.o: MeshIO.cpp MeshIO.h HalfEdge.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) MeshIO.cpp

Connection.o: Connection.cpp Connection.h HalfEdge.h CMWrapper.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) Connection.cpp

Problem.o: Problem.cpp Problem.h HalfEdge.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) Problem.cpp

CMWrapper.o: CMWrapper.cpp CMWrapper.h HalfEdge.h
	$(CPP) -c $(CPPFLAGS) $(INCLUDE) CMWrapper.cpp

clean:
	rm -f $(OBJS)
	rm -f $(TARGET)

