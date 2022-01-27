CXX=mpicxx
##CXX=mpicxx -O3 -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DUSE_OPENMP
##CXX=mpicxx -g -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DDEBUG
##CXX=mpicxx -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI
##CXX=g++ --std=c++11 -g -Wall -DDEBUG
#CXX=g++ --std=c++11 -O3 -Wall

# Before compiling, switch to gcc/8.1
#	module load gcc/8.1.0
## module load intel/oneapi.2021.3
## module load mkl/2020.0
CFLAGS= -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY  -DMKL_ILP64  -I"${MKLROOT}/include" -DUSE_MKL -DUSE_BLAS
LINKFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

dlars: dlars.C Vector.h IntVector.h Matrix.h DLARS.h
	$(CXX) $(CFLAGS) dlars.C $(LINKFLAGS) -o dlars

realclean:
	rm dlars



