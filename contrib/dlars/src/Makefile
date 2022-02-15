CXX=mpicxx
##CXX=mpicxx -O3 -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DUSE_OPENMP
##CXX=mpicxx -g -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DDEBUG
##CXX=mpicxx -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI
##CXX=g++ --std=c++11 -g -Wall -DDEBUG
#CXX=g++ --std=c++11 -O3 -Wall

## Use intel threads.
## module load intel/oneapi.2021.3
## module load mkl/2020.0
## Add -DTIMING to CFLAGS for detailed timing info output.
CFLAGS= -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY  -DMKL_ILP64  -I"${MKLROOT}/include" -DUSE_MKL -DUSE_BLAS -qopenmp -DUSE_OPENMP
#CFLAGS= -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY  -DMKL_ILP64  -I"${MKLROOT}/include" -DTIMING -qopenmp
#CFLAGSDBG= -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY  -DMKL_ILP64  -I"${MKLROOT}/include" -DTIMING -qopenmp -DUSE_OPENMP
CFLAGSDBG= -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY  
LINKFLAGS= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

opt: dlars.C Matrix.C Vector.h IntVector.h Matrix.h DLARS.h
	$(CXX) $(CFLAGS) dlars.C Matrix.C $(LINKFLAGS) -o dlars

debug: dlars.C Vector.h IntVector.h Matrix.h DLARS.h
	$(CXX) $(CFLAGSDBG) dlars.C Matrix.C $(LINKFLAGS) -o dlars

clean:
	rm dlars



