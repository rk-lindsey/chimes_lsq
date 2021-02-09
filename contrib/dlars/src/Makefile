CXX=mpicxx -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DWARN_ONLY
##CXX=mpicxx -O3 -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DUSE_OPENMP
##CXX=mpicxx -g -fopenmp -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DDEBUG
##CXX=mpicxx -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI
##CXX=g++ --std=c++11 -g -Wall -DDEBUG
#CXX=g++ --std=c++11 -O3 -Wall

# Before compiling, switch to gcc/8.1.0.
#	module load gcc/8.1.0
dlars: dlars.C Vector.h IntVector.h Matrix.h DLARS.h
	$(CXX) dlars.C -lm -o dlars

realclean:
	rm dlars



