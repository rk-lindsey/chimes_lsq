##CXX=mpicxx -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI
##CXX=mpicxx -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DDEBUG
CXX=g++ --std=c++11 -g -Wall -DDEBUG
##CXX=g++ --std=c++11 -O3 -Wall

dlars: dlars.C Vector.h IntVector.h Matrix.h DLARS.h
	$(CXX) dlars.C -lm -o dlars

realclean:
	rm dlars



