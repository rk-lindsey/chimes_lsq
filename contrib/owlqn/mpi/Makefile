## Note: -DDEBUG will significantly slow calculations.
##CXX=mpicxx -g -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI -DDEBUG
CXX=mpicxx -O3 -std=c++11 -Wall -Wno-sign-compare -DUSE_MPI
CXX_SRC=OWLQN.cpp TerminationCriterion.cpp leastSquares.cpp logreg.cpp main.cpp
CXX_HEADERS=OWLQN.h TerminationCriterion.h  leastSquares.h  logreg.h  mpiexterns.h
.o.c:
	$(CXX) -c $*.c

dowlqn : $(CXX_SRC) $(CXX_HEADERS)
	$(CXX) --std=c++11 -Wall $(CXX_SRC) -lm -o dowlqn

clean :
	rm dowlqn
