CXX=g++ -g -std=c++11 -Wall -Wno-sign-compare -DENABLE_FP_EXCEPT
CXX_SRC=OWLQN.cpp TerminationCriterion.cpp leastSquares.cpp logreg.cpp main.cpp

owlqn : OWLQN.cpp TerminationCriterion.cpp leastSquares.cpp logreg.cpp main.cpp
	$(CXX) --std=c++11 -Wall $(CXX_SRC) -lm -o owlqn

clean :
	rm owlqn
