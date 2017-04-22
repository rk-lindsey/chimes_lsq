# Choose g++ for debugging and icpc for fast execution.
#
#CXX=g++ -g -Wall -Wextra -Wuninitialized -O1#-pg
#CXX=g++ -g -Wall -Wextra 
CXX=g++ -O3
#CXX=icc -O3
#CXX=icpc -fast
CXX_SRC=splines_ls.C ../mdcode/functions.C ../mdcode/Ewald.C ../mdcode/Cheby3B.C
CXX_INCLUDES=-I ../mdcode

all : splines_ls

splines_ls : $(CXX_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_INCLUDES) -o splines_ls
	cp splines_ls ../../bin/splines_ls

%.o: %.C
	$(CXX) $(CXX_INCLUDES) -c $< -o $@

clean:
	/bin/rm -f *.o

realclean : clean
	rm -f splines_ls





