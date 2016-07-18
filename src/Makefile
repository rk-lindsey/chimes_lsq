# FYI: clean removes .o files, realclean removes executable



# Choose g++ for debugging and icpc for fast execution.
#
#CXX=g++ -Wall -pg # NOTE: pg flag is for profiling with gprof (OSX compatible)
#CXX=g++ -g -Wall -Wextra -Wuninitialized -O1#-pg
#CXX=g++ -g   # -Wall -Wextra 
CXX=g++ -O2   # -Wall -Wextra 
#CXX=g++ -O3   # -Wall -Wextra 
#CXX=icpc -fast


CXX_SRC=splines_ls.C functions.C Ewald.C #  ../mdcode/Cheby3B.C
CXX_INCLUDES=-I .

all : splines_ls

splines_ls : $(CXX_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_INCLUDES) -o splines_ls
#	cp splines_ls ../../bin/splines_ls

%.o: %.C
	$(CXX) $(CXX_INCLUDES) -c $< -o $@

clean:
	/bin/rm -f *.o

realclean : clean
	rm -f splines_ls





