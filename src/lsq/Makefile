# Choose g++ for debugging and icpc for fast execution.
#
CXX=g++ -g -Wall #-pg
#CXX=icpc -fast	
CXX_SRC=splines_ls.C functions.C
CXX_INCLUDES=

all : splines_ls

splines_ls : $(CXX_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_INCLUDES) -o splines_ls
	cp splines_ls ../../bin/splines_ls

%.o: %.C
	$(CXX) -Wno-deprecated $(CXX_INCLUDES) -c $< -o $@

clean:
	/bin/rm -f *.o a.out

realclean : clean
	rm -f splines_ls





