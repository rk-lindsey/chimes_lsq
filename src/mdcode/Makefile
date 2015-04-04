CXX=g++ -g -Wall
##CXX=icpc -fast -Wall
CXX_SRC=splines_md.C functions.C Ewald.C
CXX_INCLUDES=

splines_md: $(CXX_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_INCLUDES) -o splines_md
	cp splines_md ../../bin

all : splines_md makeparam

%.o: %.C
	$(CXX) -Wno-deprecated $(CXX_INCLUDES) -c $< -o $@

clean:
	/bin/rm -f *.o

realclean: clean
	rm -f splines_md makeparam

# The makeparam program generates spline parameters for a Lennard-Jones potential.
# This is used for testing purposes only.
makeparam: makeparam.C
	$(CXX) makeparam.C -o makeparam
