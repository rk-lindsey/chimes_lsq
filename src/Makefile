##################################
# CHOOSE COMPILATION TYPE, SET FLAGS
##################################

# Quick note on flag PESFORMAT:
# 1. gnuplot style
# 2. matlab style
# 3. python style
# Currently, only gnuplot style is supported

MYFLAGS= -D VERBOSITY=1 -D FORCECHECK=1 -D PESFORMAT=1

CXX=g++ -g  $(MYFLAGS) 
#CXX=g++ -O3  $(MYFLAGS) 
#CXX=g++ -g -Wall -Wextra -Wuninitialized -O1
#CXX=g++ -g # -Wall -Wextra
#CXX=icpc -fast -Wall



##################################
# MAKE THE LSQ CODE
##################################

CXX_LSQ_SRC=splines_ls.C functions.C Ewald.C

CXX_LSQ_INCLUDES=

house_lsq: $(CXX_LSQ_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_LSQ_INCLUDES) -o house_lsq

	
clean_lsq:
	/bin/rm -f splines_ls.o functions.o Ewald.o

realclean_lsq: clean_lsq
	rm -f house_lsq 
	
build_tests_lsq:
	cd ../test_suite-lsq
	./generate_test_suite.sh
	cd -

##################################
# MAKE THE C++ CODE
##################################

CXX_MD_SRC=splines_md.C functions.C Ewald.C

CXX_MD_INCLUDES=

house_md: $(CXX_MD_SRC:%.C=%.o)
	$(CXX) $^ $(CXX_MD_INCLUDES) -o house_md


clean_md:
	/bin/rm -f splines_md.o functions.o Ewald.o

realclean_md: clean_md
	rm -f house_md 

build_tests_md:
	cd ../test_suite-md
	./generate_test_suite.sh
	cd -
##################################
# MAKE BOTH THE LSQ AND MD CODES
##################################


all:
	make house_md
	make house_lsq
	
tests:
	make all
	build_tests_lsq
	build_tests_md	
clean:
	 rm -f *o
	 
##################################
# MAKE Travis analyzer
##################################	

travis:
	make -C ../contrib/Travis
	
clean_travis:
	make clean -C ../contrib/Travis
	
