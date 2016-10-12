##################################
# CHOOSE COMPILATION TYPE, SET FLAGS
##################################

# Compiling flags

MYFLAGS= -D VERBOSITY=1 -D FPENALTY_POWER=3.0 # -D XFORMTYPE=LEFT # -D FORCECHECK=1 -D PESFORMAT=1 -D LOG_FORCE -D LOG_PES

CXX=mpiicc  -g -D USE_MPI 	# For MPI (intel  ). Run on slurm machines with: srun -n <#procs> house_md < run_md.in
#CXX=mpicc   -g -D USE_MPI 	# For MPI (gnu    ). Run on slurm machines with: srun -n <#procs> house_md < run_md.in
#CXX=mpicxx  -g  -D USE_MPI 	# For MPI (gnu/osx). Run on slurm machines with: srun -n <#procs> house_md < run_md.in
#CXX=icc -g			# Serial, for debugging
#CXX=g++ -O3			# Serial, let compiler optimize to a very high level (see also -O2)
#
# Below: Other debugging/optimiazation options
#
#CXX=g++ -g -Wall -Wextra -Wuninitialized -O1
#CXX=g++ -g # -Wall -Wextra
#CXX=icpc -fast -Wall

%.o: %.C
	$(CXX) $(MYFLAGS) -c $< -o $@


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
# Make Travis analyzer
##################################	

travis:
	make -C ../contrib/Travis

clean_travis:
	make clean -C ../contrib/Travis

##################################
# Make packmol initial configuration generator
##################################	

packmol:
	make -C ../contrib/packmol

clean_packmol:
	make clean -C ../contrib/packmol
