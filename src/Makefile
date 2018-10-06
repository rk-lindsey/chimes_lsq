##################################
# CHOOSE COMPILATION TYPE, SET FLAGS
##################################

# For linking to LAMMPS:

# I get strange behavior when I use intel/mvapich/mpicxx... I use intel/impi/mpiicc... check which other edits were made (Besides restart and lsq)


LINK2LMPS=no
LAMMPSDIR=/usr/gapps/polymers/lammps_borax/

# Compiling flags

MYFLAGS= -D VERBOSITY=1 -D FPENALTY_POWER=3.0 

# Other options:
#
# -D XFORMTYPE=LEFT 
# -D FORCECHECK=1 
# -D PESFORMAT=1 
# -D LOG_FORCE 
# -D LOG_PES

# Compilers

# Use mpicxx in conjunction with use mvapich2-intel-2.2 for intel compiler + infiniband (fastest) MPI.
# Type these commands before using this Makefile:
# Chaos 5 commands:
# use ic-17.0.174
# use mvapich2-intel-2.2
#
# TOSS-3 commands:
# module load intel/17.0.2
# module load mvapich2/2.2
#
#CXX=mpicxx  -std=c++11  -g -pg -D USE_MPI 

#CXX=mpiicc -std=c++11  -O3 -D USE_MPI 	# For MPI (intel  ). Run on slurm machines with: srun -n <#procs> chimes_md < run_md.in
#CXX=mpicxx  -std=c++11  -g -Wall -D USE_MPI 	# For MPI (gnu    ). Run on slurm machines with: srun -n <#procs> chimes_md < run_md.in
CXX=mpicxx -std=c++11 -O3 -D USE_MPI 	# For MPI (gnu/osx). Run on slurm machines with: srun -n <#procs> chimes_md < run_md.in
#CXX=mpicxx -std=c++11 -g -D USE_MPI 	# For MPI (gnu/osx). Run on slurm machines with: srun -n <#procs> chimes_md < run_md.in
#CXX=icc -g				# Serial, for debugging
#CXX=icc -O3
#CXX=g++ -std=c++11 -g			# Serial, let compiler optimize to a very high level (see also -O2)
#CXX=g++ -std=c++11 -O3			# Serial, let compiler optimize to a very high level (see also -O2)
#
# Below: Other debugging/optimiazation options
#
#CXX=clang++ -g -std=c++11 -g -Wall
#CXX=g++ -g -std=c++11 -Wall -Wno-sign-compare -DENABLE_FP_EXCEPT

#CXX=icpc -fast -Wall

#CXX=clang -std=c++11

##################################
# MAKE THE LSQ CODE
##################################

CXX_LSQ_SRC=chimes_lsq.C functions.C Ewald.C ClassDefs.C Fcut.C Cluster.C util.C Cheby.C io_styles.C feexcept.C A_Matrix.C
CXX_LSQ_HEADER=functions.h Fcut.h Cluster.h util.h Cheby.h io_styles.h A_Matrix.h

chimes_lsq: $(CXX_LSQ_SRC) $(CXX_LSQ_HEADER)
	$(CXX) $(CXX_LSQ_SRC) $(CXX_LSQ_INCLUDES) -o chimes_lsq

clean_lsq:
	/bin/rm -f chimes_lsq.o functions.o Ewald.o

realclean_lsq: clean_lsq
	rm -f chimes_lsq 

build_tests_lsq:
	cd ../test_suite-lsq
	./generate_test_suite.sh
	cd -

##################################
# MAKE THE MD CODE
##################################

CXX_MD_SRC=chimes_md.C functions.C Ewald.C ClassDefs.C Fcut.C Cluster.C util.C Cheby.C io_styles.C feexcept.C
CXX_MD_HEADER=functions.h lmppath.h Fcut.h Cluster.h util.h Cheby.h io_styles.h 

FLAGLMP=-D LINK_LAMMPS
LAMMPS_LIB=$(LAMMPSDIR)/src/liblammps.a
LMP_COUPLE_LIB=$(LAMMPSDIR)/examples/COUPLE/library/libcouple.a
LIBS=/usr/lib64/libfftw3.a

chimes_md: $(CXX_MD_SRC) $(CXX_MD_HEADER)
ifeq ($(LINK2LMPS), yes)
	$(CXX) $(FLAGLMP) $(CXX_MD_SRC) $(LAMMPS_LIB) $(LMP_COUPLE_LIB) $(LIBS) -o chimes_md
else
	$(CXX) $(CXX_MD_SRC) -o chimes_md
endif

clean_md:
	/bin/rm -f *.o

realclean_md: clean_md
	/bin/rm -f chimes_md

build_tests_md:
	cd ../test_suite-md
	./generate_test_suite.sh
	cd -

##################################
# MAKE BOTH THE LSQ AND MD CODES
##################################


all:
	make chimes_md
	make chimes_lsq

tests:
	make all
	build_tests_lsq
	build_tests_md	
clean:
	 rm -f *.o

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

tags:
	./maketags.pl
