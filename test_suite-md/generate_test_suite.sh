#!/bin/bash


# Before running, be sure to export hosttype=<your machine's modfile>
# Take a look in ./modfiles to see your options. Otherwise, your machine's default modfiles will be used.
#
# NOTE: Make sure the right compiler is specified in the Makefile... don't use MPI for lsq.
#

########################################
# Define tests within the test suite
########################################

# Tests specifically for the MD code

# Tests for compatibility between LSQ C++/python codes with the MD code
TAG="verify-lsq-forces-"

# Iterate through the tests

echo " "
echo "SETTING UP FOR MD CODE..."
echo " "


###############################################################
#
# Determine the location of necessary files
#
###############################################################

# Common function for test script initialization.
source ../src/bash/init_vars.sh
init_test_vars
echo "Number of processors (unless you really know what you're doing, set to 23) = $NP"
echo "Warning: this script will overwrite all reference test output."
echo "Continue ? (yes/no)"

read ok_run
if test "x$ok_run" != "xyes" ; then
	echo 'Quitting'
	exit 0
fi
	

NP=36
TESTSU_BASE=`pwd -P` #`dirname $0`
SOURCE_BASE="${TESTSU_BASE}/../build/"

RUN_JOB="srun -n $NP"

if [ ! -v hosttype ] ; then
    echo "No hosttype specified"
    echo "Be sure to load modules/configure compilers by hand before running this script."
    echo "Otherwise, run with export hosttype=<host type>; ./this_script.sh"
    NP=1
    RUN_JOB=""
elif [[ "$hosttype" == "LLNL-LC" ]] ; then
    source ${TESTSU_BASE}/../modfiles/LLNL-LC.mod
    ICC=`which icc`    
    MPI=`which mpicxx`    
elif [[ "$hosttype" == "UM-ARC" ]] ; then
    source ${TESTSU_BASE}/../modfiles/UM-ARC.mod
    ICC=`which icc`    
    MPI=`which mpicxx`    
elif [[ "$hosttype" == "JHU-ARCH" ]] ; then
    source ${TESTSU_BASE}/../modfiles/JHU-ARCH.mod
    ICC=`which icc`
    MPI=`which mpicxx`   
elif [[ "$hosttype" == "UT-TACC" ]] ; then
    source ${TESTSU_BASE}/../modfiles/UT-TACC.mod
    RUN_JOB="ibrun" 
else
    echo ""
    echo "ERROR: Unknown hosttype ($hosttype) specified"
    echo ""
    echo "Valid options are:"
    for i in `ls modfiles`; do echo "	${i%.mod}"; done
    echo ""
    echo "Please run again with: export hosttype=<host type>; ./install.sh"
    echo "Or manually load modules and run with: ./this_script.sh"
    exit 0
fi

NUM_THREADS=$NP		# Number of threads for SVD decomposition
export OMP_NUM_THREADS=$NUM_THREADS    



###############################################################
#
# Determine the location of necessary files
#
###############################################################


# Common function for test script initialization.
source ../src/bash/init_vars.sh
DLARS_PATH=../contrib/dlars/src

init_test_vars
echo "NP = $NP"

if [ $# -eq 0 ] 
then
	 JOBS=$LSQ_ALL_JOBS
	 MAKE_JOBS=$LSQ_MAKE_JOBS
else
	 JOBS=$1
	 MAKE_JOBS=$2
fi


# Ensure MKL is available for DLARS

if [ ! -v hosttype ] ; then
	echo "Will not run make jobs: "
	echo "Automated DLARS compilation currently requires access to "
	echo "a Livermore Computing system"
	MAKE_JOBS=""
fi



###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ..

if ./install.sh  ; then
    echo "Compiling chimes_lsq succeeded"
else
    echo "Compiling chimes_lsq failed"
    exit 1
fi

cd -

# Determine which tests to run

if [ $# -gt 0 ] ; then
	 MD_JOBS=$1
    LSQ_FORCE_JOBS=$2
fi
echo "MD JOBS = $MD_JOBS"
echo "LSQ_FORCE_JOBS = $LSQ_FORCE_JOBS"

for i in $MD_JOBS
do

	echo " "
	echo "Running $i test..."
	
	
	if [ ! -d $i/current_output ] ; then mkdir $i/current_output ; fi
	if [ ! -d $i/correct_output ] ; then mkdir $i/correct_output ; fi

	cp $i/* $i/current_output 2> /dev/null
	cd $i/current_output
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		 if ../../../build/chimes_md run_md.in > run_md.out ; then
			  SUCCESS=1
		 else
			  echo "Chimes_md failed"
			  SUCCESS=0
		 fi
			
	else
		 if $RUN_JOB ../../../build/chimes_md run_md.in > run_md.out ; then
			  SUCCESS=1
		 else
			  echo "Chimes_md failed"
			  SUCCESS=0
		 fi
			  
	fi		

	if [[ $SUCCESS -eq 1 ]] ; then
		 cp *.* ../correct_output
	fi
	
	cd -
done


if [ -n "$LSQ_FORCE_JOBS" ] ; then

echo " "
echo "SETTING UP FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../test_suite-lsq 

for i in ${LSQ_FORCE_JOBS}
do
	./run_test_suite.sh $i # $LSQ_FORCE_JOBS
done

cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in ${LSQ_FORCE_JOBS}
do

	echo $i
	
	if [[ $i == *"lsq2"* ]] ; then
		continue
	fi 

	echo " "
	echo "Running $i test..."
	
	cd ${TAG}${i}

	if [ ! -d current_output ] ; then mkdir current_output ; fi
	if [ ! -d correct_output ] ; then mkdir correct_output ; fi
	
	# Grab the parameter and force files from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt    .
	cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
	cp ../../test_suite-lsq/$i/current_output/force.txt     .
	
	if ../../build/chimes_md run_md.in > run_md.out ; then
		 echo 'Chimes_md succeeded'
		 SUCCESS=1
	else
		 echo 'Chimes_md failed'
		 SUCCESS=0
	fi
	
	cp *.* current_output
	if [[ $SUCCESS -eq 1 ]] ; then
		 cp *.* correct_output
	fi
	
	cd -
done	


echo " "

fi # -n LSQ_JOBS
	
exit 0
