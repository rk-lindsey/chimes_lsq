# Before running, be sure to export hosttype=<your machine's modfile>
# Take a look in ./modfiles to see your options. Otherwise, your machine's default modfiles will be used.
#
# NOTE: Make sure the right compiler is specified in the Makefile... don't use MPI for lsq.
#

###############################################################
#
# Determine the location of necessary files
#
###############################################################

# Common function for test script initialization.
source ../src/bash/init_vars.sh
init_test_vars
echo "Number of processors (unless you really know what you're doing, set to 36) = $NP"
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

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

if [ $# -eq 0 ] ; then # Use default JOBS.  
	#  h2o-3bcheby2' -- gives a diff answer than old code b/c of layer bug in old code
    JOBS=$LSQ_ALL_JOBS
    MAKE_JOBS=$LSQ_MAKE_JOBS
else  # Take JOBS from command line.
    JOBS=$1
    MAKE_JOBS=$2
fi

echo ""
echo "SETTING UP FOR CHIMES_LS..."

if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
    RUN_JOB=""
fi

for i in $JOBS
do
	# Strange things happen when NP >> NF. These below 4b
	# test uses only 10 frames because otherwise the test
	# would take forever 

	if ! test_dir $i ; then continue ; fi

	cd $i

	if [ ! -d current_output ] ; then mkdir current_output ; fi
	if [ ! -d correct_output ] ; then mkdir correct_output ; fi

	cd current_output
	cp ../*.txt ../*.in ../*.xyzf ../*.dat ./ >& /dev/null	

	if $RUN_JOB ../../../build/chimes_lsq fm_setup.in > fm_setup.out ; then
		 echo "Chimes_lsq succeeded"
		 SUCCESS=1
	else
		 echo "Chimes_lsq failed"
		 SUCCESS=0
	fi	
        
	if [[ $SUCCESS -eq 1 ]] ; then
	
		# Break up  A.txt file into 95M chunks, i.e., < Github's 100M file size limit
		# Then remove the big single A.txt file
		
		split -b95M A.txt A.txt.
	
 		 cp A.txt.* b.txt params.header fm_setup.out ff_groups.map ../correct_output
	fi
	
	cd ../..
done

###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "SETTING UP FOR SVD SCRIPT..."

for i in  $JOBS
do

	if ! test_dir $i ; then
		 continue 
	fi

	cd $i/current_output
	
	
	if $RUN_LSQ_PYTHON_CODE > params.txt ; then
		 echo "LSQ code succeeded"
		 SUCCESS=1
	else
		 echo "LSQ code failed"
		 SUCCESS=0
	fi
	
	for j in params.txt force.txt
	do
	
		if [ "$j" == params.txt ]; then
			j=params.txt-tailed
			tail -n+2 params.txt > $j
		fi
		
	done
	

	cd ../..
	
	if [[ $SUCCESS -eq 1 ]] ; then
		 cp $i/current_output/* $i/correct_output
	fi
done

echo "Running Makefile jobs $MAKE_JOBS"


for job in $MAKE_JOBS ; do

    echo "Running $job"
    
    if ! test_dir $job ; then
	echo "Directory $job does not exist"
	continue 
    fi

    cd $job
	 
    if make RUN_JOB="$RUN_JOB" PYTHON=$PYTHON generate ; then
	echo "$job succeeded"
    else
	echo "$job failed"
    fi
    
    cd -
done

exit 0
