#! /bin/bash
# Usage:  With no arguments, all tests are run.  Otherwise, only tests specified on the command line are run.
#  
# Before running, be sure to export hosttype=<your machine's modfile>
# Take a look in ./modfiles to see your options. Otherwise, your machine's default modfiles will be used.
#         run_test_suite.sh 'lsq-jobs' 'make-jobs'
#         'lsq-jobs' is a list of lsq tests that are run by the script.  This argument may be an empty string.
#         'make-jobs' is a list of lsq tests that are run by makefiles.  This argument may be an empty string.
#

#######
#
# Determine which system you're running on and load the necessary modules.
# If a modfile is found, assume it is an HPC system with 36 procs. Otherwise, assume serial
#
######

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
#  Run the tests for the chimes_lsq program
#
###############################################################

echo ""
echo "VALIDATING FOR CHIMES_LSQ..."

SET_PASSED=true
SVD_PASSED=true
ALL_PASSED=true

for i in $JOBS
do

	 if ! test_dir $i ; then
		  continue 
	 fi
		  
	 PASS=true

	 
	 if [[ ! -d $i/current_output ]] ; then
		  mkdir $i/current_output
	 fi

	 cd $i/current_output

	 cp ../*.txt ../*.in ../*.xyzf ../*.dat ./ >& /dev/null
	 
	 if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
	     if ../../../build/chimes_lsq fm_setup.in > fm_setup.out ; then
		 echo 'Chimes_lsq succeeded'
		 SUCCESS=1
	     else
		 echo 'Chimes_lsq failed'
		 SUCCESS=0
	     fi
	else
	    if $RUN_JOB ../../../build/chimes_lsq fm_setup.in  | awk '!/TACC:/{print}' > fm_setup.out ; then
		echo 'Chimes_lsq succeeded'
		SUCCESS=1
	    else
		echo 'Chimes_lsq failed'
		SUCCESS=0
		PASS=false
		ALL_PASSED=false
	    fi
	fi

        # There seems to be an NFS filesystem lag in finding output files.
	cd ..
	if [[ $SUCCESS -eq 1 ]] ; then
	
		sleep 5
	
		 # First, we need to reconstruct the A.txt file in correct output
		 
		 cat correct_output/A.txt.* > correct_output/A.txt

		 for j in A.txt b.txt params.header fm_setup.out
		 do
		
			  perl ../../contrib/compare/compare.pl correct_output/$j current_output/$j > diff-$j.out
	
			  NO_DIFF_LINES=`cat diff-$j.out | wc -l`
	
			  if [ $NO_DIFF_LINES -gt 0 ] ; then
					echo " "
					echo "Differences found in $j files... "
					echo " "

					if [ "$j" = A.txt ] ; then
						 echo "	...differences are in file A.txt. Examine file by hand "
						 echo "     to determine whether differences are within reasonable"
						 echo "     amt (i.e. < 10**-10)"
						 echo "     Otherwise, see the *.diff file in ${i}/correct_output/"
					else
						 diff correct_output/$j current_output/$j
					fi
					echo " "
					PASS=false
			  fi
		 done
	
		 if [ "$PASS" = true ] ; then
			  echo "		...Test passed."
			  rm -f diff-*
		 else
			  SET_PASSED=false
			  ALL_PASSED=false
			  echo "		...Test failed."
		 fi
		 
		 rm -f correct_output/A.txt
	fi 
	cd ..
done

###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "VALIDATING FOR SVD SCRIPT ..."

for i in $JOBS
do

	if ! test_dir $i ; then
		 continue 
	fi

	PASS=true
	TECHNICAL_PASS=true

	cd $i/current_output
	rm -rf diff-*
	if $RUN_LSQ_PYTHON_CODE > params.txt ; then
		 echo "LSQ code succeeded"
		 SUCCESS=1 
	else
		 echo "LSQ code failed"
		 SUCCESS=0
		 PASS=false
		 SVD_PASSED=false
	fi

	if [[ $SUCCESS -eq 1 ]] ; then
		 for j in params.txt force.txt
		 do
	
			  if [ "$j" == params.txt ]; then
					j=params.txt-tailed
					awk '!/Date/{print}' params.txt > $j
					awk '!/Date/{print}' ../correct_output/params.txt > ../correct_output/$j
			  fi
			  
		
			  # Ignore the date when running diff...

			  diff ../correct_output/$j $j > ../diff-$j.out
	
			  NO_DIFF_LINES=`cat ../diff-$j.out | wc -l`
	
			  if [[ $NO_DIFF_LINES -gt 0 ]] ; then
					echo " "
					echo "Differences found in $j files:"
					echo " "
					if [ "$j" == "force.txt" ] ; then
						 paste ../correct_output/force.txt force.txt > check_tol.dat
						 awk 'BEGIN{tol=10^-7;any=0}{val=$1-$2;   if(sqrt(val*val)>=tol) {any++;print ("	Parameter index", $1, " differences exceeded tolerance(+/-", tol, "): ",val)}}END{if(any==0){print ("	No parameters differ by more than tol (",tol,").")}}' check_tol.dat | tee tol_status.dat
						 echo "	Check file ../diff-$j.out for any other (non-parameter) differences."
						 rm -f check_tol.dat			
					else
						 # See if the differences exist between parameters... 
						 # if so, tell the user what the tolerance is, and
						 # direct them to the diff file to look at to make
						 # sure that is the only difference
				
						 paste ../correct_output/test_suite_params.txt test_suite_params.txt > check_tol.dat
						 awk 'BEGIN{tol=10^-9;any=0}{val=$2-$4;   if(sqrt(val*val)>=tol) {any++;print ("	Parameter index", $1, " differences exceeded tolerance(+/-", tol, "): ",val)}}END{if(any==0){print ("	No parameters differ by more than tol (",tol,").")}}' check_tol.dat  | tee tol_status.dat
						 rm -f check_tol.dat
					fi
					echo " "
			
					TECHNICAL_PASS_STATUS=`grep "No" tol_status.dat`
			
					if [[ "$TECHNICAL_PASS_STATUS" != *"No"* ]]; then
						 TECHNICAL_PASS=false
					fi
					
					PASS=false
			  fi
		
			  if [ "$j" == params.txt ]; then
					rm -rf $j
					rm -rf ../correct_output/$j
			  fi
		 done
	
		 if [ "$PASS" = true ] ; then
			  echo "		...Test passed."
			  rm -f ../diff-*
		 else
			  SVD_PASSED=false
			  ALL_PASSED=false
		
			  if [ "$TECHNICAL_PASS" = true ]; then
					echo "		...Technical test pass - all parameters within tolerances, but other aspects may differ."
			  else
					echo "		...Test failed."
			  fi
		 fi

		 cd ../..
	fi
done

echo "Running Makefile jobs $MAKE_JOBS"
pwd ;

for job in $MAKE_JOBS ; do

    echo "Running $job"
    
    if ! test_dir $job ; then
	echo "Directory $job does not exist"
	continue 
    fi
	 
    cd $job
	 
    if make RUN_JOB="$RUN_JOB" PYTHON=$PYTHON all ; then
	echo "$job succeeded"
    else
	echo "$job failed"
	ALL_PASSED=FALSE
	PASSED=FALSE
    fi

   cd -
done

if   [ "$ALL_PASSED" = true ] ; then
	echo "ALL TESTS PASSED"
elif [ "$ALL_PASSED" = false ] ; then
	echo 'AT LEAST ONE EACH OF SETUP AND SVD TEST(S) FAILED...'
	echo "Check individual results above to confirm if test had technical pass."
elif [ "$SVD_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SVD SCRIPT"
elif [ "$SET_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SETUP SCRIPT"
else
	echo 'ERROR: BAD LOGIC IN TEST SUITE DRIVER (THIS SCRIPT)'
	echo "ALL/SVD PASSED: $ALL_PASSED $SVD_PASSEDs" 
fi
	
exit 0
