# NOTE: The path below needs to point to the "lsq.py" version of the lsq code.
#       If you want to use other versions, you'll probably need to modify the inputs that
#       are sent to the script (way down below)
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
echo "NP = $NP"

# Run the job with the new version of the python code (Compatible with non-generalized md code)
#
###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
rm -f *.o chimes_lsq
if make -f Makefile-TS-LSQ chimes_lsq ; then
	 echo "Compiling chimes_lsq succeeded"
else
	 echo "Compiling chimes_lsq failed"
	 exit 1
fi

mv chimes_lsq ../test_suite-lsq/
cd ../test_suite-lsq

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

if [ $# -eq 0 ] 
# Use default JOBS.  
then
#  h2o-3bcheby2' -- gives a diff answer than old code b/c of layer bug in old code
	 JOBS=$LSQ_ALL_JOBS
	 MAKE_JOBS=$LSQ_MAKE_JOBS
else
# Take JOBS from command line.
	 JOBS=$*
	 MAKE_JOBS=""
fi

echo ""
echo "SETTING UP FOR SPLINES_LS..."

if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
	 RUN_JOB=""
fi

for i in $JOBS
do

	echo " "
	echo "Running $i test..."

	cd $i
	

	if $RUN_JOB ../chimes_lsq < fm_setup.in > fm_setup.out ; then
		 echo "Chimes_lsq succeeded"
		 SUCCESS=1
	else
		 echo "Chimes_lsq failed"
		 SUCCESS=0
	fi	

	if [[ $SUCCESS -eq 1 ]] ; then
 		 cp A.txt b.txt params.header fm_setup.out ff_groups.map correct_output	
		 mv A.txt b.txt params.header fm_setup.out ff_groups.map current_output
	fi
	
	cd ..
done


###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "SETTING UP FOR SVD SCRIPT..."

for i in  $JOBS
do
	echo " "
	echo "Running $i test..."
	echo "Using command $RUN_LSQ_PYTHON_CODE"

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
	 cd $job
	 if make generate ; then
		  echo "$job succeeded"
	 else
		  echo "$job failed"
	 fi
done

exit 0
