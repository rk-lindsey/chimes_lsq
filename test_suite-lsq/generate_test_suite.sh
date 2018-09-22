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
make -f Makefile-TS-LSQ chimes_lsq
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
else
# Take JOBS from command line.
	 JOBS=$*
fi

echo ""
echo "SETTING UP FOR SPLINES_LS..."

for i in $JOBS
do

	echo " "
	echo "Running $i test..."

	cd $i
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		 ../chimes_lsq < fm_setup.in > fm_setup.out
	else
		 $RUN_JOB ../chimes_lsq < fm_setup.in > fm_setup.out
	fi	
 	cp A.txt b.txt params.header fm_setup.out ff_groups.map correct_output	
	mv A.txt b.txt params.header fm_setup.out ff_groups.map current_output
	
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
	
	$RUN_LSQ_PYTHON_CODE > params.txt
	
	for j in params.txt force.txt
	do
	
		if [ "$j" == params.txt ]; then
			j=params.txt-tailed
			tail -n+2 params.txt > $j
		fi
		
	done
	

	cd ../..
	
	cp $i/current_output/* $i/correct_output
done

exit 0
