#! /bin/sh
# NOTE: The path below needs to point to the "lsq-new-md-fmt.py" version of the lsq code.
#       If you want to use other versions, you'll probably need to modify the inputs that
#       are sent to the script (down below)
#
# Usage:  With no arguments, all tests are run.  Otherwise, only tests specified on the command line are run (LEF).

###############################################################
#
# Determine the location of necessary files
#
###############################################################

if [ $# -eq 0 ] 
then
JOBS='h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby' #  h2o-3bcheby2' -- gives a diff answer than old code b/c of layer bug in old code
else
  JOBS=$*
fi

TESTSU_BASE=`pwd -P` #`dirname $0`
SOURCE_BASE="${TESTSU_BASE}/../src/"

# Run the job with the new version of the python code (Compatible with non-generalized md code)
#
PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/lsq-new-md-fmt.py" # Path to the python code.
RUN_LSQ_PYTHON_CODE="python $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header ff_groups.map TEST_SUITE_RUN"

# Run the job with the old version of the python code (Compatible with non-generalized md code)
#
#PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/lsq.py" # Path to the python code.
#RUN_LSQ_PYTHON_CODE="python $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header TEST_SUITE_RUN"

###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
make clean_lsq; make realclean_lsq
make house_lsq; rm -f ../test_suite-lsq/house_lsq; cp house_lsq ../test_suite-lsq/; make clean_lsq; make realclean_lsq;
cd ../test_suite-lsq

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

echo ""
echo "VALIDATING FOR SPLINES_LS..."

SET_PASSED=true
SVD_PASSED=true
ALL_PASSED=true

for i in $JOBS
do

	echo " "
	echo "Running $i test..."

	PASS=true

	cd $i
	rm -rf *diff*
	../house_lsq < fm_setup.in > fm_setup.out
	
	mv A.txt b.txt params.header fm_setup.out ff_groups.map current_output

	for j in A.txt b.txt params.header fm_setup.out
	do
		
		perl ../../contrib/compare/compare.pl correct_output/$j current_output/$j > diff-$j.out
	
		NO_DIFF_LINES=`cat diff-$j.out | wc -l`
	
		if [ $NO_DIFF_LINES -gt 0 ] ; then
			echo " "
			echo "Differences found in $j files... "
			echo " "
#		OLD WAY: 			
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
	
	cd ..
done

###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "VALIDATING FOR SVD SCRIPT..."

for i in $JOBS
do

	echo " "
	echo "Running $i test..."

	PASS=true
	TECHNICAL_PASS=true

	cd $i/current_output
	rm -rf diff-*
#	python $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header ff_groups.map TEST_SUITE_RUN > params.txt
	$RUN_LSQ_PYTHON_CODE > params.txt

#	cp params.txt force.txt current_output	
	
	for j in params.txt force.txt
	do
	
		if [ "$j" == params.txt ]; then
			j=params.txt-tailed
			tail -n+2 params.txt > $j
			tail -n+2 ../correct_output/params.txt > ../correct_output/$j
		fi
		
		# Ignore the date when running diff...
		
		diff ../correct_output/$j $j | awk '!/#/{print}' > ../diff-$j.out
	
		NO_DIFF_LINES=`cat ../diff-$j.out | wc -l`
	
		if [ $NO_DIFF_LINES -gt 0 ] ; then
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
done

echo " "

if   [ "$ALL_PASSED" = true ] ; then
	echo "ALL TESTS PASSED"
elif [ "$ALL_PASSED" = false ] ; then
	echo "AT LEAST ONE EACH OF SETUP AND SVD TEST(S) FAILED..."
	echo "Check individual results above to confirm if test had technical pass."
elif [ "$SVD_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SVD SCRIPT"
elif [ "$SET_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SETUP SCRIPT"
else
	echo "ERROR: BAD LOGIC IN TEST SUITE DRIVER (THIS SCRIPT)"
fi


	
	
	
	
exit 0
