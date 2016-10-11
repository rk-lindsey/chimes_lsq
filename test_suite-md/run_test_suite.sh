#!/bin/bash

MPICOMMAND=mpirun
#MPICOMMAND="${MPICOMMAND} -stdin all"	# Only needed if compiled with mpicxx (OSX)

###############################################################
#
# Make a fresh compilation of the code
#
###############################################################
cd ../src
rm -rf *o *dSYM
#make house_lsq; rm -f ../test_suite-lsq/house_lsq; cp house_lsq ../test_suite-lsq/
make house_md;  rm -f ../test_suite-lsq/house_md;  cp house_md  ../test_suite-md/
cd ../test_suite-md

########################################
# Define tests within the test suite
########################################

# Tests specifically for the MD code

MD_TESTS[0]="h2o-2bcheby"
MD_TESTS[1]="h2o-3bcheby" 
MD_TESTS[2]="h2o-splines"
MD_TESTS[3]="generic-lj"
MD_TESTS[4]="h2o-2bcheby-genvel" 
MD_TESTS[5]="h2o-2bcheby-numpress"
MD_TESTS[6]="h2o-2bcheby-velscale"

# Now that handling of layers has changed, we no longer expect to 
# recover the same forces on atoms that were observed in the LSQ step
# For this reason, these tests have been omitted. 
#
# Eventually, I'll add a test for a system where the cutoffs are
# within 0.5 x the natural box length.

#LSQ_TESTS[0]="chon-dftbpoly"
#LSQ_TESTS[1]="h2o-2bcheby"
#LSQ_TESTS[2]="h2o-3bcheby"
#LSQ_TESTS[3]="h2o-splines"
#LSQ_TESTS[4]="h2o-invr"
#LSQ_TESTS[5]="h2o-dftbpoly"

NP=0

## Allow command line arguments of jobs to test.  MD jobs should be single-quoted in a string followed 
## by LSQ jobs single quoted. (LEF)
##

	MD_JOBS="${MD_TESTS[@]}"
	LSQ_JOBS="${LSQ_TESTS[@]}"

if [ $# -gt 0 ] ; then

	if   [ "$1" == "NP" ] ; then
		NP=$2
	elif [ "$3" == "NP" ] ; then
		MD_JOBS=$1
		LSQ_JOBS=$2
		NP=$4
	else
		MD_JOBS=$1
		LSQ_JOBS=$2
	fi
fi


# Tests for compatibility between LSQ C++/python codes with the MD code
TAG="verify-lsq-forces-"

########################################
# Iterate through the tests -- MD CODE
########################################

echo " "
echo "VALIDATING FOR MD CODE..."
echo " "

ALL_PASS=true

for i in $MD_JOBS
do

	echo " "
	echo "Running $i test..."
	
	PASS=true
	
	cd $i

	if [ $NP -eq 0 ] ; then
		$MPICOMMAND ../house_md < run_md.in > run_md.out
			
	else
		$MPICOMMAND -np $NP ../house_md < run_md.in > run_md.out
		
	fi
	
	cp *.* current_output

	
	for j in run_md.out traj.gen output.xyz 
	do
		if [[ -e current_output/$j  &&  -e correct_output/$j ]] ; then
			diff current_output/$j correct_output/$j > $j-diff.txt
			
			LINES=`wc -l $j-diff.txt | awk '{print $1}'`
			
			if [ $LINES -gt 0 ] ; then
				echo " "
				echo "		Differences found in $j files:"
				echo " "
				
				PASS=false
				ALL_PASS=false
			fi
		fi
	
		
	
	done
	
	if [ "$PASS" = true ] ; then
		echo "		...Test passed."
		rm -f ../diff-*
	else
		echo "		...Test failed."
	fi	

	
	cd ..
done

exit 0

########################################
# Iterate through the tests -- MD/LSQ CODE COMPATIBILITY
########################################


echo " "
echo "VALIDATING FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../test_suite-lsq 
./run_test_suite.sh $LSQ_JOBS
cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in $LSQ_JOBS
do

	echo " "
	echo "Running $i test..."

	PASS=true
	
	cd ${TAG}${i}
	
	# Grab the parameter file from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt .
	
	if [ $NP -eq 0 ] ; then
		../house_md < run_md.in > run_md.out		
	else
		$MPICOMMAND -np $NP ../house_md < run_md.in > run_md.out
	fi	

	cp *.* current_output

	
	for j in run_md.out traj.gen output.xyz forceout.txt
	do
		if [[ -e current_output/$j  &&  -e correct_output/$j ]] ; then
			diff current_output/$j correct_output/$j > $j-diff.txt
			
			LINES=`wc -l $j-diff.txt | awk '{print $1}'`
			
			if [ $LINES -gt 0 ] ; then
				echo " "
				echo "		Differences found in $j files:"
				echo " "
				
				PASS=false
				ALL_PASS=false
			fi
		fi

	done	
	
	if [ "$PASS" = true ] ; then
		echo "		...Test passed."
		rm -f ../diff-*
	else
		echo "		...Test failed."
	fi	

	
	cd ..
done	


echo " "

if   [ "$ALL_PASS" = true ] ; then
	echo "ALL TESTS PASSED"
elif [ "$ALL_PASS" = false ] ; then
	echo "AT LEAST ONE EACH OF MD and MD/LSQ COMPATIBILITY TEST(S) FAILED"
else
	echo "ERROR: BAD LOGIC IN TEST SUITE DRIVER (THIS SCRIPT)"
fi
	
exit 0
