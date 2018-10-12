#!/bin/bash


###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

#module load intel impi

cd ../src
rm -rf *o *dSYM chimes_md

# Common function for test script initialization.
source ../src/bash/init_vars.sh
init_test_vars
echo "NP = $NP"

if make -f Makefile-TS-MD chimes_md ; then
	 echo "Succeeded in building chimes_md"
else
	 echo "Failed to build chimes_md"
	 exit
fi

rm -f ../test_suite-lsq/chimes_md;  mv chimes_md  ../test_suite-md/
cd ../test_suite-md


########################################
# Define tests within the test suite
########################################

## Allow command line arguments of jobs to test.  MD jobs should be single-quoted in a string followed 
## by LSQ jobs single quoted. (LEF)
##

if [ $# -gt 0 ] ; then
	 MD_JOBS=$1
	 LSQ_FORCE_JOBS=$2
	 MD_MAKE_JOBS=$3
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

if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
	 RUN_JOB=""
fi

for i in $MD_JOBS
do

	echo " "
	echo "Running $i test..."
	
	PASS=true
	
	cd $i

	if $RUN_JOB ../chimes_md < run_md.in > run_md.out ; then
		 SUCCESS=1 
	else
		 echo "Chimes_md failed"
		 SUCCESS=0
		 PASS=false
		 ALL_PASS=false
	fi
	cp *.* current_output

	if [[ $SUCCESS -eq 1 ]] ; then
	
		 for j in run_md.out traj.gen output.xyz 
		 do
			  if [[ -e current_output/$j  &&  -e correct_output/$j ]] ; then
					diff current_output/$j correct_output/$j > $j-diff.txt
					
					LINES=`wc -l $j-diff.txt | awk '{print $1}'`
					
					if [ $LINES -gt 0 ] ; then
						 echo " "
						 echo "		Differences found in $j files:"
						 echo " "
						 cat $j-diff.txt

						 PASS=false
						 ALL_PASS=false
					fi
			  fi
		 done
	fi
	
	if [ "$PASS" = true ] ; then
		echo "		...Test passed."
		rm -f ../diff-*
	else
		echo "		...Test failed."
	fi	

	
	cd ..
done


########################################
# Iterate through the tests -- MD/LSQ CODE COMPATIBILITY
########################################


if [ -n "$LSQ_FORCE_JOBS" ] ; then
	 echo " "
	 echo "VALIDATING FOR LSQ/MD CODE COMPATIBILITY..."
	 echo " "
	 echo " ...Beginning by running the lsq test suite... "

	 cd ../src
	 rm -rf *o *dSYM chimes_lsq
	 cd ../test_suite-md

	 cd ../test_suite-lsq 
	 ./run_test_suite.sh "$LSQ_FORCE_JOBS"
	 
	 cd ../test_suite-md

	 echo " "
	 echo " ...Now running the force comparison tests... "
	 for i in $LSQ_FORCE_JOBS
	 do

		  if [ -d "$i" ] ; then
				echo " "
				echo "Running $i test..."
		  else
				echo "$i directory was not found"
				continue
		  fi
		  
		  PASS=true
		  
		  cd ${TAG}${i}
	
	# Grab the parameter and force files from the lsq test suite output
	
		  cp ../../test_suite-lsq/$i/current_output/params.txt    .
		  cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
		  cp ../../test_suite-lsq/$i/current_output/force.txt     .

		  if ../chimes_md < run_md.in > run_md.out ; then
				SUCCESS=1
		  else
				echo "Chimes_MD failed"
				SUCCESS=0
				PASS=false
				ALL_PASS=false
		  fi

		  cp *.* current_output

		  if [[ $SUCCESS -eq 1 ]] ; then
				for j in run_md.out forceout-labeled.txt
				do
					 if [[ -e current_output/$j  &&  -e correct_output/$j ]] ; then
						  diff current_output/$j correct_output/$j > $j-diff.txt
						  
						  LINES=`wc -l $j-diff.txt | awk '{print $1}'`
						  
						  if [ $LINES -gt 0 ] ; then
								echo " "
								echo "		Differences found in $j files:"
								echo " "

								cat $j-diff.txt

								PASS=false
								ALL_PASS=false
						  fi
					 fi
				done	
		  fi
	
		  if [ "$PASS" = true ] ; then
				echo "		...Test passed."
				rm -f ../diff-*
		  else
				echo "		...Test failed."
		  fi	
		  
		  
		  cd ..
	 done
fi

echo "PERFORMING MAKEFILE TESTS"
for i in $MD_MAKE_JOBS
do
	 cd $i ; make NP=${NP}
	 echo "Testing $i"
	 if [ $? -ne 0 ] ; then
		  echo "Test $i failed"
		  ALL_PASS=false ;
	 fi
	 echo " "
	 cd ../
done

if   [ "$ALL_PASS" = true ] ; then
	echo "ALL TESTS PASSED"
elif [ "$ALL_PASS" = false ] ; then
	echo "AT LEAST ONE EACH OF MD and MD/LSQ COMPATIBILITY TEST(S) FAILED"
else
	echo "ERROR: BAD LOGIC IN TEST SUITE DRIVER (THIS SCRIPT)"
fi
	
exit 0
