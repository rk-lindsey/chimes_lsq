#!/bin/bash


###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

# Common function for test script initialization.
source ../src/bash/init_vars.sh
init_test_vars
echo "NP = $NP"

cd ..

SOURCE_BASE="${TESTSU_BASE}/../build/"
if [ ! -f $SOURCE_BASE/chimes_lsq ] ; then
    if ./install.sh  ; then
	echo "Compiling chimes_md succeeded"
    else
	echo "Compiling chimes_md failed"
	exit 1
    fi
fi    

cd -


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

	 if ! test_dir $i ; then
		  continue 
	 fi

	 PASS=true
	
	cd $i

	if $RUN_JOB ../../build/chimes_md run_md.in > run_md.out ; then
		 SUCCESS=1 
	else
		 echo "Chimes_md failed"
		 SUCCESS=0
		 PASS=false
		 ALL_PASS=false
	fi

	if [ ! -d current_output ] ; then mkdir current_output ; fi
	
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

	 cd ../test_suite-lsq 
	 ./run_test_suite.sh "$LSQ_FORCE_JOBS"
	 
	 cd -

	 echo " "
	 echo " ...Now running the force comparison tests... "
	 for i in $LSQ_FORCE_JOBS
	 do

		  if ! test_dir $i ; then
				continue 
		  fi
		  
		  PASS=true
		  
		  if [ ! -d ${TAG}${i} ] ; then
		  	echo "MD test directory ${TAG}${i} doesn't exist"
			continue
		  fi		  
		  
		  cd ${TAG}${i}
	
		# Grab the parameter and force files from the lsq test suite output
	
		  cp ../../test_suite-lsq/$i/current_output/params.txt    .
		  cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
		  cp ../../test_suite-lsq/$i/current_output/force.txt     .

		  if ../../build/chimes_md run_md.in > run_md.out ; then
				SUCCESS=1
		  else
				echo "Chimes_MD failed"

				../chimes_md run_md.in > run_md.out
				SUCCESS=0
				PASS=false
				ALL_PASS=false
		  fi
		  
		  if [ ! -d current_output ] ; then mkdir current_output; fi

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
	 if ! test_dir $i ; then
		  continue 
	 fi
	 cd $i ; RUN_JOB=$RUN_JOB  NP=${NP} make
	 echo "Testing $i"
	 if [ $? -ne 0 ] ; then
		  echo "		...Test failed."
		  ALL_PASS=false ;
	 else
		  echo "		...Test passed."
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
