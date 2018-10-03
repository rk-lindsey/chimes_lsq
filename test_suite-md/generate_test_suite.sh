#! /usr/bin/bash

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

cd ../src
rm -rf *o *dSYM chimes_md


# Common function for test script initialization.
source ../src/bash/init_vars.sh
init_test_vars
echo "NP = $NP"

make -f Makefile-TS-MD chimes_md
rm -f ../test_suite-lsq/chimes_md;  mv chimes_md  ../test_suite-md/
cd ../test_suite-md

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
	
	cd $i
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		 if ../chimes_md < run_md.in > run_md.out ; then
			  SUCCESS=1
		 else
			  echo "Chimes_md failed"
			  SUCCESS=0
		 fi
			
	else
		 if $RUN_JOB ../chimes_md < run_md.in > run_md.out ; then
			  SUCCESS=1
		 else
			  echo "Chimes_md failed"
			  SUCCESS=0
		 fi
			  
	fi		
	cp *.* current_output
	if [[ $SUCCESS -eq 1 ]] ; then
		 cp *.* correct_output
	fi
	
	cd ..
done


if [ -n "$LSQ_FORCE_JOBS" ] ; then

echo " "
echo "SETTING UP FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../src
rm -rf *o *dSYM chimes_lsq chimes_md
cd ../test_suite-md

cd ../test_suite-lsq 
./run_test_suite.sh $LSQ_FORCE_JOBS

cd ../src
cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in ${LSQ_FORCE_JOBS}
do

	echo " "
	echo "Running $i test..."
	
	cd ${TAG}${i}
	
	# Grab the parameter and force files from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt    .
	cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
	cp ../../test_suite-lsq/$i/current_output/force.txt     .
	
	if ../chimes_md < run_md.in > run_md.out ; then
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
	
	cd ..
done	


echo " "

fi # -n LSQ_JOBS
	
exit 0
