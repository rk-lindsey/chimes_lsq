#!/bin/bash

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


# Initialize test suite parameters
source ../src/bash/init_vars.sh
init_test_vars
echo "Number of processors = $NP"
echo "Warning: this script will overwrite all reference test output."
echo "Continue ? (yes/no)"

read ok_run
if test "x$ok_run" != "xyes" ; then
	echo 'Quitting'
	exit 0
fi


# Compile the code

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
