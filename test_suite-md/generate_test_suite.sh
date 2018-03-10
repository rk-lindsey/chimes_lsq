#! /usr/bin/bash

NP=1

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
MD_TESTS[7]="h2o-4bcheby"

# Tests for compatibility between LSQ C++/python codes with the MD code
TAG="verify-lsq-forces-"

LSQ_TESTS[0]="chon-dftbpoly"	# -- DOESN'T EXIST IN ZCALC FOR MD!
LSQ_TESTS[1]="h2o-2bcheby"
LSQ_TESTS[2]="h2o-3bcheby"
LSQ_TESTS[3]="h2o-splines"
LSQ_TESTS[4]="h2o-4bcheby"
LSQ_TESTS[4]="h2o-invr" 	# -- DOESN'T EXIST IN ZCALC FOR MD!
LSQ_TESTS[5]="h2o-dftbpoly"	# -- DOESN'T EXIST IN ZCALC FOR MD!


# Iterate through the tests


echo " "
echo "SETTING UP FOR MD CODE..."
echo " "

cd ../src
rm -rf *o *dSYM chimes_md

if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
	 source /usr/local/tools/dotkit/init.sh
	 use ic-17.0.174
    use mvapich2-intel-2.2
else
    module load intel impi
fi


make -f Makefile-TS-MD chimes_md
rm -f ../test_suite-lsq/chimes_md;  mv chimes_md  ../test_suite-md/
cd ../test_suite-md

if [ $# -gt 0 ] ; then
	 MD_JOBS=$1
else
	 MD_JOBS="${MD_TESTS[@]}"
fi
echo "MD JOBS = $MD_JOBS"

if [ -n $LSQ_TESTS[0] ] ; then
	 LSQ_JOBS="${LSQ_TESTS[@]}"
fi

for i in "$MD_JOBS"
do

	echo " "
	echo "Running $i test..."
	
	cd $i
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		../chimes_md < run_md.in > run_md.out
			
	else
		srun -n $NP ../chimes_md < run_md.in > run_md.out
	fi		
		
	cp *.* current_output
	cp *.* correct_output

	
	cd ..
done


if [ -n "$LSQ_JOBS" ] ; then

echo " "
echo "SETTING UP FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../src
rm -rf *o *dSYM chimes_lsq chimes_md
cd ../test_suite-md

cd ../test_suite-lsq 
./run_test_suite.sh $LSQ_JOBS

cd ../src
cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in "${LSQ_JOBS}"
do

	echo " "
	echo "Running $i test..."
	
	cd ${TAG}${i}
	
	# Grab the parameter and force files from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt    .
	cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
	cp ../../test_suite-lsq/$i/current_output/force.txt     .
	
	../chimes_md < run_md.in > run_md.out		
	
	cp *.* current_output
	cp *.* correct_output
	
	cd ..
done	


echo " "

fi # -n LSQ_JOBS
	
exit 0
