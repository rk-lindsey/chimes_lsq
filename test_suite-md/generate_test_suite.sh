

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

# Tests for compatibility between LSQ C++/python codes with the MD code
TAG="verify-lsq-forces-"

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


# Iterate through the tests


echo " "
echo "SETTING UP FOR MD CODE..."
echo " "

ALL_PASS=true

for i in "${MD_TESTS[@]}"
do

	echo " "
	echo "Running $i test..."
	
	PASS=true
	
	cd $i
	mpirun -np 1 ../house_md < run_md.in > run_md.out		
	cp *.* current_output
	cp *.* correct_output

	
	cd ..
done


exit 0

echo " "
echo "SETTING UP FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../test_suite-lsq 
./generate_test_suite.sh
cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in "${LSQ_TESTS[@]}"
do

	echo " "
	echo "Running $i test..."

	PASS=true
	
	cd ${TAG}${i}
	
	# Grab the parameter file from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt .
	
	mpirun -np 1 ../house_md < run_md.in > run_md.out	

	cp *.* current_output
	cp *.* correct_output
	
	cd ..
done	


echo " "

	
exit 0
