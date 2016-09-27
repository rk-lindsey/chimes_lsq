###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
make clean_lsq; make realclean_lsq
make house_lsq; rm -f ../test_suite-lsq/house_lsq; cp house_lsq ../test_suite-lsq/; make clean_lsq; make realclean_lsq;
make clean_md; make realclean_md
make house_md; rm -f ../test_suite-lsq/house_md; cp house_md ../test_suite-md/; make clean_md; make realclean_md;
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

# Tests for compatibility between LSQ C++/python codes with the MD code
TAG="verify-lsq-forces-"

LSQ_TESTS[0]="chon-dftbpoly"
LSQ_TESTS[1]="h2o-2bcheby"
LSQ_TESTS[2]="h2o-3bcheby"
LSQ_TESTS[3]="h2o-splines"
LSQ_TESTS[4]="h2o-invr"
LSQ_TESTS[5]="h2o-dftbpoly"

########################################
# Iterate through the tests -- MD CODE
########################################

echo " "
echo "VALIDATING FOR MD CODE..."
echo " "

ALL_PASS=true

for i in "${MD_TESTS[@]}"
do

	echo " "
	echo "Running $i test..."
	
	PASS=true
	
	cd $i
	../house_md < run_md.in > run_md.out		
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

########################################
# Iterate through the tests -- MD/LSQ CODE COMPATIBILITY
########################################


echo " "
echo "VALIDATING FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../test_suite-lsq 
./run_test_suite.sh
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
	
	../house_md < run_md.in > run_md.out	

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
