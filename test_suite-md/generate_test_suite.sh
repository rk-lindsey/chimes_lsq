

NP=16

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

#LSQ_TESTS[0]="chon-dftbpoly"	# -- DOESN'T EXIST IN ZCALC FOR MD!
LSQ_TESTS[1]="h2o-2bcheby"
LSQ_TESTS[2]="h2o-3bcheby"
LSQ_TESTS[3]="h2o-splines"
#LSQ_TESTS[4]="h2o-invr" 	# -- DOESN'T EXIST IN ZCALC FOR MD!
#LSQ_TESTS[5]="h2o-dftbpoly"	# -- DOESN'T EXIST IN ZCALC FOR MD!


# Iterate through the tests


echo " "
echo "SETTING UP FOR MD CODE..."
echo " "

cd ../src
rm -rf *o *dSYM house_md
cp Makefile Makefile-back

cp Makefile-TS-LSQ Makefile
module load intel
make house_md
mv house_md ../test_suite-md/house_md-serial

cp Makefile-TS-MD Makefile
module load intel impi
make house_md;  
rm -f ../test_suite-lsq/house_md;  mv house_md  ../test_suite-md/
mv Makefile-back Makefile
cd ../test_suite-md


for i in "${MD_TESTS[@]}"
do

	echo " "
	echo "Running $i test..."
	
	cd $i
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		../house_md < run_md.in > run_md.out
			
	else

		if [[ "$i" == "generic-lj" || "$i" == "h2o-2bcheby-numpress" ]] ; then
		
			# Numerical pressure calcs currently only supported on 1 proc
			echo "Running in serial (numerical pressure)"
			../house_md-serial < run_md.in > run_md.out
			
		else
			srun -n $NP ../house_md < run_md.in > run_md.out
		fi
	fi		
		
	cp *.* current_output
	cp *.* correct_output

	
	cd ..
done


echo " "
echo "SETTING UP FOR LSQ/MD CODE COMPATIBILITY..."
echo " "
echo " ...Beginning by running the lsq test suite... "

cd ../src
rm -rf *o *dSYM house_lsq house_md
cp Makefile Makefile-back
cp Makefile-TS-MD-Verif Makefile
cd ../test_suite-md

cd ../test_suite-lsq 
./run_test_suite.sh $LSQ_JOBS

cd ../src
mv Makefile-back Makefile

cd ../test_suite-md

echo " "
echo " ...Now running the force comparison tests... "
for i in "${LSQ_TESTS[@]}"
do

	echo " "
	echo "Running $i test..."
	
	cd ${TAG}${i}
	
	# Grab the parameter and force files from the lsq test suite output
	
	cp ../../test_suite-lsq/$i/current_output/params.txt    .
	cp ../../test_suite-lsq/$i/current_output/ff_groups.map . 
	cp ../../test_suite-lsq/$i/current_output/force.txt     .
	
	../house_md < run_md.in > run_md.out		
	
	cp *.* current_output
	cp *.* correct_output
	
	cd ..
done	


echo " "

	
exit 0
