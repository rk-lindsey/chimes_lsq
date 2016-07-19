PATH_TO_LSQ_PY_CODE="/Users/lindsey11/Desktop/FORCE_MATCHING_VERSIONS/BECKY-LSQ-PY-REG-071816/lsq.py" # Path to the python code.


SET_PASSED=true
SVD_PASSED=true
ALL_PASSED=true

###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
make; cp splines_ls ../test_suite/; make clean; make realclean;
cd ../test_suite

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

echo ""
echo "VALIDATING FOR SPLINES_LS..."

for i in h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby
do

	echo " "
	echo "Running $i test..."

	PASS=true

	cd $i
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	mv A.txt b.txt params.header fm_setup.out current_output

	for i in A.txt b.txt params.header fm_setup.out
	do
		diff correct_output/$i current_output/$i > diff-$i.out
	
		NO_DIFF_LINES=`cat diff-$i.out | wc -l`
	
		if [ $NO_DIFF_LINES -gt 0 ] ; then
			echo " "
			echo "Differences found in $i files:"
			echo " "
			cat diff-$i.out
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

echo " "

###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "VALIDATING FOR SVD SCRIPT..."

for i in h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby 
do

	echo " "
	echo "Running $i test..."

	PASS=true

	cd $i/current_output
	rm -rf diff-*
	python $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header > params.txt
#	mv params.txt force.txt current_output
	
	for i in params.txt force.txt
	do
		if [ "$i" == params.txt ]; then
			i=params.txt-tailed
			tail -n+2 params.txt > $i
			tail -n+2 ../correct_output/params.txt > ../correct_output/$i
		fi
		
		# Ignore the date when running diff...
		
		diff ../correct_output/$i $i | awk '!/#/{print}' > ../diff-$i.out
	
		NO_DIFF_LINES=`cat ../diff-$i.out | wc -l`
	
		if [ $NO_DIFF_LINES -gt 0 ] ; then
			echo " "
			echo "Differences found in $i files:"
			echo " "
			cat ../diff-$i.out
			echo " "
			PASS=false
		fi
		
		if [ "$i" == params.txt ]; then
			rm -rf $i
			rm -rf ../correct_output/$i
		fi
	done
	
	if [ "$PASS" = true ] ; then
		echo "		...Test passed."
		rm -f ../diff-*
	else
		SVD_PASSED=false
		ALL_PASSED=false
		echo "		...Test failed."
	fi

	cd ../..
done

echo " "

if   [ "$ALL_PASSED" = true ] ; then
	echo "ALL TESTS PASSED"
elif [ "$ALL_PASSED" = false ] ; then
	echo "AT LEAST ONE EACH OF SETUP AND SVD TEST(S) FAILED"
elif [ "$SVD_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SVD SCRIPT"
elif [ "$SET_PASSED" = false ] ; then
	echo "TEST(S) FAILED FOR SETUP SCRIPT"
else
	echo "ERROR: BAD LOGIC IN TEST SUITE DRIVER (THIS SCRIPT)"
fi


	
	
	
	
exit 0


###############################################################
#
# Make a fresh compilation of the code... this is used to 
#  build the test suite files... eventually move this to its own
#  script...
#
###############################################################

cd ../src
make; cp splines_ls ../test_suite/; make clean; make realclean;
cd ../test_suite

###############################################################
#
#  Run test no. X: CHON atoms with DFTBPOLY pair potential
#
###############################################################

	echo " "
	echo "Running chon-dftbpoly test..."

	PASS=false

	cd chon-dftbpoly
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	python ../../../BECKY-LSQ-PY-REG-070516/lsq.py A.txt b.txt params.header > params.txt
	mv A.txt b.txt params.header fm_setup.out force.txt params.txt correct_output

	cd ..
	
###############################################################
#
#  Run test no. X: H2O with CHEBYSHEV pair potential (2-Body)
#
###############################################################

	echo " "
	echo "Running h2o-2bcheby test..."

	PASS=false

	cd h2o-2bcheby
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	python ../../../BECKY-LSQ-PY-REG-070516/lsq.py A.txt b.txt params.header > params.txt
	mv A.txt b.txt params.header fm_setup.out force.txt params.txt correct_output

	cd ..

###############################################################
#
#  Run test no. X: H2O with INVRSE_R pair potential
#
###############################################################

	echo " "
	echo "Running h2o-invr test..."

	PASS=false

	cd h2o-invr
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	python ../../../BECKY-LSQ-PY-REG-070516/lsq.py A.txt b.txt params.header > params.txt
	mv A.txt b.txt params.header fm_setup.out force.txt params.txt correct_output

	cd ..
	
###############################################################
#
#  Run test no. X: H2O with DFTBPOLY pair potential
#
###############################################################

	echo " "
	echo "Running h2o-dftbpoly test..."

	PASS=false

	cd h2o-dftbpoly
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	python ../../../BECKY-LSQ-PY-REG-070516/lsq.py A.txt b.txt params.header > params.txt
	mv A.txt b.txt params.header fm_setup.out force.txt params.txt correct_output

	cd ..	
	
	
	
	
###############################################################
#
#  Run test no. X: H2O with DFTBPOLY pair potential
#
###############################################################

	echo " "
	echo "Running h2o-splines test..."

	PASS=false

	cd h2o-splines
	rm -rf diff-*
	../splines_ls < fm_setup.in > fm_setup.out
	python ../../../BECKY-LSQ-PY-REG-070516/lsq.py A.txt b.txt params.header > params.txt
	mv A.txt b.txt params.header fm_setup.out force.txt params.txt correct_output

	cd ..	
		
	
	

