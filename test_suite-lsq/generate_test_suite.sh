# NOTE: The path below needs to point to the "lsq-new-md-fmt.py" version of the lsq code.
#       If you want to use other versions, you'll probably need to modify the inputs that
#       are sent to the script (way down below)

###############################################################
#
# Determine the location of necessary files
#
###############################################################

TESTSU_BASE=`pwd -P` #`dirname $0`
SOURCE_BASE="${TESTSU_BASE}/../src/"

# Run the job with the new version of the python code (Compatible with non-generalized md code)
#
PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/lsq-new-md-fmt.py" # Path to the python code.

###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
make house_lsq; cp house_lsq ../test_suite-lsq/; make clean_lsq; make realclean_lsq;
cd ../test_suite-lsq

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

echo ""
echo "SETTING UP FOR SPLINES_LS..."

for i in h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby
do

	echo " "
	echo "Running $i test..."

	cd $i
	
	../house_lsq < fm_setup.in > fm_setup.out
 	cp A.txt b.txt params.header fm_setup.out ff_groups.map correct_output	
	mv A.txt b.txt params.header fm_setup.out ff_groups.map current_output
	
	cd ..
done


###############################################################
#
#  Run the tests for the SVD script... 
#
###############################################################

echo "SETTING UP FOR SVD SCRIPT..."

for i in  h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby 
do
	echo " "
	echo "Running $i test..."

	cd $i/current_output
	
	python $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header ff_groups.map TEST_SUITE_RUN > params.txt
	
	
	for j in params.txt force.txt
	do
	
		if [ "$j" == params.txt ]; then
			j=params.txt-tailed
			tail -n+2 params.txt > $j
			tail -n+2 ../correct_output/params.txt > ../correct_output/$j
		fi
		
	done
	

	cd ../..
	
	cp $i/current_output/* $i/correct_output
done

exit 0
