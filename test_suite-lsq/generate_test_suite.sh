# NOTE: The path below needs to point to the "lsq-new-md-fmt.py" version of the lsq code.
#       If you want to use other versions, you'll probably need to modify the inputs that
#       are sent to the script (way down below)
#
# NOTE: Make sure the right compiler is specified in the Makefile... don't use MPI for lsq.
#

###############################################################
#
# Determine the location of necessary files
#
###############################################################


if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
	 source /usr/local/tools/dotkit/init.sh
	 use ic-17.0.174
    use mvapich2-intel-2.2
else
    module load intel impi
fi

# Number of MPI tasks.
NP=16

# Number of threads for SVD decomposition
NUM_THREADS=16

# SVD regularization factor.

EPS_FAC=1.0e-5 # 1.0E-5 is the old default value... should match value used in gen test suite script.  1.0e-09

TESTSU_BASE=`pwd -P` #`dirname $0`
SOURCE_BASE="${TESTSU_BASE}/../src/"

# Intel parallel python - supports thread parallelism.
PYTHON=/collab/usr/global/tools/intel/chaos_5_x86_64_ib/python-2.7.10/bin/python
# Default python
# PYTHON=python

# Run the job with the new version of the python code (Compatible with non-generalized md code)
#
PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/lsq-new-md-fmt.py" # Path to the python code.

###############################################################
#
# Make a fresh compilation of the code
#
###############################################################

cd ../src
rm -f *o chimes_lsq
make -f Makefile-TS-LSQ chimes_lsq
mv chimes_lsq ../test_suite-lsq/
cd ../test_suite-lsq

###############################################################
#
#  Run the tests for the splines_ls program
#
###############################################################

if [ $# -eq 0 ] 
# Use default JOBS.  
then
#  h2o-3bcheby2' -- gives a diff answer than old code b/c of layer bug in old code
JOBS='h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby h2o-2bcheby2 h2o-3bcheby2 h2o-3bcheby3 par-ewald' 
else
# Take JOBS from command line.
  JOBS=$*
fi

echo ""
echo "SETTING UP FOR SPLINES_LS..."

for i in $JOBS
do

	echo " "
	echo "Running $i test..."

	cd $i
	
	if [[ $NP -eq 0 || $NP -eq 1 ]] ; then
		 ../chimes_lsq < fm_setup.in > fm_setup.out
	else
		 srun -n $NP ../chimes_lsq < fm_setup.in > fm_setup.out
	fi	
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

for i in  $JOBS
do
	echo " "
	echo "Running $i test..."

	cd $i/current_output
	export OMP_NUM_THREADS=$NUM_THREADS	
	$PYTHON $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header ff_groups.map $EPS_FAC TEST_SUITE_RUN > params.txt
	
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
