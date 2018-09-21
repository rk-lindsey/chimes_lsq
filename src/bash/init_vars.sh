#! /bin/bash

function init_test_vars
# Initialize variables for test suite runs.
{
	 TESTSU_BASE=`pwd -P` #`dirname $0`
	 SOURCE_BASE="${TESTSU_BASE}/../src/"
	 # Intel parallel python - supports thread parallelism.
	 #PYTHON=/collab/usr/global/tools/intel/chaos_5_x86_64_ib/python-2.7.10/bin/python
	 # Default python
	 PYTHON=python
	 # Run the job with the new version of the python code (Compatible with non-generalized md code)
	 #
	 PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/lsq.py" # Path to the python code.
	 EPS_FAC=1.0e-5 # 1.0E-5 is the old default value... should match value used in gen test suite script.  1.0e-09
	 # SVD regularization factor.
	 RUN_LSQ_PYTHON_CODE="$PYTHON $PATH_TO_LSQ_PY_CODE A.txt b.txt params.header ff_groups.map ${EPS_FAC} TEST_SUITE_RUN"

	 LSQ_ALL_JOBS='h2o-splines h2o-invr h2o-dftbpoly chon-dftbpoly h2o-2bcheby h2o-3bcheby h2o-2bcheby2 h2o-3bcheby2 h2o-3bcheby3 par-ewald h2o-4bcheby test_4atoms test_4atoms.2' 

	 MD_JOBS='h2o-2bcheby h2o-3bcheby h2o-4bcheby h2o-splines generic-lj h2o-2bcheby-genvel h2o-2bcheby-numpress h2o-3bcheby-numpress h2o-4bcheby-numpress h2o-2bcheby-velscale h2o-3bcheby3 '

	 MD_MAKE_JOBS='verify-invert verify-translate verify-scramble h2o-4bcheby-numforce verify-relabel verify-relabel.2'


	 LSQ_FORCE_JOBS='h2o-2bcheby h2o-3bcheby h2o-splines h2o-4bcheby'

	 if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
		  source /usr/local/tools/dotkit/init.sh
		  use ic-17.0.174
		  use mvapich2-intel-2.2
		  NP=24
		  RUN_JOB="srun -n $NP"
	 elif [ "$SYS_TYPE" == "toss_3_x86_64_ib" ] ; then
		  module load impi/2018.0
		  module load intel/18.0.1
		  NP=36
		  RUN_JOB="srun -n $NP"
	 else
		  NP=1
		  RUN_JOB=""
	 fi
	 # Number of threads for SVD decomposition
	 NUM_THREADS=$NP
	 export OMP_NUM_THREADS=$NUM_THREADS	
}
