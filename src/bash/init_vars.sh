#! /bin/bash

function init_test_vars
# Initialize variables for test suite runs.
{
    TESTSU_BASE=`pwd -P` #`dirname $0`

    SOURCE_BASE="${TESTSU_BASE}/../src/"
    # Intel parallel python - supports thread parallelism.
    #PYTHON=/collab/usr/global/tools/intel/chaos_5_x86_64_ib/python-2.7.10/bin/python
    # Default python

    PYTHON=python2.7 #/usr/tce/bin/python
    # Run the job with the new version of the python code (Compatible with non-generalized md code)
    #

    PATH_TO_LSQ_PY_CODE="${SOURCE_BASE}/chimes_lsq.py" # Path to the python code.

    # SVD regularization factor.
    EPS_FAC=1.0e-5 # 1.0E-5 is the old default value... should match value used in gen test suite script.  1.0e-09

    RUN_LSQ_PYTHON_CODE="$PYTHON $PATH_TO_LSQ_PY_CODE --eps ${EPS_FAC} --test_suite true"

    LSQ_ALL_JOBS='h2o-2bcheby  
		  h2o-3bcheby 
		  h2o-3bcheby2 
		  h2o-3bcheby3 
		  h2o-4bcheby
		  h2o-4bcheby-small
		  h2o-invr
		  fitstrs
		  fitstrs3b
		  test_4atoms 
		  test_4atoms.2 
		  special3b 
		  special4b
                  split_files
		  stress-and-ener-4b
		  stress-and-ener-2b1
		  stress-and-ener-2b2
		  stress-and-ener-2b3
		  two-traj-files
		  one-file-three-trajs
		  nonorth2'


    LSQ_MAKE_JOBS='lsq2'

	MD_JOBS='carbon-penalty
		 h2o-2bcheby
		 h2o-2bcheby-genvel
		 h2o-3bcheby 
		 h2o-4bcheby 
		 hn3-published
		 generic-lj
		 h2o-3bcheby-numpress 
		 h2o-4bcheby-numpress 
		 h2o-3bcheby3 
		 nonorth
		 npt-lj
		 npt-berend
		 npt-berend-aniso
		 nvt-berend
		 serial-chimes
		 serial-chimes-3b2
		 serial-chimes4b
		 small-lj
		 special3b 
		 special4b
		 chebyfix'

    MD_MAKE_JOBS='verify-invert 
		verify-translate 
		verify-scramble 
		h2o-4bcheby-numforce 
		verify-relabel 
		verify-relabel.2'

	LSQ_FORCE_JOBS='h2o-3bcheby 
		h2o-4bcheby 
		special3b 
		special4b
                stress-and-ener-2b1'

	 if [ "$SYS_TYPE" == "chaos_5_x86_64_ib" ] ; then
	     source /usr/local/tools/dotkit/init.sh
	     use ic-17.0.174
	     use mvapich2-intel-2.2
	     NP=24
	     RUN_JOB="srun -n $NP"
	 elif [ "$SYS_TYPE" == "toss_3_x86_64_ib" ] ; then
	     module load cmake/3.14.5
	     module load impi/2018.0
	     module load intel/18.0.1
	     module load python/2.7.16
	     NP=36
	     RUN_JOB="srun -n $NP"
	 elif [ "$SYS_TYPE" == "toss_3_x86_64" ] ; then
	     module load cmake/3.14.5	     
	     module load impi/2018.0
	     module load intel/18.0.1
	     module load python/2.7.16		  
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

function test_dir
# Test to see that a directory exists, and print an informational message.
{
	if [[ -d "$1" ]] ; then
		 echo " "
		 echo "Running $1 test..."
		 return 0
	else
		 echo " "
		 echo "$1 directory was not found"
		 return 1
	fi
}
