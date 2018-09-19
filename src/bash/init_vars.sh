#! /bin/bash

function init_test_vars
# Initialize variables for test suite runs.
{
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
}
