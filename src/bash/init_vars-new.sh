#! /bin/bash

function init_test_vars
# Initialize variables for test suite runs.
{
    TESTSU_BASE=`pwd -P` #`dirname $0`

    SOURCE_BASE="${TESTSU_BASE}/../src/"

    PYTHON=python3 #/usr/tce/bin/python
 
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
		  stress-and-ener-4b
                  stress-and-ener-4b2
		  stress-and-ener-2b1
		  stress-and-ener-2b2
		  stress-and-ener-2b3
		  two-traj-files
		  one-file-three-trajs
		  nonorth2'


    LSQ_MAKE_JOBS='lsq2
                   tatb'

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
		 chebyfix
		 h2o-4bcheby-numstress
		 lj-stress
                 triclinic
                 ti-ortho
                 ti-nonortho'

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
                stress-and-ener-2b1
                stress-and-ener-4b
                stress-and-ener-4b2'
	 
	 NP=36
	 RUN_JOB="srun -n $NP"
	 RUN_JOB="srun -N 1 -n $NP -t 01:00:00 -ppdebug -A iap"
	 
	 if [ ! -v hosttype ] ; then
	     echo "No hosttype specified"
	     echo "Be sure to load modules/configure compilers by hand before running this script!"
	     NP=1
	     RUN_JOB=""
	 elif [[ "$hosttype" == "LLNL-LC" ]] ; then
	     source ${TESTSU_BASE}/../modfiles/LLNL-LC.mod
	     ICC=`which icc`	
	     MPI=`which mpicxx`    
	 elif [[ "$hosttype" == "UM-ARC" ]] ; then
	     source ${TESTSU_BASE}/../modfiles/UM-ARC.mod
	     ICC=`which icc`	
	     MPI=`which mpicxx`    
	 elif [[ "$hosttype" == "JHU-ARCH" ]] ; then
	     source ${TESTSU_BASE}/../modfiles/JHU-ARCH.mod
	     ICC=`which icc`
	     MPI=`which mpicxx`   
	 elif [[ "$hosttype" == "UT-TACC" ]] ; then
	     source ${TESTSU_BASE}/../modfiles/UT-TACC.mod
	 else
	     echo ""
	     echo "ERROR: Unknown hosttype ($hosttype) specified"
	     echo ""
	     echo "Valid options are:"
	     for i in `ls modfiles`; do echo "   ${i%.mod}"; done
	     echo ""
	     echo "Please run again with: export hosttype=<host type>; ./install.sh"
	     echo "Or manually load modules and run with: ./install.sh"
	     NP=1
	     RUN_JOB=""
	     exit 0
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
