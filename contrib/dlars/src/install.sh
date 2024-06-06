#! /bin/bash
# Compile the dlars application on LC TOSS3 computer systems.
module load intel/19.0.4
module load mkl/2020.0

DEBUG=${1-0}
VERBOSE=${2-0}

make clean  2>&1 /dev/null

if [[ $DEBUG == '1' ]] ; then
    make debug
else
    make opt
fi
