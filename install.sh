#!/bin/bash

echo ""
echo "***"
echo "*** Executing chimes_lsq install.sh"
echo "***" 
echo ""

# Builds all relevant chimes_calculator executables/library files
# Run with:
# ./install.sh
# or
# ./install.sh <debug option (0 or 1)> <install prefix (full path)> <verbosity option (0 or 1 or 2 or 3)> <MPI option (0 or 1)>

DEBUG=${1-0}  # False (0) by default; if false, compiles with -O3, otherwise, uses -g
PREFX=${2-""} # Empty by default
VERBO=${3-1}  # Verbosity set to 1 by default, 0 gives minimal output, 3 gives minimal with DEBUG_CHEBY output, 4 gives all output
DOMPI=${4-1}  # Compile with MPI support by default


echo "Attempting to perform a fresh install"
echo "Imports directory will be deleted and re-cloned/installed. Proceed? (y/n)"
read PROCEED
if [[ "$PROCEED" == "n" ]] ; then
    echo 'Will use pre-existing imports directory'
else
    echo "Are you sure? Imports directory will be deleted! (y/n)"
    read PROCEED
    if [[ "$PROCEED" == "y" ]] ; then
        rm -rf imports
    fi
fi


# Determine computing environment and attempt to load module files automatically

lochost=`hostname`
hosttype=""

if [[ $lochost == *"arc-ts.umich.edu"* ]]; then
    hosttype=UM-ARC
elif [[ $lochost == *"quartz"* ]]; then
    hosttype=LLNL-LC
fi

# Load module files and configure compilers

ICC=`which g++` # Default option in case impi unavailable

if [[ "$hosttype" == "LLNL-LC" ]] ; then
    source modfiles/LLNL-LC.mod
    ICC=`which icc`
elif [[ "$hosttype" == "UM-ARC" ]] ; then
    source modfiles/UM-ARC.mod
    ICC=`which icc`
fi

echo "Detected hosttype $hosttype"
module list

MPI=`which mpicxx`


# Grab and install dependencies

if [[ ! -d imports ]] ; then
    ./clone-all.sh
fi

# Compile dlars if mpi compilers are available on a HPC platform

if [[ ! -z $hosttype ]] ; then
    cd contrib/dlars/src
    
    if [[ "$hosttype" == "LLNL-LC" ]] ; then
        make
    elif [[ "$hosttype" == "UM-ARC" ]] ; then
        make CXX=mpiicpc
    fi    

    cd - 1&>/dev/null
fi

# Compile molanal


cd contrib/molanal/src
make molanal.new
cd ../../..


# Clean up previous installation,

./uninstall.sh $PREFX


# Move into build directory

mkdir build
cd build


# Generate cmake flags

my_flags="-DCMAKE_CXX_COMPILER=${ICC}"

if [ ! -z $PREFX ] ; then
        my_flags="-DCMAKE_INSTALL_PREFIX=${PREFX}"
fi

if [ $DEBUG -eq 1 ] ;then
    my_flags="${my_flags} -Wall -DCMAKE_BUILD_TYPE=Debug"
else
    my_flags="${my_flags} -DCMAKE_BUILD_TYPE=Release"
fi

if   [ $VERBO -eq 0 ] ;then
        my_flags="${my_flags} -DVERBOSITY=0" 
    my_flags="${my_flags} -DDEBUG_CHEBY=0" 
elif [ $VERBO -eq 1 ] ; then
        my_flags="${my_flags} -DVERBOSITY=1" 
    my_flags="${my_flags} -DDEBUG_CHEBY=0" 
elif [ $VERBO -eq 2 ] ; then
        my_flags="${my_flags} -DVERBOSITY=0" 
    my_flags="${my_flags} -DDEBUG_CHEBY=1"     
elif [ $VERBO -eq 3 ] ; then
        my_flags="${my_flags} -DVERBOSITY=1" 
    my_flags="${my_flags} -DDEBUG_CHEBY=1"     
fi

if [ $DOMPI -eq 1 ] ;then
        my_flags="${my_flags} -DUSE_MPI=1" 
    my_flags="${my_flags} -DMPI_CXX_COMPILER=${MPI}"
else
        my_flags="${my_flags} -DUSE_MPI=0" 
fi

echo "compiling with flags: $my_flags"


# Setup, make and install

cmake $my_flags ..
make

cp ../src/chimes_lsq.py .
cp ../src/post_proc_chimes_lsq.py .

if [ $DOMPI -eq 1 ] ;then
 
 # Create some executables for the ALD

    cp chimes_md chimes_md-mpi
    cp chimes_lsq chimes_lsq.tmp

    my_flags=`echo $my_flags | awk '{for(i=1;i<=NF; i++){if($i~"DUSE_MPI=1"){$i="-DUSE_MPI=0"}}{print}}'`

    cmake $my_flags ..
    make

    cp chimes_md chimes_md-serial
    cp chimes_md-mpi chimes_md
    mv chimes_lsq.tmp chimes_lsq
fi

if [ ! -z $PREFX ] ; then
        make install
        cp src/chimes_lsq.py $PREFX
    cp src/post_proc_chimes_lsq.py $PREFX
fi

cd ..

echo ""
echo "***"
echo "*** chimes_lsq install complete"
echo "***"
echo ""      
