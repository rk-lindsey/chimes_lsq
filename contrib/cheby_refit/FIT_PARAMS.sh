#!/bin/bash

# Example of  a shell script that could be used to run a fit/refit 


CODE_BASE="/Users/lindsey11/Desktop/FORCE_MATCHING_VERSIONS/CURRENT_WORKING/"

# Do the C++/Python LSQ fitting as usual (rmin = 1.0)

echo "Running lsq C++"

${CODE_BASE}/src/house_lsq < fm_setup.in > fm_setup.out

echo "Running lsq Python"

${CODE_BASE}/src/lsq-new-md-fmt.py A.txt b.txt params.header ff_groups.map > params.txt

# Get the pes scan for the FF  (Note we only have one pair type here)

echo "Getting scan"

${CODE_BASE}/src/house_md < pes_scan.in > /dev/null

echo "New parameters (with rmin = 0) are:"

python ${CODE_BASE}/contrib/cheby_refit/cheby_fitting.py 5 1.0 3.15 1.25 0.25 0.02 2b_Cheby_Pot-CC.dat 2 True | tee refit.txt



