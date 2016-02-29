Test directory to do force matching of ninth-order polynomial to determine DFTB Erep. 
Nir Goldman, 2/27/2016
==
STEP 1: Determine needed values of forces from Erep. This is done by 
first running DFTB for a series of configurations with the 
contribution from Erep zeroed out. These forces are then subtracted 
from the user supplied forces from DFT to yield the fitting data set
for force matching. In this example, the configurations and forces 
are supplied in input_openmx_1.5gcc.xyzf and the forces to be 
fit to are saved in input_dftb.xyzf.

STEP 2: Perform force matching by computing the optimal polynomial
coefficients through SVD, as in other examples. PLEASE NOTE: the
DFTB forces are left unit atomic units, not kcal/mol A as in other 
examples. As a result, the RMS error should be multiplied by the 
appropriate unit conversion factor in order to compare to other 
force matching results. 

The makefile command is also set up to perform a check on the forces 
from the SVD calculation by computing the Erep forces with the 
resulting coefficients in a separate code. 

This step will also create the corresponding O-O.skf, O-H.skf, 
H-O.skf, and H-H.skf files for future calculations with DFTB+.
