These are test problems for force matching MD.

The test problems are as follows:

gather-force: Collect forces from OpenMX DFT code into a single
input.xyzf file for force matching.

h2o-lsq: Perform least-squares fitting of the forces in input.xyzf, using both
	 spline and Chebyshev polynomial functions.

h2o-md: Perform molecular dynamics using potential parameters derived
from force matching.  Both the spline and Chebyshev potentials are
tested.

h2o-nonlin: Perform nonlinear optimization of the overcoordination
parameters used in force matching.

lsq2: Perform least-squares fitting of a larger set of forces.

openmx: Generate a series of OpenMX input files for MD snapshots.  The
MD snapshots are taken from a DFTB-format .gen file.  The .gen
file is used as a simple lowest common denominator trajectory file.
The configurations could actually come from force-matched MD, DFT, etc.

test-force: Test the consistency of forces output by the lsq.py
least-squares fitting script with those calculated by spline_md.C.

