#! /usr/bin/python
#
# Compares the force between two files.  Used to compare forceout.txt from chimes_md to the
# force.txt file produced by lsq-new-md-fmt.py
#
import sys
import math

af = open(sys.argv[1], "r").readlines()
bf = open(sys.argv[2], "r").readlines()

err = 0

nf = min( len(af), len(bf) )

print "Number of force components = ", nf

for i in range (0,nf):
    aff = float(af[i])
    bff = float(bf[i])
    err += (aff - bff) * (aff - bff)

err /= nf
err = math.sqrt(err)

print "RMS Error = %7.5f kcal/mol-Ang." % err

