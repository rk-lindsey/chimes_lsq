# Takes as input a number of frames and a .xyz or .xyzf file, and breaks it apart 
# into frames. 

import sys
import math

# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
FRAMES = int(sys.argv[1])

#########

# What is the input .xyz file?
XYZFILE = open(sys.argv[2],"r")

#########

ZEROES = len(str(FRAMES))+1


for f in xrange(FRAMES):

	# Generate the output filename

	TAG = ""
	for i in xrange(ZEROES):
		if f == 0:
			if f+1 < pow(10.0,i):
				for j in xrange(ZEROES-i):
					TAG += "0"
				TAG += `f`
				break
			
		elif f < pow(10.0,i):
			for j in xrange(ZEROES-i):
				TAG += "0"
			TAG += `f`
			break

	OUTFILE  = sys.argv[2]
	TESTER   = OUTFILE [0:-4]
	TESTER   = TESTER  [-1]

	if TESTER == ".":
		OUTFILE  = OUTFILE[0:-5] + "_#" + TAG + ".xyzf" 
	else:
		OUTFILE  = OUTFILE[0:-4] + "_#" + TAG + ".xyz" 

	OFSTREAM = open(OUTFILE,"w")

	# Read the first line to get the number of atoms in the frame,
	# print back out to the xyzf file
	
	ATOMS = XYZFILE.readline()
	
	OFSTREAM.write(ATOMS)
	
	ATOMS = ATOMS.split()
	ATOMS = int(ATOMS[0])
	
	# Read/print the comment line
	
	OFSTREAM.write(XYZFILE.readline())
	
	# Now, read/print each atom line
	
	for j in xrange(ATOMS):

		OFSTREAM.write(XYZFILE.readline())
		
	OFSTREAM.close()
	


