# Takes as input a single-freame .xyzf file and a forceout.txt file, replaces forces in .xyzf 
# file with those forces minus those in the forceout file 

import sys

#print sys.argv

# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
FRAMES = int(sys.argv[1])

# What is the input .xyz file?
XYZFILE = open(sys.argv[2],"r")

# What is the input force file?
FRCILE = open(sys.argv[3],"r")

# What is the outputfile
OUTFILE  = sys.argv[2]
OUTFILE  = OUTFILE[0:-5] + "_subtracted.xyzf" 
OFSTREAM = open(OUTFILE,"w")

#print "Will process " + `FRAMES` + " Frames."

#########

for i in xrange(FRAMES):
	
	# Read the first line to get the number of atoms in the frame,
	# print back out to the xyzf file
	
	ATOMS = XYZFILE.readline()
	
	OFSTREAM.write(ATOMS)
	
	ATOMS = ATOMS.split()
	ATOMS = int(ATOMS[0])
	
	# Read/print the box lengths (Ignore the frame line in the input xyz file
	OFSTREAM.write(XYZFILE.readline())
	
	# Now, read/print each atom line
	
	for j in xrange(ATOMS):
	
		# The input forces are (or at least should be) in kcal/mol... lsq c++ expects them in hartree/bohr.
		# Convert units to hartree/bohr before printing
		
		
		ORIG_A  = XYZFILE.readline().split()
		COORD_X = float(ORIG_A[1])
		COORD_Y = float(ORIG_A[2])
		COORD_Z = float(ORIG_A[3])
		ORIG_X  = float(ORIG_A[4])
		ORIG_Y  = float(ORIG_A[5])
		ORIG_Z  = float(ORIG_A[6])
		ORIG_A  = ORIG_A[0]
		
		X = float(FRCILE.readline().rstrip('\n'))
		Y = float(FRCILE.readline().rstrip('\n'))
		Z = float(FRCILE.readline().rstrip('\n'))
	
		OFSTREAM.write(ORIG_A + "\t" + `COORD_X` + "\t" + `COORD_Y` + "\t" + `COORD_Z ` + "\t" + 
		 `ORIG_X - X` + "\t" + `ORIG_Y - Y` + "\t" + `ORIG_Z - Z` + "\n")

	
