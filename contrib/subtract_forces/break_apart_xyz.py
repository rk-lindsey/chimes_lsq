# Takes as input a number of frames and a .xyz or .xyzf file, and breaks it apart 
# into frames. 

# Optional: break into chunks of (3rd arg)

import sys
import math

print "WARNING: Converting forces from Hartree/bohr to simulation units (kca/mol/Ang)"


# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
FRAMES = int(sys.argv[1])

#########

# What is the input .xyz file?
XYZFILE = open(sys.argv[2],"r")

CHUNK_LEN = 1
if len(sys.argv) == 4:
	CHUNK_LEN = int(sys.argv[3])

#########

ZEROES = len(str(FRAMES))+1


for f in xrange(FRAMES):

	if f%CHUNK_LEN == 0:
	
		if f > 1:
			OFSTREAM.close()
			FRSTREAM.close()

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
		FORCES   = sys.argv[2]
		TESTER   = OUTFILE [0:-4]
		TESTER   = TESTER  [-1]

		if TESTER == ".":
			FORCES   = OUTFILE[0:-5] + "_FORCES_#" + TAG + ".xyzf" 
			OUTFILE  = OUTFILE[0:-5] + "_#"        + TAG + ".xyzf" 
			
		else:
			FORCES   = OUTFILE[0:-4] + "_FORCES_#" + TAG + ".xyz" 
			OUTFILE  = OUTFILE[0:-4] + "_#"        + TAG + ".xyz" 
			

		OFSTREAM = open(OUTFILE,"w")
		FRSTREAM = open(FORCES,"w")
	
	

	# Read the first line to get the number of atoms in the frame,
	# print back out to the xyzf file
	
	ATOMS = XYZFILE.readline()
	
	OFSTREAM.write(ATOMS)
	
	ATOMS = ATOMS.split()
	ATOMS = int(ATOMS[0])
	
	# Read/print the comment line

	OFSTREAM.write( XYZFILE.readline())

	# Now, read/print each atom line
	
	for j in xrange(ATOMS):
	
		LINE = XYZFILE.readline()

		OFSTREAM.write(LINE)
		
		LINE = LINE.split()

		if len(LINE)>4:
			FRSTREAM.write(`float(LINE[4])*(627.50960803*1.889725989)` + '\n')
			FRSTREAM.write(`float(LINE[5])*(627.50960803*1.889725989)` + '\n')
			FRSTREAM.write(`float(LINE[6])*(627.50960803*1.889725989)` + '\n')
			
#			FRSTREAM.write(LINE[4]  + '\n')
#			FRSTREAM.write(LINE[5]  + '\n')
#			FRSTREAM.write(LINE[6]  + '\n')

	
		

	


