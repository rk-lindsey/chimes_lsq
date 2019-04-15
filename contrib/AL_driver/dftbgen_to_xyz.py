import sys


def dftbgen_to_xyz(*argv):
	#NOTE: Assumes an orthorhombic box

	# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
	FRAMES = int(argv[0])

	# What is the input file?
	IFSTREAM = open(argv[1],"r")

	SKIP = 1
	if len(argv) == 3:
		SKIP = int(argv[2])

	# What is the outputfile
	OUTFILE  = argv[1]
	OUTFILE  = OUTFILE[0:-4] + ".xyz" # replace ".gen" with ".xyz"
	OFSTREAM = open(OUTFILE,"w")

	BOXFILE  = argv[1]
	BOXFILE  = BOXFILE[0:-4] + ".box" # replace ".gen" with ".xyz"
	BOXSTREAM = open(BOXFILE,"w")


	for i in xrange(FRAMES):
		
		# Read the first line to get the number of atoms in the frame
		
		ATOMS = IFSTREAM.readline()
		ATOMS = ATOMS.split()
		ATOMS = int(ATOMS[0])
		
		# Read the next line to get the atom types
		
		SYMBOLS = IFSTREAM.readline()
		SYMBOLS = SYMBOLS.split()
		
		# Print the header bits of the xyz file
		
		if (i+1)%SKIP == 0:
		
			OFSTREAM.write(`ATOMS` + '\n')
			OFSTREAM.write("Frame " + `i+1` + '\n')
		
		# Now read/print all the atom lines in the present frame
		
		for j in xrange(ATOMS):
		
			LINE = IFSTREAM.readline()
			LINE = LINE.split()
			
			# Replace the atom type index with a chemical symbol
		
			for k in xrange(len(SYMBOLS)):
				if k+1 == int(LINE[1]):
					LINE[1] = SYMBOLS[k]
					break
					
			# Print out the line
			
			if (i+1)%SKIP == 0:
			
				OFSTREAM.write(' '.join(LINE[1:len(LINE)]) + '\n')
			
		# Finally, read the box lengths... assume cubic
		
		LINE = IFSTREAM.readline()	# Cell angles?
		
		LINE = IFSTREAM.readline().split()	
		X = LINE[0]
		
		LINE = IFSTREAM.readline().split()	
		Y = LINE[1]
		
		LINE = IFSTREAM.readline().split()	
		Z = LINE[2]	
		
		if (i+1)%SKIP == 0:
		
			BOXSTREAM.write(X + " " + Y + " " + Z + '\n')
			
	return




