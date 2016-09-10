import sys

#NOTE: Assumes an orthorhombic box

# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
FRAMES = int(sys.argv[1])

# What is the input file?
IFSTREAM = open(sys.argv[2],"r")

# What is the outputfile
OUTFILE  = sys.argv[2]
OUTFILE  = OUTFILE[0:-4] + ".xyz" # replace ".gen" with ".xyz"
OFSTREAM = open(OUTFILE,"w")

BOXFILE  = sys.argv[2]
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
		
		OFSTREAM.write(' '.join(LINE[1:len(LINE)]) + '\n')
		
	# Finally, read the box lengths... assume cubic
	
	LINE = IFSTREAM.readline()	# Cell angles?
	
	LINE = IFSTREAM.readline().split()	
	X = LINE[0]
	
	LINE = IFSTREAM.readline().split()	
	Y = LINE[1]
	
	LINE = IFSTREAM.readline().split()	
	Z = LINE[2]	
	
	
	BOXSTREAM.write(X + " " + Y + " " + Z + '\n')
