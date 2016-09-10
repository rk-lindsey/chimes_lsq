import sys

# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
FRAMES = int(sys.argv[1])

# What is the input .xyz file?
XYZFILE = open(sys.argv[2],"r")

# What is the input force file?
FRCILE = open(sys.argv[3],"r")

# What is the input boxlen file?
BOXFILE = open(sys.argv[4],"r")

# What is the outputfile
OUTFILE  = sys.argv[2]
OUTFILE  = OUTFILE[0:-4] + ".xyzf" # replace ".xyz" with ".xyzf"
OFSTREAM = open(OUTFILE,"w")

#########

for i in xrange(FRAMES):
	
	# Read the first line to get the number of atoms in the frame,
	# print back out to the xyzf file
	
	ATOMS = XYZFILE.readline()
	
	OFSTREAM.write(ATOMS)
	
	ATOMS = ATOMS.split()
	ATOMS = int(ATOMS[0])
	
	# Read/print the box lengths (Ignore the frame line in the input xyz file
	XYZFILE.readline()
	OFSTREAM.write(BOXFILE.readline())
	
	# Now, read each atom line, and append the x,y, and z forces before printing
	
	for j in xrange(ATOMS):
	
		# The input forces are (or at least should be) in kcal/mol... lsq c++ expects them in hartree/bohr.
		# Convert units to hartree/bohr before printing
		
		X = float(FRCILE.readline().rstrip('\n'))/627.50960803/1.889725989
		Y = float(FRCILE.readline().rstrip('\n'))/627.50960803/1.889725989
		Z = float(FRCILE.readline().rstrip('\n'))/627.50960803/1.889725989
	
		OFSTREAM.write(XYZFILE.readline().rstrip('\n') + " " + `X `+ " " +  `Y` + " " + `Z` + "\n")
	
