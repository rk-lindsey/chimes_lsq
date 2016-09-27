import sys

# USEAGE:  python this_script.py <OUTCAR FILE> <POSCAR MAPPING FILE>
# OUTPUTS: output_vasp.xyzf
# NOTE:    Assumes an orthorhombic box

ATOMS = 0;
FIRST_ATOMS   = True
FIRST_BOXDIMS = True
XDIM = 0.0
YDIM = 0.0
ZDIM = 0.0

XCOORD = []
YCOORD = []
ZCOORD = []

XFORCE = []
YFORCE = []
ZFORCE = []


# Read all needed stuff from OUTCAR file

OFSTREAM = open("output_vasp.xyzf",'w')

OUTCAR_STREAM = open(sys.argv[1])

while 1:
	LINE = OUTCAR_STREAM.readline()

	if FIRST_ATOMS and "POSCAR" in LINE:	# How many atoms are in the system?

		ATOMS = int(LINE.split()[6])
		FIRST_ATOMS = False
		OFSTREAM.write(`ATOMS` + '\n')
		
	if FIRST_BOXDIMS and "LATTYP" in LINE:	# What are the box dimensions?
		ZDIM = float(OUTCAR_STREAM.readline().split()[2])
		XDIM = float(OUTCAR_STREAM.readline().split()[2])*ZDIM
		YDIM = float(OUTCAR_STREAM.readline().split()[2])*ZDIM
		
		OFSTREAM.write(`XDIM` + " " + `YDIM` + " " + `ZDIM` + '\n')
		
		FIRST_BOXDIMS = False
		
			
	if "TOTAL-FORCE" in LINE:
		OUTCAR_STREAM.readline()	# "----"
		
		for i in xrange(ATOMS):
		
			LINE = OUTCAR_STREAM.readline().split()
		
			XCOORD.append(LINE[0])
			YCOORD.append(LINE[1])
			ZCOORD.append(LINE[2])
			
			# Convert from eV/A to Hartree/Bohr
		
			XFORCE.append(float(LINE[3])/27.2114/1.889725989)
			YFORCE.append(float(LINE[4])/27.2114/1.889725989)
			ZFORCE.append(float(LINE[5])/27.2114/1.889725989)
			
		OUTCAR_STREAM.close()
		break
			

# Read all needed stuff from mapper file

MMMD_IDX = []
VASP_IDX = []
ATOM_TYP = []

with open(sys.argv[2]) as MAP_STREAM:		# What is the name of the OUTCAR file?
	for LINE in MAP_STREAM:
		TEMP = LINE.split()
		MMMD_IDX.append(int(TEMP[0]))
		VASP_IDX.append(int(TEMP[1]))
		ATOM_TYP.append(TEMP[2])

MAP_STREAM.close()
		
# Print out the xyz file

for i in xrange(len(MMMD_IDX)):
	OFSTREAM.write(ATOM_TYP[i] + " " + XCOORD[VASP_IDX[i]] + " " + YCOORD[VASP_IDX[i]] + " " + ZCOORD[VASP_IDX[i]] + " " + `XFORCE[VASP_IDX[i]]` + " " + `YFORCE[VASP_IDX[i]]` + " " + `ZFORCE[VASP_IDX[i]]` + "\n")

OFSTREAM.close()
			
			
		
		
		
		
		
		
