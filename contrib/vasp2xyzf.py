#!/opt/local/bin/python
import sys

# Usage is one of: 
#	python <this script> <OUTCAR file>
#	python <this script> <OUTCAR file> <frame skip> <optional args>
# Optional args are:
#	STRESS or ALLSTR which print diagonal or all stress tensor components, respectively
#	ENERGY which prints overall frame energy
# The frame skip arg tell the script to only output every frame skip'th frame

# Read the input file

SKIP   = 1
OUTCAR = sys.argv[1]
INCL_STRESS = False
INCL_ALLSTR = False
INCL_ENERGY = False


# If requested, read how frequently to print frames.. print every SKIP'th frame

if len(sys.argv) >= 3:
	SKIP = int(sys.argv[2])
	
if len(sys.argv) >= 4:
	if sys.argv[3] == "STRESS":
		INCL_STRESS = True
	if sys.argv[3] == "ALLSTR":
		INCL_ALLSTR = True
	if sys.argv[3] == "ENERGY":
		INCL_ENERGY = True
		
if len(sys.argv) >= 5:
	if sys.argv[4] == "STRESS":
		INCL_STRESS = True
	if sys.argv[4] == "ALLSTR":
                INCL_ALLSTR = True
	if sys.argv[4] == "ENERGY":
		INCL_ENERGY = True		


# Open up file, read some "header" type information

ATOM_TYPES = []
BOXLENS    = []

TMP0 = None
TMP1 = None
TMP2 = None
TIMESTEP = None
TEMPERAT = None

EXIT_CONDITION = False

IFSTREAM = open(OUTCAR, 'r')

A = None
B = None
C = None

POT_ENER = 0.0

DETECT_LENGTHS = -1

while not EXIT_CONDITION:
	
	LINE = IFSTREAM.readline()
	
	if not LINE:
		
		print "	!ERROR: End of file reached: ",OUTCAR
		exit()
	
	if "free  energy   TOTEN" in LINE:
	
		POT_ENER = LINE.split()[4]
		if "***" in POT_ENER:
			POT_ENER = "10000000.0"
		else:
			POT_ENER = float(POT_ENER)*23.06 # Convert from eV (VASP) to kcal/mol (house_md)
		
		EXIT_CONDITION = True	

	# Get the box lengths
	
	if "direct lattice vectors" in LINE:
		DETECT_LENGTHS = 0
		continue
	
	if DETECT_LENGTHS == 0:
		A = LINE.split()[0:3]
		DETECT_LENGTHS += 1
		continue
		
	if DETECT_LENGTHS == 1:
		B = LINE.split()[0:3]
		DETECT_LENGTHS += 1
		continue
		
	if DETECT_LENGTHS == 2:
		C = LINE.split()[0:3]
		DETECT_LENGTHS = -1
		continue	
		
	# Get the number of atoms of each type
	
	if "ions per type" in LINE:
		TMP2 = LINE.split()[4:len(LINE.split())+1]
	
	LINE = LINE.split()
	
	# Get a list of atom types
	
	if len(LINE) > 0 and LINE[0] == "POSCAR:":
		TMP1 = len(LINE)+1
		TMP1 = LINE[1:TMP1+1]
		
	# Get the temperature and timestep
		
	if len(LINE) > 0 and LINE[0] == "POTIM":
		TIMESTEP = LINE[2].rstrip(';')	
		
	if len(LINE) > 0 and LINE[0] == "TEBEG":
		TEMPERAT = LINE[2].rstrip(';')	
	
BOXLENS.append(' '.join(A))
BOXLENS.append(' '.join(B))
BOXLENS.append(' '.join(C))
	
		
IFSTREAM.close()
		
# Generate a list of atoms

TOTAL_ATOMS = 0

for i in xrange(len(TMP2)):
	for j in xrange(int(TMP2[i])):
		ATOM_TYPES.append(TMP1[i])
		TOTAL_ATOMS += 1
TOTAL_FRAMES = 0
PRINTED_FRAMES = 0

print ""
print "================================================================"
print ""
print "Atom types                       : " + ' '.join(TMP1)
print "Number of atoms of each type     : " + ' '.join(TMP2)
print "Total atoms                      : " + `TOTAL_ATOMS` 
print "Simulation temperature           : " + TEMPERAT + " (K)"
print "Simulation timestep              : " + TIMESTEP + " (fs)"
print "Frame skip set to                : " + `SKIP`
print "Time between generated xyz frames: " + `round((SKIP)*float(TIMESTEP),5)` + " (fs)"

if INCL_STRESS:
	print "Include xx, yy, and zz stresses  : true"
else:
	print "Include xx, yy, and zz stresses  : false"
if INCL_ALLSTR:
	print "Include all stresses             : true"
else:
	print "Include all stresses             : false"
if INCL_ENERGY:
	print "Include energies                 : true"
else:	
	print "Include energies                 : false"


# Start printing the .xyz file

OFSTREAM = open(OUTCAR + ".xyzf",'w')

IFSTREAM = open(OUTCAR, 'r')

START_READ = 0

START_TENS = 0
TENS_FRAMES = 0

READ_ATOMS = 0

TENSOR_XX = 0
TENSOR_YY = 0
TENSOR_ZZ = 0
TENSOR_XY = 0
TENSOR_YZ = 0
TENSOR_ZX = 0


for LINE in IFSTREAM:

	if not LINE:
		
		print "ERROR: End of file reached"
		exit()

	if "in kB" in LINE:
		if TOTAL_FRAMES%SKIP == 0: # Then we should print this frame -- Convert kbar to GPa
			LINE = LINE.split()
			TENSOR_XX = float(LINE[2])/10.0
			TENSOR_YY = float(LINE[3])/10.0
			TENSOR_ZZ = float(LINE[4])/10.0
			TENSOR_XY = float(LINE[5])/10.0
			TENSOR_YZ = float(LINE[6])/10.0
			TENSOR_ZX = float(LINE[7])/10.0
			LINE = " ".join(LINE)
			
	if "free  energy   TOTEN" in LINE:
	
		POT_ENER = LINE.split()[4]
		if "***" in POT_ENER:
			POT_ENER = "10000000.0"
		else:
			POT_ENER = float(POT_ENER)*23.06 # Convert from eV (VASP) to kcal/mol (house_md)			
	
	if "POSITION" in LINE:		# Ignore current (title) line and the "----" line that follows
		TOTAL_FRAMES += 1
			
		START_READ = 1
		continue
	
	if START_READ == 1:		# Ignore current (title) line and the "----" line that follows
		if TOTAL_FRAMES%SKIP == 0: # Then we should print this frame
			START_READ = 2
			PRINTED_FRAMES += 1

			# Print the xyz file headers
		
			OFSTREAM.write(`TOTAL_ATOMS` + '\n')
			
			OFSTREAM.write(BOXLENS[0] + " " + BOXLENS[1] + " " + BOXLENS[2].rstrip(')'))
			
			if INCL_STRESS:
				OFSTREAM.write(" " +`TENSOR_XX `+ " " + `TENSOR_YY` + " " + `TENSOR_ZZ`)
			elif INCL_ALLSTR:
				OFSTREAM.write(" " +`TENSOR_XX `+ " " + `TENSOR_YY` + " " + `TENSOR_ZZ` + " " + `TENSOR_XY` + " " + `TENSOR_YZ` + " " + `TENSOR_ZX`) 
			if INCL_ENERGY:
				OFSTREAM.write(" " + `POT_ENER`)
			OFSTREAM.write('\n')

		continue
	
	if START_READ == 2 and READ_ATOMS < TOTAL_ATOMS:
		
		# Coords are already in angstrom, which is good, but the fmatch 
		# code expects forces in Hartree/bohr, while vasp outputs in eV/Angstr
		
		LINE = LINE.split()
		
		if len(LINE)<6:
			LINE = LINE[0:3]
			LINE.append("10000000.0")
			LINE.append("10000000.0")
			LINE.append("10000000.0")
		
		if "***" in LINE[3]:
			LINE[3] = "10000000.0"
		if "***" in LINE[4]:
			LINE[4] = "10000000.0"
		if "***" in LINE[5]:
			LINE[5] = "10000000.0"						
		
		LINE[3] = `float(LINE[3])/1.8897161646320724/27.2114`
		LINE[4] = `float(LINE[4])/1.8897161646320724/27.2114`
		LINE[5] = `float(LINE[5])/1.8897161646320724/27.2114`
		LINE = ' '.join(LINE)
		OFSTREAM.write(ATOM_TYPES[READ_ATOMS] + " " + LINE + '\n')
		READ_ATOMS += 1
		
	elif START_READ == 2:
		START_READ = 0
		READ_ATOMS = 0


print "Counted total frames             : " + `TOTAL_FRAMES`
print "Printed total frames             : " + `PRINTED_FRAMES`
print ""
print "================================================================"
print ""		
		
	
		
