import numpy.polynomial.chebyshev as cheby
import numpy.linalg as linalg
import numpy as np
import math
import sys

# Useage is: 
#
# python cheby_fitting.py  <VAR1 > <VAR2 > <VAR3 > <VAR4 > <VAR5 > <VAR6 > <VAR7 > <VAR8 > <VAR9 > <VAR10> <VAR11> <VAR12>
#
#
#
# <VAR1 >: Desired fitting order/start point
# <VAR2 >: Original r_min
# <VAR3 >: Original r_max
# <VAR4 >: Original Morse lambda
# <VAR5 >: Threshold to search for slope change (i.e. from 0 to <VAR5>)
# <VAR6 >: Printing resolution of input PES file
# <VAR7 >: Input PES file
# <VAR8 >: Case (1 == re-fit only between original r_min and r_max; 2 == re-fit between some user-defined r_min and r_max)
# <VAR9 >: Only used if Case 2 is requested: Original simulation temperature
# <VAR10>: Only used if Case 2 is requested: A multiple of kBT to fit to. If (-), will ignore, and extrapolate out to zero, otherwise, will extrapolate as far as needed to have repulsive potential equal to this porduct of kBT
# <VAR11>: Enforce the user-specified order (true) or find the optimal order during fittig process (false)?
# <VAR12>: Range to transform distances to. 1 == -1 to 1; 2 == -1 to 0; 3 == 0 to +1. In general, choosing 1 will result in a larger order, when <VAR11> is true.
#
#
# Note: Suggested RMIN to achieve target KT takes fcut function into consideration.

##################################################################################################
##################################################################################################
###############################		BEGIN SCRIPT	##########################################
##################################################################################################
##################################################################################################



# Define text printing styles/tags

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
# Set up the output log file

LOGSTREAM = open("refit_cheby.log",'w')    
    
    
# Define the desired range of the Cheby fit

RANGE = 1 # -1 to 1
#RANGE = 2 # -1 to 0
#RANGE = 3 #  0 to 1     

# Definen remaining variables

#CASE = 1		# 1. only fit between orig rmin/rmax... 2. Fit from 0 to rmax
ENFORCE_ORDER = False	# If false, a best order will be recommended. Otherwise, the user specified order will be enforced
KT_TARG = 1


ORDER = int  (sys.argv[1]) # 10
RMIN  = float(sys.argv[2]) # 1.0
RMAX  = float(sys.argv[3]) # 4.5
LAMBDA= float(sys.argv[4]) # 1.25
THRESH= float(sys.argv[5]) # 0.25
SPACNG= float(sys.argv[6]) # 0.02 # How far apart are scan data points?
SCANFL=       sys.argv[7]  # HH-NF-25-orig.dat

OUTFL = "refit_scan.dat"

CASE = int (sys.argv[8])

if CASE == 2:
	TEMPER= float(sys.argv[9]) # What temperature was the simulation run at? (K) 
	KT_TARG=float(sys.argv[10]) # How many KT do we want out potential to go up to? (sets how rmin gets modified for a type 2 calculation) .. Either 1000KT or 100KT seem safe-ish

ENFORCE_ORDER = sys.argv[len(sys.argv)-2].lower()
RANGE         = int(sys.argv[len(sys.argv)-1])

print ""
print "Initial order:    " + `ORDER `
print "Initial r_min:    " + `RMIN  `
print "r_max:            " + `RMAX  `
print "Morse lambda:     " + `LAMBDA`
print "Refit threshold:  " + `THRESH`
print "PES scan spacing: " + `SPACNG`
print "Input scan file:  " +  SCANFL
print "Output scan file: " +  OUTFL 

LOGSTREAM.write(""+"\n")
LOGSTREAM.write("Initial order:    " + `ORDER `+"\n")
LOGSTREAM.write("Initial r_min:    " + `RMIN  `+"\n")
LOGSTREAM.write("r_max: 	   " + `RMAX  `+"\n")
LOGSTREAM.write("Morse lambda:     " + `LAMBDA`+"\n")
LOGSTREAM.write("Refit threshold:  " + `THRESH`+"\n")
LOGSTREAM.write("PES scan spacing: " + `SPACNG`+"\n")
LOGSTREAM.write("Input scan file:  " +  SCANFL +"\n")
LOGSTREAM.write("Output scan file: " +  OUTFL  +"\n")

if CASE == 2:
	print "Simulation temp:  " + `TEMPER` + " K"
	LOGSTREAM.write("Simulation temp:  " + `TEMPER` + " K"+"\n")

	if KT_TARG > 0:
		print "Case:             " + `CASE` + " ...modifying r_min to achieve extrapolated energies of at least " + `KT_TARG` + " k_B.T"
		LOGSTREAM.write("Case:             " + `CASE` + " ...modifying r_min to achieve extrapolated energies of at least " + `KT_TARG` + " k_B.T"+"\n")
	else:
		print "Case:             " + `CASE` + " ...resetting r_min to zero."
		LOGSTREAM.write("Case:             " + `CASE` + " ...resetting r_min to zero."+"\n")
print ""
LOGSTREAM.write(""+"\n")

x 	= []
y 	= []
xform_x = []


# Read in the PES scan

IFSTREAM = open(SCANFL,'r').readlines()

for i in xrange(len(IFSTREAM)):
    if len(IFSTREAM[i].split()) >0:
	    x.append(float(IFSTREAM[i].split()[0]))
	    y.append(float(IFSTREAM[i].split()[1]))


# Compute the slopes for all points within RMIN+THRESH

SLOPES = []
XVALUS = []

for i in xrange(len(x)):
	if x[i] < RMIN+THRESH:
		if i == 0:
			x_o = x[i]
			y_o = y[i]
		elif i >= 1:
			SLOPES.append((y[i] - y_o)/(x[i] - x_o))
			XVALUS.append(i)	# Use the large x rather than midpoint to encourage values outside of inflection point
			
			if SLOPES[i-1] > 0:
				print "ERROR: Found a (+)slope within 0 - RMIN+THRESH > 0!"
				print "       All slopes in this domain must be (-)       "
				print "       Choose a different THRESH and try again.    "
				print "       Base THRES on PES scan shape near rmin      "
				
				LOGSTREAM.write("ERROR: Found a (+)slope within 0 - RMIN+THRESH > 0!"+"\n")
				LOGSTREAM.write("	All slopes in this domain must be (-)	    "+"\n")
				LOGSTREAM.write("	Choose a different THRESH and try again.    "+"\n")
				LOGSTREAM.write("	Base THRES on PES scan shape near rmin      "+"\n")
				
				exit()
			
			x_o = x[i]
			y_o = y[i]
			
if len(SLOPES) < 2:
	print "ERROR: Need at least 3 data points"
	LOGSTREAM.write("ERROR: Need at least 3 data points"+"\n")
	exit(0)		

# Search for the minimum of all slopes. This is where the turn-over begins...
# Set extrapolation fitting points based off this location 

MIN_X = 0
MIN_S = 0	    
	 
FOUND_MIN = False 	 
	 
for i in xrange(len(SLOPES)):
	if i == 0:
		MIN_X = XVALUS[i]
		MIN_S = SLOPES[i]
	else:
		if SLOPES[i] < MIN_S:
			MIN_X = XVALUS[i]
			MIN_S = SLOPES[i]

if x[MIN_X] == x[1] and MIN_S < SLOPES[1]:
	print bcolors.WARNING +  bcolors.BOLD + "All slopes within provided threshold exhibit increasing steepness." + bcolors.ENDC
	print bcolors.WARNING +  bcolors.BOLD + "No refitting of existing data points required" + bcolors.ENDC
	
	LOGSTREAM.write("All slopes within provided threshold exhibit increasing steepness."+"\n")
	LOGSTREAM.write("No refitting of existing data points required"+"\n")
	
	MIN_X = 0
else:			
	print bcolors.WARNING +  bcolors.BOLD + "Found min slope at x = " + `x[MIN_X]` + bcolors.ENDC
	LOGSTREAM.write("Found min slope at x = " + `x[MIN_X]`+"\n")
print ""
LOGSTREAM.write(""+"\n")

START_EXTRAP = MIN_X

if MIN_X == 0:	# Do it with the first 3 data points		
	START_EXTRAP = 0
else:	
	START_EXTRAP += 1    

		
# If it does, we will replace these data points with START_EXTRAP and do an exponential
# fit of the data (i.e. y = A*exp(-bx) ==> ln(y) = -bx + ln(y)

x_extr = x[START_EXTRAP:START_EXTRAP+3]
y_extr = y[START_EXTRAP:START_EXTRAP+3]

# Shift the y values to ensure none have a negative value (problematic for the log)

NEG = y_extr[0]

for i in xrange(len(y_extr)-1):
	if y_extr[i+1] < NEG:
		NEG = y_extr[i+1]
		
if NEG <=0:
	NEG = abs(NEG) + 3	
	for i in xrange(len(y_extr)):
		y_extr[i] += NEG
else:
	NEG = 0

print  bcolors.WARNING +  bcolors.BOLD + "Using the following points for exponential extrapolation: " + bcolors.ENDC
print  bcolors.WARNING +  bcolors.BOLD + "x_extr: " + `x_extr` + bcolors.ENDC
print  bcolors.WARNING +  bcolors.BOLD + "y_extr: " + `y_extr` + bcolors.ENDC
print ""

LOGSTREAM.write("Using the following points for exponential extrapolation: "+"\n")
LOGSTREAM.write("x_extr: " + `x_extr` +"\n")
LOGSTREAM.write("y_extr: " + `y_extr` +"\n")
LOGSTREAM.write(""+"\n")


for i in xrange(len(y_extr)):
	y_extr[i] = math.log(y_extr[i])

x_extr = np.array(x_extr)
y_extr = np.array(y_extr)

matr = np.vstack([x_extr, np.ones(len(x_extr))]).T

b,A = linalg.lstsq(matr,y_extr)[0]
b = -1.0/b
A = math.exp(A)

print "Extrapolation based on A.exp(-x/b)... found parameters:"
print "A: " + `A`
print "b: " + `b`
print ""

LOGSTREAM.write( "Extrapolation based on A.exp(-x/b)... found parameters:"+"\n")
LOGSTREAM.write( "A: " + `A`+"\n")
LOGSTREAM.write( "b: " + `b`+"\n")
LOGSTREAM.write( ""+"\n")

BOLTZK  = 1.9872036/1000.0 # In Kcal units
#KT_TARG = 1000 # Either 1000KT or 100KT recommended -- Set at input
if CASE == 2:
	KT_VALU = BOLTZK*TEMPER





# CASE 1: ONLY REFIT WITHIN ORIGINAL RMIN


if CASE == 1:

	print bcolors.WARNING +  bcolors.BOLD + "Refitting according to CASE 1 (only between r_min and r_max)" + bcolors.ENDC
	print ""
	
	LOGSTREAM.write("Refitting according to CASE 1 (only between r_min and r_max)"+"\n")
	LOGSTREAM.write(""+"\n")

	# Now replace the values with new extrapolated values 

	for i in xrange(START_EXTRAP):
		y[i] = A*math.exp(-1.0*x[i]/b) - NEG
	
	# Now re-fit the cheby, and print it out 

	for i in xrange(len(x)): # Morse transformation of pair distance
	
		# Case one: Fit on range -1:1 (When used alone)
		xform_x.append ((math.exp(-1.0*x[i]/LAMBDA) - 0.5 * (math.exp(-1.0*RMIN/LAMBDA) + math.exp(-1.0*RMAX/LAMBDA))) / (0.5 * (math.exp(-1.0*RMIN/LAMBDA) - math.exp(-1.0*RMAX/LAMBDA))))	
		
		if RANGE == 2:
		
			# Case 2: Fit on range -1 to 0 (When used in combo with case 1)
			xform_x[i] = xform_x[i]/2.0 - 0.5
			
		elif RANGE == 3:	
		
			# Case 3: Fit on range 0 to 1 (When used in combo with case 1)
			xform_x[i] = xform_x[i]/2.0 + 0.5		
else:	

	RMIN_NEW = 0

	if KT_TARG > 0:
	
		#"Goal seek" to estimate rmin needed to achieve target KT multiple
	
		test = 0 
		while True:
			
			CHECK = pow((1-test/RMAX),3.0)*A*math.exp(-1.0*test/b) - (KT_TARG*KT_VALU)

			if   CHECK < 0 and test == 0:
				RMIN_NEW = test
				RMIN_NEW = float(int(RMIN_NEW*100.0)/100.0)
				print  bcolors.WARNING +  bcolors.BOLD +"Warning: requested KT product requires r_min < 0." + bcolors.ENDC
				print  bcolors.WARNING +  bcolors.BOLD +"      ...setting RMIN to " + `RMIN_NEW ` + bcolors.ENDC
				
				LOGSTREAM.write("Warning: requested KT product requires r_min < 0."+"\n")
				LOGSTREAM.write("      ...setting RMIN to " + `RMIN_NEW `+"\n")
				
			elif CHECK < 0:
				RMIN_NEW = test
				RMIN_NEW = float(int(RMIN_NEW*100.0)/100.0)
				print "Setting RMIN to " + `RMIN_NEW `+ ", will get up to " + `pow((1-test/RMAX),3.0)*A*math.exp(-1.0*test/b)/KT_VALU` + " k_B.T " 
				LOGSTREAM.write("Setting RMIN to " + `RMIN_NEW `+ ", will get up to " + `pow((1-test/RMAX),3.0)*A*math.exp(-1.0*test/b)/KT_VALU` + " k_B.T " +"\n")
				break
			else:
				test += SPACNG
			
	
		#RMIN_NEW = -b*math.log(KT_TARG*KT_VALU/A) # = 0.0
		
	else:
		RMIN_NEW = 0.0

	print bcolors.WARNING +  bcolors.BOLD + "Refitting according to CASE 2 (between " + `RMIN_NEW` + " and r_max)" + bcolors.ENDC
	print ""
	
	LOGSTREAM.write("Refitting according to CASE 2 (between " + `RMIN_NEW` + " and r_max)"+"\n")
	LOGSTREAM.write(""+"\n")
	
	

	# Now replace the values with new extrapolated values 

	for i in xrange(START_EXTRAP):
		y[i] = A*math.exp(-1.0*x[i]/b)
	
	for i in reversed(xrange(int((RMIN-RMIN_NEW)/SPACNG))):
		if SPACNG*i >= RMIN:
			break
		x.insert(0,(SPACNG*i+RMIN_NEW))
		y.insert(0,A*math.exp(-1.0*(SPACNG*i+RMIN_NEW)/b))


	# Do the morse transform on distances

	for i in xrange(len(x)): # Morse transformation of pair distance

		# Case one: Fit on range -1:1 (When used alone)
		xform_x.append ((math.exp(-1.0*x[i]/LAMBDA) - 0.5 * (math.exp(-1.0*RMIN_NEW/LAMBDA) + math.exp(-1.0*RMAX/LAMBDA))) / (0.5 * (math.exp(-1.0*RMIN_NEW/LAMBDA) - math.exp(-1.0*RMAX/LAMBDA))))
		
		if RANGE == 2:
		
			# Case 2: Fit on range -1 to 0 (When used in combo with case 1)
			xform_x[i] = xform_x[i]/2.0 - 0.5
			
		elif RANGE == 3:	
		
			# Case 3: Fit on range 0 to 1 (When used in combo with case 1)
			xform_x[i] = xform_x[i]/2.0 + 0.5		
		
		
	RMIN = RMIN_NEW	



# Do the fitting. We will run with a threshold in R_SQ_WORST of 20.

R_SQ_WORST = 25.0
ORDER     -= 1 # Only because of while loop

while R_SQ_WORST > 20.0:

	ORDER += 1
	
	full_params = cheby.chebfit(xform_x,y,range(1,ORDER+1),len(x)*2e-16, True)
	params = full_params[0]

	fitted = cheby.chebval(xform_x,params)
	
	if len(full_params[1][0]) > 0:
	
		R_SQ_WORST = full_params[1][0][0] # sum of squared residuals of the least squares fit rank
	else:
		R_SQ_WORST = -1
	
	
	print "Order/R^2: " + `ORDER` + " " + `R_SQ_WORST`
	LOGSTREAM.write("Order/R^2: " + `ORDER` + " " + `R_SQ_WORST`+"\n")
	
	if ENFORCE_ORDER == "true":
		break

# Print out the results


print ""		
print bcolors.WARNING +  bcolors.BOLD + "New Cheby order and r_min are: " + `ORDER` + " " + `RMIN`  + bcolors.ENDC
print ""

LOGSTREAM.write(""+"\n")
LOGSTREAM.write("New Cheby order and r_min are: " + `ORDER` + " " + `RMIN`+"\n")
LOGSTREAM.write(""+"\n")

if   RANGE == 1:
	print bcolors.WARNING +  bcolors.BOLD + "Cheby values fit over range -1 to 1" + bcolors.ENDC
	LOGSTREAM.write("Cheby values fit over range -1 to 1"+"\n")
elif RANGE == 2:
	print bcolors.WARNING +  bcolors.BOLD + "Cheby values fit over range -1 to 0" + bcolors.ENDC
	LOGSTREAM.write("Cheby values fit over range -1 to 0"+"\n")
else:
	print bcolors.WARNING +  bcolors.BOLD + "Cheby values fit over range 0 to 1" + bcolors.ENDC
	LOGSTREAM.write("Cheby values fit over range 0 to 1"+"\n")
print ""		
print bcolors.WARNING +  bcolors.BOLD + "Sum of squared residuals of the least squares fit rank: " + `R_SQ_WORST`  + bcolors.ENDC
print ""


LOGSTREAM.write(""+"\n")
LOGSTREAM.write("Sum of squared residuals of the least squares fit rank: " + `R_SQ_WORST`+"\n")
LOGSTREAM.write(""+"\n")

if R_SQ_WORST > 20.0:
	if ENFORCE_ORDER:
		print "WARNING: User define order has been enforced...   "
		print "         Good fit not guaranteed when R_SQ_WORST  "
		print "         is larger than ~20.0. Consider increasing"
		print "         order or allowing order to auto-adjust   "
		print ""
		
		LOGSTREAM.write("WARNING: User define order has been enforced...   "+"\n")
		LOGSTREAM.write("	  Good fit not guaranteed when R_SQ_WORST  "+"\n")
		LOGSTREAM.write("	  is larger than ~20.0. Consider increasing"+"\n")
		LOGSTREAM.write("	  order or allowing order to auto-adjust   "+"\n")
		LOGSTREAM.write("")

		

print bcolors.OKBLUE + "Found the following parameters: "  + bcolors.ENDC
print ""

LOGSTREAM.write("Found the following parameters: "+"\n")
LOGSTREAM.write(""+"\n")


for i in xrange(len(params)):
	if i > 0:
		print bcolors.OKBLUE +  `i-1` +  "	" + `params[i]`   + bcolors.ENDC
		LOGSTREAM.write(`i-1` +  "	" + `params[i]`+"\n")
		
OUTSTREAM  = open(OUTFL,'w')

for i in xrange(len(x)):
		#OUTSTREAM.write( `x[i]` + " " + `xform_x[i]` + " " + `y[i]` + " " + `cheby.chebval(xform_x[i],params)` + '\n' )	
		OUTSTREAM.write( `x[i]` + " " + `xform_x[i]` + " " + `y[i]` + " " + `cheby.chebval(xform_x[i],params)` + " " + `cheby.chebval(xform_x[i],params)*pow((1-x[i]/RMAX),3.0)` + '\n' )	
OUTSTREAM.close()


print ""
LOGSTREAM.write(""+"\n")
	
TEST_MAX = (math.exp(-1.0*(RMAX)/LAMBDA) - 0.5 * (math.exp(-1.0*RMIN/LAMBDA) + math.exp(-1.0*RMAX/LAMBDA))) / (0.5 * (math.exp(-1.0*RMIN/LAMBDA) - math.exp(-1.0*RMAX/LAMBDA)))
TEST_MIN = (math.exp(-1.0*(RMIN)/LAMBDA) - 0.5 * (math.exp(-1.0*RMIN/LAMBDA) + math.exp(-1.0*RMAX/LAMBDA))) / (0.5 * (math.exp(-1.0*RMIN/LAMBDA) - math.exp(-1.0*RMAX/LAMBDA)))

if RANGE == 2:

	# Case 2: Fit on range -1 to 0 (When used in combo with case 1)
	
	TEST_MAX = TEST_MAX/2.0 - 0.5
	TEST_MIN = TEST_MIN/2.0 - 0.5
	
elif RANGE == 3:	

	# Case 3: Fit on range 0 to 1 (When used in combo with case 1)
	
		TEST_MAX = TEST_MAX/2.0 + 0.5
		TEST_MIN = TEST_MIN/2.0 + 0.5	

#print bcolors.OKBLUE + "Note: Potential has value " + `cheby.chebval(TEST_MAX,params)` + " at r_max (" + `RMAX` + ")"  + bcolors.ENDC
#print bcolors.OKBLUE + "Note: Potential has value " + `cheby.chebval(TEST_MIN,params)` + " at r_min (" + `RMIN` + ")"  + bcolors.ENDC

print ""
print "Run complete."
print ""

LOGSTREAM.write(""+"\n")
LOGSTREAM.write("Run complete."+"\n")
LOGSTREAM.write(""+"\n")
