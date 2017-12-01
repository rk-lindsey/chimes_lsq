# ! /usr/bin/python

import sys
import numpy
import scipy.linalg
import math as m

from numpy import *
from numpy.linalg import lstsq
from datetime import *

TEST_SUITE_RUN = False
DO_WEIGHTING   = False
WEIGHTS        = None

if(len(sys.argv) < 6):
    print "Usage is: ./lsq.py A.txt b.txt params.header ff_groups.map eps_fac"
    print "or "
    print "./lsq.py A.txt b.txt params.header ff_groups.map eps_fac TESTING"
    print "or "
    print "./lsq.py A.txt b.txt params.header ff_groups.map eps_fac WEIGHTS weight_file.txt" 
    sys.exit()
    
if(len (sys.argv) == 7): 
    # Then this run is being used for the test suite... 
    # print out parameters without any fanciness so tolerances can be checked  
    TEST_SUITE_RUN = "do"

if(len (sys.argv) == 8): 
    # Then this run is being used for the test suite... 
    # print out parameters without any fanciness so tolerances can be checked  
    WEIGHTS=open(sys.argv[7],"r").readlines()
    DO_WEIGHTING = True

#################################
#################################
#           BEGIN CODE
#################################
#################################


#################################
#   Process input, setup output
#################################

af=open(sys.argv[1],"r").readlines()
bf=open(sys.argv[2],"r").readlines()
hf=open(sys.argv[3],"r").readlines()
eps_fac = float(sys.argv[5]) ;

nlines=len(af)
nlines2=len(bf)
if ( nlines != nlines2 ):
    print "Error: the number of lines in the input files do not match\n"
    exit(1) 

np=len(af[0].split())
print "! Date ", date.today() ;
print "!"
print "! Number of variables = ", np
print "! Number of equations    = ", nlines

A=zeros((nlines,np))
b=zeros((nlines))

for i in range(0,nlines):
    af[i]=af[i].split()
    
    for n in range(0,np):
        A[i][n]=float(af[i][n])
        
    b[i]=float(bf[i])
    
    if DO_WEIGHTING:
        WEIGHTS[i] = float(WEIGHTS[i])

af = 0
bf = 0

#################################
# Apply weighting to A and b
#################################

if DO_WEIGHTING:
    weightedA = dot(diag(WEIGHTS),A)
    weightedb = dot(diag(WEIGHTS),b)
    

#################################
# Do the SVD, process output
#################################

if DO_WEIGHTING:
    #U,D,VT=numpy.linalg.svd(weightedA)
    U,D,VT=scipy.linalg.svd(weightedA,overwrite_a=True)
    Dmat=array((transpose(weightedA)))
else:
    #U,D,VT=numpy.linalg.svd(A)
    U,D,VT=scipy.linalg.svd(A,overwrite_a=True)
    Dmat=array((transpose(A)))  

dmax = 0.0

for i in range(0,len(Dmat)):
    if ( abs(D[i]) > dmax ) :
        dmax = abs(D[i])
        
    for j in range(0,len(Dmat[i])):
        Dmat[i][j]=0.0

# Cut off singular values based on fraction of maximum value as per 
# numerical recipes (LEF).
eps=eps_fac * dmax
nvars = 0
for i in xrange(0,len(D)):
    if abs(D[i]) > eps:
        Dmat[i][i]=1.0/D[i]
        nvars += 1
        
print "! eps (= eps_fac*dmax) =  ", eps        

x=dot(transpose(VT),Dmat)

if DO_WEIGHTING:
    x=dot(x,dot(transpose(U),weightedb))
else:
    x=dot(x,dot(transpose(U),b))

y=dot(A,x)
Z=0.0

# Put calculated forces in force.txt
yfile = open("force.txt", "w")
for a in range(0,len(b)):
    Z=Z+(y[a]-b[a])**2.0
    yfile.write("%13.6e\n"% y[a]) 


bic = float(nlines) * log(Z/float(nlines)) + float(nvars) * log(float(nlines))

print "! RMS force error = " , sqrt(Z/float(nlines))
print "! max abs variable = ",  max(abs(x))
print "! SVD regularization factor = ", eps_fac
print "! number of SVD fitting vars = ", nvars
print "! Bayesian Information Criterion =  ", bic
print "!"


####################################
# Actually process the header file...
####################################

BREAK_COND = False

# Figure out whether we have triplets and/or quadruplets

for i in range(0, len(hf)):
    sys.stdout.write(hf[i])
    TEMP = hf[i].split()
    if len(TEMP)>3:
        if (TEMP[2] == "TRIPLETS:"):
            for j in range(i, len(hf)):
                TEMP = hf[j].split()
                if len(TEMP)>3:
                    if (TEMP[2] == "QUADRUPLETS:"):
                        sys.stdout.write(hf[j])
                        BREAK_COND = True
                        break
        if (BREAK_COND):
             break

# 1. Figure out what potential type we have

POTENTIAL = hf[7].split()
POTENTIAL = POTENTIAL[1]

print ""

print "PAIR " + POTENTIAL + " PARAMS \n"

# 2. Figure out how many coeffs each atom type will have

CASE_BY_CASE = False

SNUM_2B = 0
SNUM_3B = 0
SNUM_4B = 0

if POTENTIAL == "CHEBYSHEV" or POTENTIAL == "DFTBPOLY":
    TMP = hf[7].split()
    
    if len(TMP) >= 4:
        SNUM_3B = int(TMP[3])
        
        if len(TMP) >= 5:
            SNUM_4B = int(TMP[4])
        
        SNUM_2B = int(TMP[2])
        
        #print "Expecting:"
        #print SNUM_2B
        #print SNUM_3B
        #print SNUM_4B        

elif POTENTIAL == "INVRSE_R":
    TMP = hf[6].split()
    SNUM_2B = int(TMP[2])
else:
    CASE_BY_CASE = True # We'll need to do it per atom pair type.	

# 3. Print out the parameters

FIT_COUL = hf[1].split()
FIT_COUL = FIT_COUL[1]

FIT_POVER = hf[3].split()
FIT_POVER = FIT_POVER[1]

USE_POVER = hf[2].split()
USE_POVER = USE_POVER[1]

ATOM_TYPES_LINE=9

TOTAL_ATOM_TYPES = hf[ATOM_TYPES_LINE].split()

TOTAL_ATOM_TYPES = int(TOTAL_ATOM_TYPES[2])

ATOM_PAIRS_LINE=11+TOTAL_ATOM_TYPES+2

TOTAL_PAIRS =  hf[ATOM_PAIRS_LINE].split()
TOTAL_PAIRS = int(TOTAL_PAIRS[2])

# ... figure out which line "ATOM PAIR TRIPLETS: " is on 

ATOM_TRIPS_LINE=0
if FIT_POVER == "true" or USE_POVER == "true":
    ATOM_TRIPS_LINE=ATOM_PAIRS_LINE+2+TOTAL_PAIRS+2+TOTAL_PAIRS+2
else:
    ATOM_TRIPS_LINE=ATOM_PAIRS_LINE+2+TOTAL_PAIRS+2

TEST = hf[ATOM_TRIPS_LINE].split()
if len(TEST) == 3 and TEST[2] == "CUBIC":
	ATOM_TRIPS_LINE += 2
	
TEST = hf[ATOM_TRIPS_LINE].split()
if TEST[0] == "FCUT":
	ATOM_TRIPS_LINE += 2	
	
TEST = hf[ATOM_TRIPS_LINE].split()
if len(TEST) == 4 and TEST[2] == "S_MINIM:":
    if TEST[3] == "ALL":
        ATOM_TRIPS_LINE += 2
if len(TEST) == 5 and TEST[2] == "S_MINIM:":	
        ATOM_TRIPS_LINE += 2 + int(TEST[4])   
        
TEST = hf[ATOM_TRIPS_LINE].split()
if len(TEST) == 4 and TEST[2] == "S_MAXIM:":
    if TEST[3] == "ALL":
        ATOM_TRIPS_LINE += 2
if len(TEST) == 5 and TEST[2] == "S_MAXIM:":	
        ATOM_TRIPS_LINE += 2 + int(TEST[4])       



# Figure out which line "ATOM PAIR TRIPLETS" is on

TOTAL_TRIPS =  hf[ATOM_TRIPS_LINE].split()
TOTAL_TRIPS = int(TOTAL_TRIPS[3])

# ... figure out which line "ATOM PAIR QUADRUPLETS: " is on 

ATOM_QUADS_LINE = 0
TOTAL_QUADS     = 0

#if TOTAL_TRIPS == 0:
ATOM_QUADS_LINE = 1 + ATOM_TRIPS_LINE
TOTAL_QUADS = hf[ATOM_TRIPS_LINE+1].split()

if len(TOTAL_QUADS)>0:
    TOTAL_QUADS = int(TOTAL_QUADS[3])

A1 = ""
A2 = ""

P1 = ""
P2 = ""
P3 = ""

# PAIRS, AND CHARGES

# Figure out how many 3B parameters there are

SNUM_3B   = 0
ADD_LINES = 0

if TOTAL_TRIPS > 0:
    for t in xrange(0, int(TOTAL_TRIPS)):

        P1 = hf[ATOM_TRIPS_LINE+2+ADD_LINES].split()
        
        #print P1
        
        SNUM_3B +=  int(P1[4])
        
        TOTL = P1[6]
        
        ADD_LINES += 4
        
        for i in xrange(0,int(TOTL)):
            ADD_LINES += 1

for i in range(0,TOTAL_PAIRS):
    
    A1 = hf[ATOM_PAIRS_LINE+2+i+1].split()
    A2 = A1[2]
    A1 = A1[1]
    
    print "PAIRTYPE PARAMS: " + `i` + " " + A1 + " " + A2 + "\n"
    
    if CASE_BY_CASE:  
        
        # Figure out individual SNUM_2B
        
        MIN = hf[ATOM_PAIRS_LINE+2+i+1].split()
        MAX = float(MIN[4])
        DEL = float(MIN[5])
        MIN = float(MIN[3])
        
        SNUM_2B = int((2+m.floor((MAX - MIN)/DEL))*2)
        
    for j in range(0, int(SNUM_2B)):
        print `j` + " " + `x[i*SNUM_2B+j]`
    
    if FIT_COUL == "true":
        print "q_" + A1 + " x q_" + A2 + " " + `x[TOTAL_PAIRS*SNUM_2B + SNUM_3B + i]`
    
    print " "

# TRIPLETS

ADD_LINES = 0
ADD_PARAM = 0

COUNTED_TRIP_PARAMS = 0

if TOTAL_TRIPS > 0:
    print "TRIPLET " + POTENTIAL + " PARAMS \n"
    
    TRIP_PAR_IDX = 0
    
    for t in xrange(0, int(TOTAL_TRIPS)):
    
        PREV_TRIPIDX = 0

        P1 = hf[ATOM_TRIPS_LINE+2+ADD_LINES].split()
        
        UNIQ = P1[4]
        TOTL = P1[6]
        
        P2 = P1[2]
        P3 = P1[3]
        P1 = P1[1]
        
        print "TRIPLETTYPE PARAMS: " +`t` + " " + P1 + " " + P2 + " " + P3 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL + "\n"
        print "     index  |  powers  |  equiv index  |  param index  |       parameter       "
        print "   ----------------------------------------------------------------------------"
        
        ADD_LINES += 2
        
        if(t>0):
            ADD_PARAM += 1
        
        for i in xrange(0,int(TOTL)):
            ADD_LINES += 1
            LINE       = hf[ATOM_TRIPS_LINE+2+ADD_LINES].rstrip('\n')
            LINE_SPLIT = LINE.split()

            print LINE + " " + `x[TOTAL_PAIRS*SNUM_2B + TRIP_PAR_IDX+int(LINE_SPLIT[5])]`

        TRIP_PAR_IDX += int(UNIQ)
        COUNTED_TRIP_PARAMS += TRIP_PAR_IDX
            
        print ""
        
        ADD_LINES += 2

if TOTAL_QUADS > 0:
    print "QUADRUPLET " + POTENTIAL + " PARAMS \n"
    
    QUAD_PAR_IDX = 0
    
    for t in xrange(0, int(TOTAL_QUADS)):
    
        PREV_QUADIDX = 0
                
        P1 = hf[ATOM_QUADS_LINE+2].split()
        
        UNIQ = P1[8]
        TOTL = P1[10]
        
        P2 = P1[2]
        P3 = P1[3]
        P4 = P1[4]
        P5 = P1[5]
        P6 = P1[6]
        P1 = P1[1]
        
        print "QUADRUPLETYPE PARAMS: " +`t` + " " + P1 + " " + P2 + " " + P3 + " " + P4 + " " + P5 + " " + P6 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL + "\n"
        print "     index  |  powers  |  equiv index  |  param index  |       parameter       "
        print "   ----------------------------------------------------------------------------"
        
        ADD_LINES += 2
        
        if(t>0):
            ADD_PARAM += 1
        
        for i in xrange(0,int(TOTL)):
            ADD_LINES += 1
            LINE       = hf[ATOM_QUADS_LINE+2+ADD_LINES].rstrip('\n')
            LINE_SPLIT = LINE.split()

            print LINE + " " + `x[TOTAL_PAIRS*SNUM_2B + COUNTED_TRIP_PARAMS + QUAD_PAR_IDX+int(LINE_SPLIT[8])]`
            #print LINE + " " + `x[TOTAL_PAIRS*SNUM_2B + COUNTED_TRIP_PARAMS]`

        QUAD_PAR_IDX += int(UNIQ)
            
        print ""
        
        ADD_LINES += 2        

if FIT_POVER == "true":
	print "P OVER: " + `x[len(x)-1]`
	
mapsfile=open(sys.argv[4],"r").readlines()

print ""

for i in range(0,len(mapsfile)):
	print mapsfile[i].rstrip('\n')
	
print ""
			
print "ENDFILE"		
		
		
if TEST_SUITE_RUN == "do":
	test_suite_params=open("test_suite_params.txt","w")		
	for i in range(0,len(x)):
		phrase = `i` + " " + `x[i]` + '\n'
		test_suite_params.write(phrase)
	test_suite_params.close()
		
# OLD WAY:
#for i in range(0,len(x)):
#    print i,x[i]




















