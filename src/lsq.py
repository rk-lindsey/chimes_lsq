#! /usr/bin/python

import sys
import numpy
import scipy.linalg
import math as m

from numpy import *
from numpy.linalg import lstsq
from datetime import *

def is_number(s):
# Test if s is a number.
    try:
        float(s)
        return True
    except ValueError:
        return False



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
    WEIGHTS= numpy.genfromtxt(sys.argv[7],dtype='float')
    DO_WEIGHTING = True

#################################
#################################
#           BEGIN CODE
#################################
#################################


#################################
#   Process input, setup output
#################################

# Use genfromtxt to avoid parsing large files.
A = numpy.genfromtxt(sys.argv[1], dtype='float')
nlines = A.shape[0] 
np = A.shape[1] 

b = numpy.genfromtxt(sys.argv[2], dtype='float') 
nlines2 = b.shape[0] 

hf=open(sys.argv[3],"r").readlines()
eps_fac = float(sys.argv[5]) ;

if ( nlines != nlines2 ):
    print "Error: the number of lines in the input files do not match\n"
    exit(1) 

if np > nlines:
    print "Error: number of variables > number of equations"
    exit(1)


print "! Date ", date.today() ;
print "!"
print "! Number of variables = ", np
print "! Number of equations    = ", nlines

#A=zeros((nlines,np))
#b=zeros((nlines))

if DO_WEIGHTING:
    if ( WEIGHTS.shape[0] != nlines ):
        print "Wrong number of lines in WEIGHTS file"
        exit(1)

#################################
# Apply weighting to A and b
#################################

if DO_WEIGHTING:
    weightedA = dot(diag(WEIGHTS),A)
    weightedb = dot(diag(WEIGHTS),b)
    

#################################
# Do the SVD, process output
#################################

try:
    if DO_WEIGHTING:
        #U,D,VT=numpy.linalg.svd(weightedA)
        # OK to overwrite weightedA.  It is not used to calculate y (predicted forces) below.
        U,D,VT=scipy.linalg.svd(weightedA,overwrite_a=True)
        Dmat=array((transpose(weightedA)))
    else:
        #U,D,VT=numpy.linalg.svd(A)
        # Do not overwrite A.  It is used to calculate y (predicted forces) below.
        U,D,VT=scipy.linalg.svd(A,overwrite_a=False)
        Dmat=array((transpose(A)))  
except LinAlgError:
    sys.stderr.write("SVD algorithm failed")
    exit(1)


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






############

# Zero-out certain params.
#
# Assume: 12 unique 2-b and 111 unqiue 3-b params
# Assume: Only keep 1000 4B params with largest magnitude
#
#strt  = 12+111-1 # 122
#tmp_x = list(x[strt:])
#tmp_i = []
#
#for i in xrange(len(tmp_x)):
#	tmp_x[i] = abs(tmp_x[i])
#	tmp_i.append(i+122)
#	
# Sort parameters in order of decreasing parameter magnitude
#
#tmp_i.sort(key=dict(zip(tmp_i,tmp_x)).get,reverse=True) # sort the indexes
#
#tmp_x.sort(reverse=True)				# Sort the parameters
#
# only keep largest 108 parameters
#
#for i in xrange(len(tmp_x)):
#
#
#	if i+1 <= 1000:
#		#print `tmp_i[i]` + "	" + `tmp_x[i]` + "	" + `x[tmp_i[i]]` + ":	KEEP"
#		tmp_x[i] = 1.0 
#	else:
#		#print `tmp_i[i]` + "	" + `tmp_x[i]` + "	" + `x[tmp_i[i]]` + ":	LOSE"
#		tmp_x[i] = -1.0
#		
#for i in xrange(len(tmp_x)):
#
#	if tmp_x[i] == -1.0:
#		x[tmp_i[i]] = 0.0

		
	

	









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
# Find the ATOM_TRIPS_LINE and ATOM_QUADS_LINE
# Find the TOTAL_TRIPS and TOTAL_QUADS

ATOM_TRIPS_LINE = 0
ATOM_QUADS_LINE = 0
TOTAL_TRIPS = 0
TOTAL_QUADS = 0

for i in range(0, len(hf)):
    print hf[i].rstrip('\n')
    TEMP = hf[i].split()
    if len(TEMP)>3:
        if (TEMP[2] == "TRIPLETS:"):
            TOTAL_TRIPS = TEMP[3]
            ATOM_TRIPS_LINE = i
            
            for j in range(i, len(hf)):
                TEMP = hf[j].split()
                if len(TEMP)>3:
                    if (TEMP[2] == "QUADRUPLETS:"):
                        print hf[j].rstrip('\n')
                        TOTAL_QUADS = TEMP[3]
                        ATOM_QUADS_LINE = j
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

if POTENTIAL == "CHEBYSHEV" :
    TMP = hf[7].split()
    
    if len(TMP) >= 4:
        SNUM_2B = int(TMP[2])
        
        #print "Expecting:"
        #print SNUM_2B
        #print SNUM_3B
        #print SNUM_4B        
elif POTENTIAL == "DFTBPOLY" :
    TMP = hf[7].split()
    SNUM_2B = int(TMP[2])

elif POTENTIAL == "INVRSE_R":
    TMP = hf[7].split()
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

A1 = ""
A2 = ""

P1 = ""
P2 = ""
P3 = ""

# PAIRS, AND CHARGES

# Figure out how many 3B parameters there are

SNUM_3B   = 0
ADD_LINES = 0
COUNTED_COUL_PARAMS = 0 

if TOTAL_TRIPS > 0:
    for t in xrange(0, int(TOTAL_TRIPS)):

        P1 = hf[ATOM_TRIPS_LINE+3+ADD_LINES].split()

        if P1[4] != "EXCLUDED:":
            SNUM_3B +=  int(P1[4])
        
            TOTL = P1[6]
            ADD_LINES += 5
        
            for i in xrange(0, int(TOTL)):
                ADD_LINES += 1
            
# Figure out how many 4B parameters there are

SNUM_4B   = 0
ADD_LINES = 0

if TOTAL_QUADS > 0:
    for t in xrange(0, int(TOTAL_QUADS)):

        P1 = hf[ATOM_QUADS_LINE+3+ADD_LINES].split()
        
        #print "QUAD HEADER", P1
        if P1[7] != "EXCLUDED:":

            SNUM_4B +=  int(P1[7])
        
            TOTL = P1[9]
        
            ADD_LINES += 5

            for i in xrange(0,int(TOTL)):
                ADD_LINES += 1

#print "TOTAL 4B PARAMETERS ", SNUM_4B            

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
        print "q_" + A1 + " x q_" + A2 + " " + `x[TOTAL_PAIRS*SNUM_2B + SNUM_3B + SNUM_4B + i]`
        COUNTED_COUL_PARAMS += 1
    
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

        print "TRIPLETTYPE PARAMS:"
        print "  " + hf[ATOM_TRIPS_LINE+2+ADD_LINES].rstrip() 

        P1 = hf[ATOM_TRIPS_LINE+3+ADD_LINES].split()
        
        #print "HEADER: ", P1

        V0 = P1[1] 
        V1 = P1[2]
        V2 = P1[3]

        if P1[4] == "EXCLUDED:" :
            print "   PAIRS: " + V0 + " " + V1 + " " + V2 + " EXCLUDED:"
            ADD_LINES += 1
        else:
            UNIQ = P1[4]
            TOTL = P1[6].rstrip() 

            print "   PAIRS: " + V0 + " " + V1 + " " + V2 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL 
            print "     index  |  powers  |  equiv index  |  param index  |       parameter       "
            print "   ----------------------------------------------------------------------------"
        
            ADD_LINES += 3
        
            if(t>0):
                ADD_PARAM += 1
        
            for i in xrange(0,int(TOTL)):
                ADD_LINES += 1
                LINE       = hf[ATOM_TRIPS_LINE+2+ADD_LINES].rstrip('\n')
                LINE_SPLIT = LINE.split()
                    
                print LINE + " " + `x[TOTAL_PAIRS*SNUM_2B + TRIP_PAR_IDX+int(LINE_SPLIT[5])]`

            TRIP_PAR_IDX += int(UNIQ)
            COUNTED_TRIP_PARAMS += int(UNIQ)
                #print "COUNTED_TRIP_PARAMS", COUNTED_TRIP_PARAMS
                
        print ""
        
        ADD_LINES += 2

ADD_LINES = 0

COUNTED_QUAD_PARAMS = 0
if TOTAL_QUADS > 0:
    print "QUADRUPLET " + POTENTIAL + " PARAMS \n"
    
    QUAD_PAR_IDX = 0

    for t in xrange(int(TOTAL_QUADS)):
    
        PREV_QUADIDX = 0
                
        #print "ATOM_QUADS_LINE " + str(ATOM_QUADS_LINE+2+ADD_LINES)

        P1 = hf[ATOM_QUADS_LINE+2+ADD_LINES].split()

        #print "P1 " + P1[1] + P1[2] + P1[3] + P1[4] + P1[5] + P1[6]
        
        print "QUADRUPLETYPE PARAMS: " 
        print "  " + hf[ATOM_QUADS_LINE+2+ADD_LINES].rstrip() 

        P1 = hf[ATOM_QUADS_LINE+3+ADD_LINES].split()

        #print P1 

        V0 = P1[1] 
        V1 = P1[2]
        V2 = P1[3]
        V3 = P1[4] 
        V4 = P1[5]
        V5 = P1[6]

        #print "UNIQUE: ", str(UNIQ)
        if P1[7] == "EXCLUDED:" :
            print "   PAIRS: " + V0 + " " + V1 + " " + V2 + " " + V3 + " " + V4 + " " + V5 + " EXCLUDED: " 
            ADD_LINES += 1

        else:
            UNIQ = P1[7]
            TOTL = P1[9].rstrip() 

            print "   PAIRS: " + V0 + " " + V1 + " " + V2 + " " + V3 + " " + V4 + " " + V5 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL 
            print "     index  |  powers  |  equiv index  |  param index  |       parameter       "
            print "   ----------------------------------------------------------------------------"
        
            ADD_LINES += 3
        
            if(t>0):
                ADD_PARAM += 1
        
            for i in xrange(0,int(TOTL)):
                ADD_LINES += 1
                LINE       = hf[ATOM_QUADS_LINE+2+ADD_LINES].rstrip('\n')
                LINE_SPLIT = LINE.split()

                UNIQ_QUAD_IDX = int(LINE_SPLIT[8])
                    #print 'UNIQ_QUAD_IDX', str(UNIQ_QUAD_IDX)

                print LINE + " " + `x[TOTAL_PAIRS*SNUM_2B + COUNTED_TRIP_PARAMS + QUAD_PAR_IDX + UNIQ_QUAD_IDX]`

            QUAD_PAR_IDX += int(UNIQ)
            COUNTED_QUAD_PARAMS += int(UNIQ)

        #print "COUNTED_QUAD_PARAMS", COUNTED_QUAD_PARAMS
            
        print ""
        
        ADD_LINES += 2
        

if FIT_POVER == "true":
    print "P OVER: " + `x[len(x)-1]`
    OVERCOORD_PARAMS = 1
else:
    OVERCOORD_PARAMS = 0
	
mapsfile=open(sys.argv[4],"r").readlines()

print ""

for i in range(0,len(mapsfile)):
	print mapsfile[i].rstrip('\n')
	
print ""
			
print "ENDFILE"		

total_params = TOTAL_PAIRS * SNUM_2B + COUNTED_TRIP_PARAMS + COUNTED_QUAD_PARAMS + COUNTED_COUL_PARAMS + OVERCOORD_PARAMS 
if total_params != len(x) :
    sys.stderr.write( "Error in counting parameters\n") 
    sys.stderr.write("len(x) " + str(len(x)) + "\n") 
    sys.stderr.write("TOTAL_PAIRS " + str(TOTAL_PAIRS) + "\n") 
    sys.stderr.write("SNUM_2B" + str(SNUM_2B) + "\n") 
    sys.stderr.write("COUNTED_TRIP_PARAMS " + str(COUNTED_TRIP_PARAMS) + "\n") 
    sys.stderr.write("COUNTED_QUAD_PARAMS " + str(COUNTED_QUAD_PARAMS) + "\n")
    sys.stderr.write("COUNTED_COUL_PARAMS " + str(COUNTED_COUL_PARAMS) + "\n")
    sys.stderr.write("OVERCOORD_PARAMS " + str(OVERCOORD_PARAMS) + "\n")
    exit(1)
		
if TEST_SUITE_RUN == "do":
	test_suite_params=open("test_suite_params.txt","w")		
	for i in range(0,len(x)):
		phrase = `i` + " " + `x[i]` + '\n'
		test_suite_params.write(phrase)
	test_suite_params.close()


	
#for i in xrange(0,len(D)):
#    if abs(D[i]) > eps:
#        print `i` + " " + `x[i]` + " KEPT" 
#    else:
#        print `i` + " " + `x[i]` + " REMOVED" 














