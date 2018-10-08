
import sys
import numpy
import scipy.linalg
import math as m
import subprocess
import os

import argparse

from numpy import *
from numpy.linalg import lstsq
from datetime import *
from subprocess import call

parser = argparse.ArgumentParser(description='Least-squares force matching based on output of chimes_lsq')


def is_number(s):
# Test if s is a number.
    try:
        float(s)
        return True
    except ValueError:
        return False


def write_matrix_market(mat, f):
# Print an array in matrix market format.
# Matrix market uses column-oriented storage.
# Only dense 1 and 2d arrays supported.    
    ff = open(f, "w")
    ff.write("%%MatrixMarket matrix array real general\n")
    ff.write("%%Comment\n") 

    if len(mat.shape) == 1:
        ff.write(str(mat.shape[0]) + " " + "1\n")
        for i in range(0, mat.shape[0]):
            ff.write(str(mat[i])+"\n")
    elif len(mat.shape) == 2:
        ff.write(str(mat.shape[0]) + " " + str(mat.shape[1]) + "\n")
        for j in range(0, mat.shape[1]):
            for i in range(0, mat.shape[0]):
                ff.write(str(mat[i,j])+"\n")


def read_matrix_market(f):
# Read an array in matrix market format.
# Matrix market uses column-oriented storage.
# Only dense 1 and 2d arrays supported.    
    ff = open(f, "r")
    dat=ff.readlines()

    sz=dat[1].split()
    szi = [int(sz[0]), int(sz[1])]

    print "ARRAY SIZE = " + str(szi[0]) + " " + str(szi[1])

    if szi[0] == 1:
        ar = numpy.zeros(szi[1], dtype=float)
        print "SHAPE: " + str(ar.shape[0])
        for j in range(0, szi[1]):
            ar[j] = float(dat[j+2])
    else:
        print "UNSUPPORTED ARRAY TYPE"
        sys.exit(1)
        
    return ar

def count_nonzero_vars(x):
    np = 0
    for i in xrange(0, len(x)):
        if ( abs(x[i]) > 1.0e-05 ):
            np = np + 1
    return np

def fit_owlqn(A,b):
## Use the OWLQN fitting algorithm
    nvars = np
    print '! OWLQN algorithm for LASSO used'
    print '! OWLQN alpha = ' + str(alpha_val)
    write_matrix_market(A, 'Amm.txt')
    write_matrix_market(b, 'bmm.txt')
    path=sys.path[0]
    exepath=path[:-3] + "contrib/owlqn/source/owlqn"
    if os.path.exists(exepath):
        command = exepath + " Amm.txt bmm.txt " + str(alpha_val) + " fit.txt -ls -tol 1.0e-05 -m 20"
        print "Running " + command

        os.system(command)
        x = read_matrix_market("fit.txt")
        return x
    else:
        print exepath + " does not exist"
        sys.exit(1)
    
# Arguments supported by the lsq code.
parser.add_argument("--algorithm", default='svd', help='fitting algorithm')
parser.add_argument("--A", default='A.txt',help='A (derivative) matrix') 
parser.add_argument("--b", default='b.txt',help='b (force) file')
parser.add_argument("--header", default='params.header',help='parameter file header')
parser.add_argument("--map", default='ff_groups.map',help='parameter file map')
parser.add_argument("--eps", type=float, default=1.0e-05,help='svd regularization')
parser.add_argument("--weights", default="None",help='weight file')
parser.add_argument("--test-suite", type=bool, default=False,help='output for test suite')
parser.add_argument("--alpha", type=float, default=1.0e-04,help='Lasso or ridge regularization')
parser.add_argument("--folds",type=int, default=4,help="Number of CV folds")

args        = parser.parse_args()

# Fitting algorithm to use.
algorithm   = args.algorithm

# Matrix of force derivatives (a. la. A x = b )
Afile       = args.A

# Matrix of forces
bfile       = args.b

# Force field header
header_file = args.header

# Pair mapping
map_file    = args.map

# SVD regularization
eps_fac     = args.eps

# Weights to use for each force
weights_file = args.weights

# Is this running in the test suite.
test_suite_run = args.test_suite

if ( weights_file == "None" ):
    DO_WEIGHTING = False 
else:
    DO_WEIGHTING = True

if (algorithm.find("lasso") >= 0 or algorithm.find("ridge") >= 0
    or algorithm.find("lassolars") >= 0 ) :
    from sklearn import linear_model

# Regularization value.
alpha_val = args.alpha

alpha_ar = [1.0e-06, 3.2e-06, 1.0e-05, 3.2e-05, 1.0e-04, 3.2e-04, 1.0e-03, 3.2e-03]

# Number of folds for cross-validation.
folds = args.folds

max_iter_num = 100000

test_suite_run = False

    # Then this run is being used for the test suite... 
    # print out parameters without any fanciness so tolerances can be checked  
if ( DO_WEIGHTING ):
    WEIGHTS= numpy.genfromtxt(weights_file,dtype='float')

#################################
#################################
#           BEGIN CODE
#################################
#################################


#################################
#   Process input, setup output
#################################

# Use genfromtxt to avoid parsing large files.
A = numpy.genfromtxt(Afile, dtype='float')
nlines = A.shape[0] 
np = A.shape[1] 

b = numpy.genfromtxt(bfile, dtype='float') 
nlines2 = b.shape[0] 

hf=open(header_file,"r").readlines()

if ( nlines != nlines2 ):
    print "Error: the number of lines in the input files do not match\n"
    exit(1) 

if np > nlines:
    print "Error: number of variables > number of equations"
    exit(1)


print "! Date ", date.today() 
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

weightedA = None
weightedb = None

if DO_WEIGHTING:

    # This way requires too much memory for long A-mat's
    # to avoid a memory error, we will do it the slow way instead:
    #
    #weightedA = dot(diag(WEIGHTS),A)
    #weightedb = dot(diag(WEIGHTS),b)
    
    weightedA = numpy.zeros((A.shape[0],A.shape[1]),dtype=float)
    weightedb = numpy.zeros((A.shape[0],),dtype=float)
    
    for i in xrange(A.shape[0]):     # Loop over rows (atom force components)
        for j in xrange(A.shape[1]): # Loop over cols (variables in fit)
	    weightedA[i][j] = A[i][j]*WEIGHTS[i]
	    weightedb[i] = b[i]*WEIGHTS[i]
	    
    

#################################
# Do the SVD, process output
#################################

if algorithm == 'svd':
    print '! svd algorithm used'
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
    print "! SVD regularization factor = ", eps_fac

    x=dot(transpose(VT),Dmat)

    if DO_WEIGHTING:
        x=dot(x,dot(transpose(U),weightedb))
    else:
        x=dot(x,dot(transpose(U),b))

elif algorithm == 'ridge':
    print '! ridge regression used'
    reg = linear_model.Ridge(alpha=alpha_val,fit_intercept=False)

    # Fit the data.
    reg.fit(A,b)
    
    x = reg.coef_
    nvars = np
    print '! Ridge alpha = ' + str(alpha_val)

elif algorithm == 'ridgecv':
    reg = linear_model.RidgeCV(alphas=alpha_ar,fit_intercept=False,cv=folds)
    reg.fit(A,b)
    print '! ridge CV regression used'
    print '! ridge CV alpha = ' + str(reg.alpha_)
    x = reg.coef_
    nvars = np

elif algorithm == 'lasso':
    print '! Lasso regression used'
    print '! Lasso alpha = ' + str(alpha_val)
    reg = linear_model.Lasso(alpha=alpha_val,fit_intercept=False,max_iter=max_iter_num)
    reg.fit(A,b)
    x = reg.coef_
    np = count_nonzero_vars(x)
    nvars = np

elif algorithm == 'lassocv':
    reg = linear_model.LassoCV(alphas=alpha_ar,fit_intercept=False,cv=folds)
    reg.fit(A,b)
    print '! Lasso CV regression used'
    print '! Lasso CV alpha = ' + str(reg.alpha_)
    x = reg.coef_
    np = count_nonzero_vars(x)
    nvars = np

elif algorithm == 'lassolars':
    print '! LARS implementation of LASSO used'
    print '! LASSO alpha = ', alpha_val
    
    if DO_WEIGHTING:
    	reg = linear_model.LassoLars(alpha=alpha_val,fit_intercept=False,fit_path=False,verbose=True,max_iter=max_iter_num, copy_X=False)
    	reg.fit(weightedA,weightedb)
    else:
    	reg = linear_model.LassoLars(alpha=alpha_val,fit_intercept=False,fit_path=False,verbose=True,max_iter=max_iter_num)
    	reg.fit(A,b)

    x = reg.coef_[0]
    np = count_nonzero_vars(x)
    nvars = np

elif algorithm == 'owlqn':
    x = fit_owlqn(A, b)
    np = count_nonzero_vars(x)
    nvars = np
else:
    print "Unrecognized fitting algorithm" 
    exit(1)

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
print "! number of fitting vars = ", nvars
print "! Bayesian Information Criterion =  ", bic
if args.weights !="None":
    print '! Using weighting file: ',weights_file
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
SNUM_4B = 0

if POTENTIAL == "CHEBYSHEV" or POTENTIAL == "DFTBPOLY":
    TMP = hf[7].split()
    
    if len(TMP) >= 4:
        if len(TMP) >= 5:
            SNUM_4B = int(TMP[4])
            
        SNUM_2B = int(TMP[2])
        
        #print "Expecting:"
        #print SNUM_2B
        #print SNUM_3B
        #print SNUM_4B        

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
    
mapsfile=open(map_file,"r").readlines()

print ""

for i in range(0,len(mapsfile)):
    print mapsfile[i].rstrip('\n')
    
print ""

total_params = TOTAL_PAIRS * SNUM_2B + COUNTED_TRIP_PARAMS + COUNTED_QUAD_PARAMS + COUNTED_COUL_PARAMS + OVERCOORD_PARAMS 
if (total_params != len(x)) and (total_params != (len(x)-1)) :
    sys.stderr.write( "Error in counting parameters")
    exit(1)

if ( total_params + 1 == len(x) ):
    print "ENERGY OFFSET: " + str(x[len(x)-1])
    
if test_suite_run:
    test_suite_params=open("test_suite_params.txt","w")		
    for i in range(0,len(x)):
        phrase = `i` + " " + `x[i]` + '\n'
        test_suite_params.write(phrase)
        test_suite_params.close()

print "ENDFILE"		


# OLD WAY:
#for i in range(0,len(x)):
#    print i,x[i]

                















