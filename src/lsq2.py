
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

def main():
    parser = argparse.ArgumentParser(description='Least-squares force matching based on output of chimes_lsq')
    # Arguments supported by the lsq code.
    # Decided to use underscore rather than - for parameter word separator character.
    parser.add_argument("--A", default='A.txt',help='A (derivative) matrix') 
    parser.add_argument("--algorithm", default='svd', help='fitting algorithm')
    parser.add_argument("--alpha", type=float, default=1.0e-04,help='Lasso or ridge regularization')
    parser.add_argument("--b", default='b.txt',help='b (force) file')
    parser.add_argument("--beta",type=float, default=0.0, help='DOWLQN L2 regularization')
    parser.add_argument("--cores", type=int, default=8, help='DOWLQN/DLARS number of cores')
    parser.add_argument("--eps", type=float, default=1.0e-05,help='svd regularization')
    parser.add_argument("--folds",type=int, default=4,help="Number of CV folds")
    parser.add_argument("--header", default='params.header',help='parameter file header')
    parser.add_argument("--map", default='ff_groups.map',help='parameter file map')
    parser.add_argument("--memory", type=int, default=20, help='OWLQN or DOWLQN Hessian memory storage') ;
    parser.add_argument("--nodes", type=int, default=1, help='DOWLQN/DLARS number of nodes')
    parser.add_argument("--normalize", type=str2bool, default=False, help='Normalize DLARS calculation')
    parser.add_argument("--read_output", type=str2bool, default=False, help='Read output from previous DLARS run')
    parser.add_argument("--restart", type=str2bool, default=False, help='Use DOWLQN restart file')
    parser.add_argument("--split_files",    type=str2bool, default=False, help='LSQ code has split A matrix output.  Works with DOWLQN or DLARS.')
    parser.add_argument("--test_suite", type=str2bool, default=False,help='output for test suite')
    parser.add_argument("--tol", type=float, default=1.0e-05, help='OWLQN or DOWLQN tolerance')
    parser.add_argument("--weights", default="None",help='weight file')
    parser.add_argument("--active", default=False, help='is this a DLARS/DLASSO run from the active learning driver?')

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

    # Do we restart calculations (currently only for DOWLQN)
    restart = args.restart

    # Do we read output from a previous DLARS run ?
    read_output = args.read_output
    
    if ( weights_file == "None" ):
        DO_WEIGHTING = False 
    else:
        DO_WEIGHTING = True
	
    is_active = args.active

    # Algorithms requiring sklearn.
    sk_algos = ["lasso", "ridge", "lassolars", "lars", "lasso_split", "ridgecv"] ;
        
    if algorithm in sk_algos:
        from sklearn import linear_model
        from sklearn import preprocessing

    # Regularization values
    alpha_val = args.alpha
    beta_val  = args.beta
    
    # Tolerance and memory for OWLQN/DOWLQN
    tol = args.tol
    memory = args.memory

    # Parameters for parallel run.
    num_cores = args.cores
    num_nodes = args.nodes
    split_files = args.split_files

    normalize = args.normalize

    alpha_ar = [1.0e-06, 3.2e-06, 1.0e-05, 3.2e-05, 1.0e-04, 3.2e-04, 1.0e-03, 3.2e-03]

    # Number of folds for cross-validation.
    folds = args.folds

    max_iter_num = 100000
    
    test_suite_run = False

    # Then this run is being used for the test suite... 
    # print out parameters without any fanciness so tolerances can be checked  
    
    if ( DO_WEIGHTING and not split_files ):
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
    
    if (is_active and not split_files): # Then this is an AL driver run with dlars/dlasso ... we do NOT split A-matrix
    
        A      = numpy.zeros((1,1),dtype=float)
        b      = numpy.genfromtxt(bfile, dtype='float') 
	np     = "undefined"
	nlines = b.shape[0]

    elif ( (not split_files) and (not read_output) ) :
        A = numpy.genfromtxt(Afile, dtype='float')
        nlines = A.shape[0] 
        np = A.shape[1] 

        b = numpy.genfromtxt(bfile, dtype='float') 
        nlines2 = b.shape[0] 

        if ( nlines != nlines2 ):
            print "Error: the number of lines in the input files do not match\n"
            exit(1) 

            if np > nlines:
                print "Error: number of variables > number of equations"
                exit(1)

    else:
        dimf = open("dim.0000.txt", "r") ;
        line = next(dimf) 
        dim = (int(x) for x in line.split())
        (np, nstart, nend, nlines) = dim
        # Dummy A and b matrices.  NOT read in.
        A = numpy.zeros((1,1),dtype=float)
        b = numpy.genfromtxt(bfile, dtype='float') 
        
    hf = open(header_file,"r").readlines()

    print "! Date ", date.today() 
    print "!"
    print "! Number of variables = ", np
    print "! Number of equations    = ", nlines

    #A=zeros((nlines,np))
    #b=zeros((nlines))

    if DO_WEIGHTING and not split_files:
        if ( WEIGHTS.shape[0] != nlines ):
            print "Wrong number of lines in WEIGHTS file"
            exit(1)

    #################################
    # Apply weighting to A and b
    #################################

    weightedA = None
    weightedb = None

    if DO_WEIGHTING and not split_files and not is_active:

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

    elif algorithm == 'lasso_split':
        nfit = 4
        x = fit_lasso_split(alpha_val, A, b, nfit)
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
        x = fit_owlqn(A, b, alpha_val, beta_val, tol,memory)
        np = count_nonzero_vars(x)
        nvars = np

    elif algorithm == 'dowlqn':
        # Note: x and y are both set by fit_dowlqn
        x,y = fit_dowlqn(A, b, num_nodes, num_cores, alpha_val, beta_val, tol, memory, split_files, restart, weights_file)
        np = count_nonzero_vars(x)
        nvars = np

    elif algorithm == 'dlars' or algorithm == 'dlasso' :

        x,y = fit_dlars(num_nodes, num_cores, alpha_val, split_files, algorithm, read_output, weights_file, normalize, Afile, bfile, is_active)
        np = count_nonzero_vars(x)
        nvars = np
        
    else:

        print "Unrecognized fitting algorithm" 
        exit(1)

    # If split_files, A is not read in.        
    # This conditional should really be set by the algorithm, since many set  y themselves...    
    if ( (not split_files) and (not read_output) and (not is_active) ) :
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

    N_ENER_OFFSETS = int(hf[9].split()[2])

## Parameter count could be off by natom_types, if energies are included in the fit
    if (total_params != len(x)) and (len(x) != (total_params+N_ENER_OFFSETS)) :
        sys.stderr.write( "Error in counting parameters\n") 
        sys.stderr.write("len(x) " + str(len(x)) + "\n") 
        sys.stderr.write("TOTAL_PAIRS " + str(TOTAL_PAIRS) + "\n") 
        sys.stderr.write("SNUM_2B " + str(SNUM_2B) + "\n") 
        sys.stderr.write("COUNTED_TRIP_PARAMS " + str(COUNTED_TRIP_PARAMS) + "\n") 
        sys.stderr.write("COUNTED_QUAD_PARAMS " + str(COUNTED_QUAD_PARAMS) + "\n")
        sys.stderr.write("COUNTED_COUL_PARAMS " + str(COUNTED_COUL_PARAMS) + "\n")
        sys.stderr.write("OVERCOORD_PARAMS " + str(OVERCOORD_PARAMS) + "\n")
        exit(1)


    if len(x) == (total_params+N_ENER_OFFSETS):
        print "NO ENERGY OFFSETS: ", N_ENER_OFFSETS
    
        for i in xrange(N_ENER_OFFSETS):
            print "ENERGY OFFSET " + `i+1` + " " + str(x[total_params+i])

    if test_suite_run:
        test_suite_params=open("test_suite_params.txt","w")		
        for i in range(0,len(x)):
            phrase = `i` + " " + `x[i]` + '\n'
            test_suite_params.write(phrase)
            test_suite_params.close()

    print "ENDFILE"		
    return 0

    # OLD WAY:
    #for i in range(0,len(x)):
    #    print i,x[i]


def is_number(s):
# Test if s is a number.
    try:
        float(s)
        return True
    except ValueError:
        return False


def str2bool(v):
## Convert string to bool.  Used in command argument parsing.
    if v.lower() in ('yes', 'true', 't', 'y'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n'):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected")
                       
        
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
            ff.write("%21.16e\n" % (mat[i]))
    elif len(mat.shape) == 2:
        ff.write(str(mat.shape[0]) + " " + str(mat.shape[1]) + "\n")
        for j in range(0, mat.shape[1]):
            for i in range(0, mat.shape[0]):
                ff.write("%21.16e\n" % (mat[i,j]))


def read_matrix_market(f):
# Read an array in matrix market format.
# Matrix market uses column-oriented storage.
# Only dense 1 and 2d arrays supported.    
    ff = open(f, "r")
    dat=ff.readlines()

    sz=dat[1].split()
    szi = [int(sz[0]), int(sz[1])]

    #print "ARRAY SIZE = " + str(szi[0]) + " " + str(szi[1])

    if szi[0] == 1:
        ar = numpy.zeros(szi[1], dtype=float)
        #print "SHAPE: " + str(ar.shape[0])
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

def fit_owlqn(A,b,alpha_val,beta_val,tol,memory):
## Use the OWLQN fitting algorithm
    print '! OWLQN algorithm for LASSO used'
    print '! OWLQN alpha = ' + str(alpha_val)
    write_matrix_market(A, 'Amm.txt')
    write_matrix_market(b, 'bmm.txt')
    path=os.path.dirname(os.path.abspath(__file__)) # "sys.path[0]" converts to lowercase. We need the actual path casing
    exepath=path[:-3] + "contrib/owlqn/source/owlqn"
    if os.path.exists(exepath):
        command = ( exepath
        + " Amm.txt bmm.txt "
        + str(alpha_val)
        + " fit.txt -ls -tol "
        + str(tol)
        + " -m "
        + str(memory)
        + " >& owlqn.log" )
        print("! Running " + command + "\n")

        os.system(command)
        x = read_matrix_market("fit.txt")
        return x
    else:
        sys.stderr.write(exepath + " does not exist\n");
        sys.exit(1)


def fit_dowlqn(A,b,num_nodes, num_cores, alpha_val, beta_val, tol, memory, split_files, restart, weights_file):
## Use the Distributed OWLQN fitting algorithm.  Returns both the solution x and
## the estimated force vector A * x, which is read from Ax.txt.    
    print '! DOWLQN algorithm for LASSO used'
    print '! DOWLQN alpha = ' + str(alpha_val)
    path=os.path.dirname(os.path.abspath(__file__))
    dowlqn_file = path[:-3] + "contrib/owlqn/mpi/dowlqn"
    exepath = "srun -N " + str(num_nodes) + " -n " + str(num_cores) + " "
    exepath = exepath + dowlqn_file

    if ( weights_files != "None" ):
        print "Error: dowlqn does not currently support weighting"
        exit(1)
        
    if os.path.exists(dowlqn_file):
        if ( split_files ) :
            command = ("{0} A.txt b.txt {1} fit.txt -ls -tol {2} -m {3} -l2weight {4} -s".format(exepath, alpha_val, tol, memory, beta_val))
        else:
            write_matrix_market(A, 'Amm.txt')
            write_matrix_market(b, 'bmm.txt')
            command = ("{0} Amm.txt bmm.txt {1} fit.txt -ls -tol {2} -m {3} -l2weight {4}".format(exepath, alpha_val, tol, memory, beta_val))
        if ( restart ) :
            command = command + " -init restart.txt "
        command = command +  " >& dowlqn.log"
        print("! DOWLQN run: " + command + "\n")

        if ( os.system(command) != 0 ) :
            print(command + " failed")
            sys.exit(1)
            
        x = read_matrix_market("fit.txt")
        y = read_matrix_market("Ax.txt")
        return x,y
    else:
        print exepath + " does not exist"
        sys.exit(1)


def fit_dlars(num_nodes, num_cores, alpha_val, split_files, algorithm, read_output, weights_file, normalize, Afile, bfile, is_active):

    # Use the Distributed LARS/LASSO fitting algorithm.  Returns both the solution x and
    # the estimated force vector A * x, which is read from Ax.txt.    

    if algorithm == 'dlasso' :
        print '! DLARS code for LASSO used'
    elif algorithm == 'dlars' :
        print '! DLARS code for LARS used'
    else:
        print "Bad algorithm in fit_dlars:" + algorithm
        exit(1)
    print '! DLARS alpha = ' + str(alpha_val)

    if not read_output:
    
    	# Use the following if lsq2 is in a full ChIMES directory

        path=os.path.dirname(os.path.abspath(__file__))
        dlars_file = path[:-3] + "contrib/dlars/src/dlars"
	
	# Otherwise, use the following to hard-code the path:
	
	#dlars_file = "/my/hardcoded/path/"
	
        exepath = "srun -N " + str(num_nodes) + " -n " + str(num_cores) + " "
	
        exepath = exepath + dlars_file
	
        if os.path.exists(dlars_file):
	
            command = None
	        
            if  is_active: # Then we're using the parallel driver and files are named differently

		command = exepath + " " + Afile + " " + bfile + " dim.txt --lambda=" + `alpha_val`
		normalize = True

            else:
                command = ("{0} A.txt b.txt dim.txt --lambda={1}".format(exepath, alpha_val))

            if ( split_files ) :
                command = command + " --split_files"
            if ( algorithm == 'dlars' ):
                command = command + " --algorithm=lars"
            elif ( algorithm == 'dlasso' ):
                command = command + " --algorithm=lasso"

            if ( weights_file != 'None' ):
                command = command + " --weights=" + weights_file

            if ( normalize ):
                command = command + " --normalize=y" 
            else:
                command = command + " --normalize=n" 
                
            command = command +  " >& dlars.log"
	    
            print("! DLARS run: " + command + "\n")

            if ( os.system(command) != 0 ) :
                print(command + " failed")
                sys.exit(1)
        else:
            print exepath + " does not exist"
            sys.exit(1)
    else:
        print "! Reading output from prior DLARS calculation"

    x = numpy.genfromtxt("x.txt", dtype='float')
    y = numpy.genfromtxt("Ax.txt", dtype='float') 
    return x,y

def fit_lasso_split(alpha_val, A, b, nsplit):
## Lasso Regression based on splitting the data n ways and averaging
## the coefficients.    
## This doesn't work well !    
    print '! Lasso Split Data regression used'
    print '! Lasso alpha = ' + str(alpha_val)

    step = A.shape[0]/nsplit
    print '! Dividing data into ' + str(nsplit) + ' parts'
    xavg = numpy.zeros(A.shape[1])
    for i in range(0,nsplit):
        if ( i != nsplit - 1 ) :
            Asplit = A[i * step: (i+1) * step]
            bsplit = b[i * step: (i+1) * step]
        else:
            Asplit = A[i * step: ]
            bsplit = b[i * step: ]
        print "bsplit shape = "
        print bsplit.shape
        print "Asplit shape = "
        print Asplit.shape

    	reg = linear_model.LassoLars(alpha=alpha_val,fit_intercept=False,fit_path=False,verbose=True,max_iter=max_iter_num)
        reg.fit(Asplit,bsplit)
        print reg.coef_[0]
        x = reg.coef_[0]
        print "X " + str(i) + "="
        print x
        xavg = xavg + x
    xavg = divide(xavg, nsplit)
    return xavg
    

# Python magic to allow having a main function definition.    
if __name__ == "__main__":
    main()
    














