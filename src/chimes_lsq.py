import sys
import numpy
import scipy.linalg
import math as m
import subprocess
import os
import argparse

from numpy        import *
from numpy.linalg import lstsq
from numpy.linalg import LinAlgError
from datetime     import *
from subprocess   import call



#############################################
#############################################
# Main
#############################################
#############################################

def main():
    
    loc = os.path.dirname(os.path.abspath(__file__))
    
    
    #############################################
    # Define arguments supported by the lsq code
    #############################################
        
    parser = argparse.ArgumentParser(description='Least-squares force matching based on output of chimes_lsq')

    parser.add_argument("--A",                    type=str,      default='A.txt',         help='A (derivative) matrix') 
    parser.add_argument("--algorithm",            type=str,      default='svd',           help='fitting algorithm')
    parser.add_argument("--dlasso_dlars_path",    type=str,      default=loc+'/../contrib/dlars/src/', help='Path to DLARS and/or DLASSO solver')
    parser.add_argument("--alpha",                type=float,    default=1.0e-04,         help='Lasso regularization')
    parser.add_argument("--b",                    type=str,      default='b.txt',         help='b (force) file')
    parser.add_argument("--cores",                type=int,      default=8,               help='DLARS number of cores')
    parser.add_argument("--eps",                  type=float,    default=1.0e-05,         help='svd regularization')
    parser.add_argument("--fast_svd",             type=str2bool, default=False,           help='solve SVD using normal equation')
    parser.add_argument("--header",               type=str,      default='params.header', help='parameter file header')
    parser.add_argument("--map",                  type=str,      default='ff_groups.map', help='parameter file map')
    parser.add_argument("--nodes",                type=int,      default=1,               help='DLARS number of nodes')
    parser.add_argument("--mpistyle",             type=str,      default="srun",          help='Command used to run an MPI job, e.g. srun, ibrun, mpriun, etc')
    parser.add_argument("--normalize",            type=str2bool, default=False,           help='Normalize DLARS calculation')
    parser.add_argument("--read_output",          type=str2bool, default=False,           help='Read output from previous DLARS run')
    parser.add_argument("--restart_dlasso_dlars", type=str,      default="",              help='Determines whether dlasso or dlars job will be restarted. Argument is the restart file name ')
    parser.add_argument("--split_files",          type=str2bool, default=False,           help='LSQ code has split A matrix output.  Works DLARS.')
    parser.add_argument("--test_suite",           type=str2bool, default=False,           help='output for test suite')
    parser.add_argument("--weights",              type=str,      default="None",          help='weight file')
    parser.add_argument("--active",               type=str2bool, default=False,           help='is this a DLARS/DLASSO run from the active learning driver?')
    parser.add_argument("--folds",                type=int,      default=4,               help="Number of CV folds")
    
    # Actually parse the arguments

    args        = parser.parse_args()
    
    dlasso_dlars_path = args.dlasso_dlars_path

    #############################################
    # Import sklearn modules, if needed
    #############################################

    # Algorithms requiring sklearn.
    sk_algos = ["lasso", "ridge", "lassolars", "lars", "ridgecv"] ;

    if args.algorithm in sk_algos:
        from sklearn import linear_model
        from sklearn import preprocessing
        
    #############################################
    # Read weights, if used
    #############################################        
        
        
    if ( args.weights == "None" ):
        DO_WEIGHTING = False 
    else:
        DO_WEIGHTING = True
        if ( not args.split_files ):
            WEIGHTS= numpy.genfromtxt(args.weights,dtype='float')

    #################################
    #   Process A and b matrices, sanity check weight dimensions
    #################################

    # Use genfromtxt to avoid parsing large files. Note that the AL driver does not use split matrices
    
    if (args.active  and not args.split_files) or ((args.algorithm == "dlasso") and not args.split_files): 
    
        A      = numpy.zeros((1,1),dtype=float)
        b      = numpy.genfromtxt(args.b, dtype='float') 
        np     = "undefined"
        nlines = b.shape[0]

    elif ( (not args.split_files) and (not args.read_output) ) :
        A       = numpy.genfromtxt(args.A , dtype='float')
        nlines  = A.shape[0] 
        np      = A.shape[1] 
        b       = numpy.genfromtxt(args.b, dtype='float') 
        nlines2 = b.shape[0] 

        if ( nlines != nlines2 ):
            print ("Error: the number of lines in the input files do not match\n")
            exit(1) 

            if np > nlines:
                print ("Error: number of variables > number of equations")
                exit(1)
    else:
        
        if not args.read_output:
            dimf = open("dim.0000.txt", "r") ;
            line = next(dimf) 
            dim  = (int(x) for x in line.split())
            A    = numpy.zeros((1,1),dtype=float)           # Dummy A matrix - NOT read in.
            b    = numpy.genfromtxt(args.b, dtype='float')  # Dummy b matrix - NOT read in.
            (np, nstart, nend, nlines) = dim
        else:
            b      = numpy.genfromtxt(args.b, dtype='float') 
            np     = "undefined"
            nlines = b.shape[0]
            
    # Sanity check weight dimensions        
    
    if DO_WEIGHTING and not args.split_files:
        if ( WEIGHTS.shape[0] != nlines ):
            print ("Wrong number of lines in WEIGHTS file")
            exit(1)  

    
    #################################
    # Apply weighting to A and b
    #################################

    weightedA = None
    weightedb = None

    if DO_WEIGHTING and not args.split_files and not args.active  and not (args.algorithm == "dlasso"):

        # This way requires too much memory for long A-mat's
        # to avoid a memory error, we will do it the slow way:

        weightedA = numpy.zeros((A.shape[0],A.shape[1]),dtype=float)
        weightedb = numpy.zeros((A.shape[0],),dtype=float)

        for i in range(A.shape[0]):     # Loop over rows (atom force components)
            for j in range(A.shape[1]): # Loop over cols (variables in fit)
                weightedA[i][j] = A[i][j]*WEIGHTS[i]
                weightedb[i]    = b[i]   *WEIGHTS[i]


    ################################                
    #  Header for output
    ################################
    
    print ("! Date ", date.today())
    print ("!")

    if np != "undefined" :
        print ("! Number of variables            = ", np)

    print ("! Number of equations            = ", nlines)

    
                
    #################################
    # Solve the matrix equation
    #################################

    if args.algorithm == 'svd':
        
        # Modify A and b matrix to reduce A matrix dimension during calculation
        # now 
        if DO_WEIGHTING:
            if args.fast_svd:
                weightedATA = dot(transpose(weightedA), weightedA)
                weightedATb = dot(transpose(weightedA), weightedb)
                eps_sq = args.eps * args.eps
            else:
                pass
        else:
            if args.fast_svd:
                ATA = dot(transpose(A), A)
                ATb = dot(transpose(A), b)
                eps_sq = args.eps * args.eps
            else:
                pass        
        
        # Make the scipy call
        
        print ('! svd algorithm used')
        try:
            if DO_WEIGHTING: # Then it's OK to overwrite weightedA.  It is not used to calculate y (predicted forces) below.
                if args.fast_svd:
                    U,D,VT = scipy.linalg.svd(weightedATA,overwrite_a=True)
                    Dmat   = array((transpose(weightedATA)))
                else:
                    U,D,VT = scipy.linalg.svd(weightedA,overwrite_a=True)
                    Dmat   = array((transpose(weightedA)))
            else:            #  Then do not overwrite A.  It is used to calculate y (predicted forces) below.
                if args.fast_svd:
                    U,D,VT = scipy.linalg.svd(ATA,overwrite_a=False)
                    Dmat   = array((transpose(ATA)))
                else:
                    U,D,VT = scipy.linalg.svd(A,overwrite_a=False)
                    Dmat   = array((transpose(A))) 
        except LinAlgError:
            sys.stderr.write("SVD algorithm failed")
            exit(1)
            
        # Process output

        dmax = 0.0

        for i in range(0,len(Dmat)):
            if ( abs(D[i]) > dmax ) :
                dmax = abs(D[i])

            for j in range(0,len(Dmat[i])):
                Dmat[i][j]=0.0

        # Cut off singular values based on fraction of maximum value as per numerical recipes.
        
        if args.fast_svd:
            eps = eps_sq * dmax
        else:
            eps=args.eps * dmax
        nvars = 0

        for i in range(0,len(D)):
            if abs(D[i]) > eps:
                Dmat[i][i]=1.0/D[i]
                nvars += 1

        print ("! eps (= args.eps*dmax)          =  %11.4e" % eps)        
        print ("! SVD regularization factor      = %11.4e" % args.eps)

        x=dot(transpose(VT),Dmat)

        if DO_WEIGHTING:
            if args.fast_svd:
                x = dot(x,dot(transpose(U),weightedATb))
            else:
                x = dot(x,dot(transpose(U),weightedb))
        else:
            if args.fast_svd:
                x = dot(x,dot(transpose(U),ATb))
            else:
                x = dot(x,dot(transpose(U),b))

    elif args.algorithm == 'ridge':
        print ('! ridge regression used')
        reg = linear_model.Ridge(alpha=args.alpha,fit_intercept=False)

        # Fit the data.
        reg.fit(A,b)

        x = reg.coef_
        nvars = np
        print ("! Ridge alpha = %11.4e" % args.alpha)

    elif args.algorithm == 'ridgecv':
        alpha_ar = [1.0e-06, 3.2e-06, 1.0e-05, 3.2e-05, 1.0e-04, 3.2e-04, 1.0e-03, 3.2e-03]
        reg = linear_model.RidgeCV(alphas=alpha_ar,fit_intercept=False,cv=args.folds)
        reg.fit(A,b)
        print ('! ridge CV regression used')
        print ("! ridge CV alpha = %11.4e"  % reg.alpha_)
        x = reg.coef_
        nvars = np

    elif args.algorithm == 'lasso':
        
        # Make the sklearn call
        
        print ('! Lasso regression used')
        print ('! Lasso alpha = %11.4e' % args.alpha)
        reg   = linear_model.Lasso(alpha=args.alpha,fit_intercept=False,max_iter=100000)
        reg.fit(A,b)
        x     = reg.coef_
        np    = count_nonzero_vars(x)
        nvars = np

    elif args.algorithm == 'lassolars':
        
        # Make the sklearn call
        
        print ('! LARS implementation of LASSO used')
        print ('! LASSO alpha = %11.4e' % args.alpha)

        if DO_WEIGHTING:
            reg = linear_model.LassoLars(alpha=args.alpha,fit_intercept=False,fit_path=False,verbose=True,max_iter=100000, copy_X=False)
            reg.fit(weightedA,weightedb)
        else:
            reg = linear_model.LassoLars(alpha=args.alpha,fit_intercept=False,fit_path=False,verbose=True,max_iter=100000)
            reg.fit(A,b)
        x       = reg.coef_[0]
        np      = count_nonzero_vars(x)
        nvars   = np

    elif args.algorithm == 'dlars' or args.algorithm == 'dlasso' :
        
        # Make the DLARS or DLASSO call

        x,y = fit_dlars(dlasso_dlars_path, args.nodes, args.cores, args.alpha, args.split_files, args.algorithm, args.read_output, args.weights, args.normalize, args.A , args.b ,args.restart_dlasso_dlars, args.mpistyle)
        np = count_nonzero_vars(x)
        nvars = np
        
    else:

        print ("Unrecognized fitting algorithm") 
        exit(1)

    #################################
    # Process output from solver(s)
    #################################

    # If split_files, A is not read in ...This conditional should really be set by the algorithm, since many set  y themselves...  
      
    if ( (not args.split_files) and (not args.read_output) and (not args.active ) and (args.algorithm != "dlasso") ):
        y=dot(A,x)
        
    Z=0.0

    # Put calculated forces in force.txt
    
    yfile = open("force.txt", "w")
    
    for a in range(0,len(b)):
        Z = Z + (y[a] - b[a]) ** 2.0
        yfile.write("%13.6e\n"% y[a]) 

    bic = float(nlines) * log(Z/float(nlines)) + float(nvars) * log(float(nlines))

    #############################################
    # Setup output
    #############################################
    
    print ("! RMS force error                = %11.4e" % sqrt(Z/float(nlines)))
    print ("! max abs variable               = %11.4e" %  max(abs(x)))
    print ("! number of fitting vars         = ", nvars)
    print ("! Bayesian Information Criterion = %11.4e" % bic)
    if args.weights !="None":
        print ('! Using weighting file:            ',args.weights)
    print ("!")

    ####################################
    # Actually process the header file...
    ####################################

    hf = open(args.header ,"r").readlines()
    
    BREAK_COND = False
    
    EXCL_2B = []
    EXCL_1B = []

    # Figure out whether we have triplets and/or quadruplets
    # Find the ATOM_TRIPS_LINE and ATOM_QUADS_LINE
    # Find the TOTAL_TRIPS and TOTAL_QUADS

    ATOM_TRIPS_LINE = 0
    ATOM_QUADS_LINE = 0
    TOTAL_TRIPS = 0
    TOTAL_QUADS = 0

    for i in range(0, len(hf)):
        print (hf[i].rstrip('\n'))
        TEMP = hf[i].split()
        
        if "EXCLD2B" in hf[i]:
            line = hf[i].split()
            if (len(line) == 2) and (line[1] == "false"):
            	EXCL_2B = []
            else:
            	EXCL_2B = list(map(int, line[1:]))	
	    
	    
        if "EXCLD1B" in hf[i]:
            line = hf[i].split()
            if (len(line) == 2) and (line[1] == "false"):
            	EXCL_1B = []
            else:
            	EXCL_1B = list(map(int, line[1:]))	    

        if len(TEMP)>3:
            if (TEMP[2] == "TRIPLETS:"):
                TOTAL_TRIPS = TEMP[3]
                ATOM_TRIPS_LINE = i

                for j in range(i, len(hf)):
                    TEMP = hf[j].split()
                    if len(TEMP)>3:
                        if (TEMP[2] == "QUADRUPLETS:"):
                            print (hf[j].rstrip('\n'))
                            TOTAL_QUADS = TEMP[3]
                            ATOM_QUADS_LINE = j
                            BREAK_COND = True
                            break
            if (BREAK_COND):
                 break

    # 1. Figure out what potential type we have

    POTENTIAL = hf[7].split()
    POTENTIAL = POTENTIAL[1]

    print ("")

    print ("PAIR " + POTENTIAL + " PARAMS \n")

    # 2. Figure out how many coeffs each atom type will have

    SNUM_2B = 0
    SNUM_4B = 0

    if POTENTIAL == "CHEBYSHEV":
        
        TMP = hf[7].split()


        if len(TMP) >= 4:
            if len(TMP) >= 5:
                SNUM_4B = int(TMP[4])

            SNUM_2B = int(TMP[2])  
 

    # 3. Print out the parameters

    FIT_COUL = hf[1].split()
    FIT_COUL = FIT_COUL[1]

    ATOM_TYPES_LINE  = 9
    TOTAL_ATOM_TYPES = hf[ATOM_TYPES_LINE].split()
    TOTAL_ATOM_TYPES = int(TOTAL_ATOM_TYPES[2])
    ATOM_PAIRS_LINE  = ATOM_TYPES_LINE+2+TOTAL_ATOM_TYPES+2
    TOTAL_PAIRS      = hf[ATOM_PAIRS_LINE].split()
    TOTAL_PAIRS      = int(TOTAL_PAIRS[2])
    
    # Remove excluded 2b interactions from accounting
    
    # TOTAL_PAIRS -= len(EXCL_2B) 

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

    if (int(TOTAL_TRIPS) > 0):
        for t in range(0, int(TOTAL_TRIPS)):

            P1 = hf[ATOM_TRIPS_LINE+3+ADD_LINES].split()

            if P1[4] != "EXCLUDED:":
                SNUM_3B +=  int(P1[4])

                TOTL = P1[6]
                ADD_LINES += 5

                for i in range(0, int(TOTL)):
                    ADD_LINES += 1

    # Figure out how many 4B parameters there are

    SNUM_4B   = 0
    ADD_LINES = 0

    if (int(TOTAL_QUADS) > 0):
        for t in range(0, int(TOTAL_QUADS)):

            P1 = hf[ATOM_QUADS_LINE+3+ADD_LINES].split()

            #print "QUAD HEADER", P1
            if P1[7] != "EXCLUDED:":

                SNUM_4B +=  int(P1[7])

                TOTL = P1[9]

                ADD_LINES += 5

                for i in range(0,int(TOTL)):
                    ADD_LINES += 1

    effective_pair = 0

    for i in range(0,TOTAL_PAIRS):

        A1 = hf[ATOM_PAIRS_LINE+2+i+1].split()
        A2 = A1[2]
        A1 = A1[1]

        #print ("PAIRTYPE PARAMS: " + `i` + " " + A1 + " " + A2 + "\n")
        print ("PAIRTYPE PARAMS: " + str(i) + " " + A1 + " " + A2 + "\n")
	
        if i in EXCL_2B:
            for j in range(0, int(SNUM_2B)):
                print ("%3d %21.13e" % (j,0.0))	
        else:
            for j in range(0, int(SNUM_2B)):
                print ("%3d %21.13e" % (j,x[effective_pair*SNUM_2B+j]))
            effective_pair += 1    

        if FIT_COUL == "true":
            if len(EXCL_2B) > 0:
                print("ERROR: cannot use hierarchical learning with charge fitting.")
                exit(0)		
            print ("q_%s x q_%s %21.13e" % (A1,A2,x[TOTAL_PAIRS*SNUM_2B + SNUM_3B + SNUM_4B + i]))
            COUNTED_COUL_PARAMS += 1

        print (" ")

    # TRIPLETS

    ADD_LINES = 0
    ADD_PARAM = 0

    COUNTED_TRIP_PARAMS = 0

    if (int(TOTAL_TRIPS) > 0):
        print ("TRIPLET " + POTENTIAL + " PARAMS \n")

        TRIP_PAR_IDX = 0

        for t in range(0, int(TOTAL_TRIPS)):

            PREV_TRIPIDX = 0

            print ("TRIPLETTYPE PARAMS:")
            print ("  " + hf[ATOM_TRIPS_LINE+2+ADD_LINES].rstrip())

            P1 = hf[ATOM_TRIPS_LINE+3+ADD_LINES].split()

            #print "HEADER: ", P1

            V0 = P1[1] 
            V1 = P1[2]
            V2 = P1[3]

            if P1[4] == "EXCLUDED:" :
                print ("   PAIRS: " + V0 + " " + V1 + " " + V2 + " EXCLUDED:")
                ADD_LINES += 1
            else:
                UNIQ = P1[4]
                TOTL = P1[6].rstrip() 

                print ("   PAIRS: " + V0 + " " + V1 + " " + V2 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL)
                print ("     index  |  powers  |  equiv index  |  param index  |       parameter       ")
                print ("   ----------------------------------------------------------------------------")

                ADD_LINES += 3

                if(t>0):
                    ADD_PARAM += 1

                for i in range(0,int(TOTL)):
                    ADD_LINES += 1
                    LINE       = hf[ATOM_TRIPS_LINE+2+ADD_LINES].rstrip('\n')
                    LINE_SPLIT = LINE.split()

                    print ("%s %21.13e" % (LINE, x[effective_pair*SNUM_2B + TRIP_PAR_IDX+int(LINE_SPLIT[5])]))

                TRIP_PAR_IDX += int(UNIQ)
                COUNTED_TRIP_PARAMS += int(UNIQ)
                #print "COUNTED_TRIP_PARAMS", COUNTED_TRIP_PARAMS

            print ("")

            ADD_LINES += 2

    ADD_LINES = 0
    
    # QUADS    

    COUNTED_QUAD_PARAMS = 0
    if (int(TOTAL_QUADS) > 0):
        print ("QUADRUPLET " + POTENTIAL + " PARAMS \n")

        QUAD_PAR_IDX = 0

        for t in range(int(TOTAL_QUADS)):

            PREV_QUADIDX = 0

            #print "ATOM_QUADS_LINE " + str(ATOM_QUADS_LINE+2+ADD_LINES)

            P1 = hf[ATOM_QUADS_LINE+2+ADD_LINES].split()

            #print "P1 " + P1[1] + P1[2] + P1[3] + P1[4] + P1[5] + P1[6]

            print ("QUADRUPLETYPE PARAMS: " )
            print ("  " + hf[ATOM_QUADS_LINE+2+ADD_LINES].rstrip() )

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
                print ("   PAIRS: " + V0 + " " + V1 + " " + V2 + " " + V3 + " " + V4 + " " + V5 + " EXCLUDED: ")
                ADD_LINES += 1

            else:
                UNIQ = P1[7]
                TOTL = P1[9].rstrip() 

                print ("   PAIRS: " + V0 + " " + V1 + " " + V2 + " " + V3 + " " + V4 + " " + V5 + " UNIQUE: " + UNIQ + " TOTAL: " + TOTL)
                print ("     index  |  powers  |  equiv index  |  param index  |       parameter       ")
                print ("   ----------------------------------------------------------------------------")

                ADD_LINES += 3

                if(t>0):
                    ADD_PARAM += 1

                for i in range(0,int(TOTL)):
                    ADD_LINES += 1
                    LINE       = hf[ATOM_QUADS_LINE+2+ADD_LINES].rstrip('\n')
                    LINE_SPLIT = LINE.split()

                    UNIQ_QUAD_IDX = int(LINE_SPLIT[8])
                    #print 'UNIQ_QUAD_IDX', str(UNIQ_QUAD_IDX)

                    print ("%s %21.13e" % (LINE,x[effective_pair*SNUM_2B + COUNTED_TRIP_PARAMS + QUAD_PAR_IDX + UNIQ_QUAD_IDX]))

                QUAD_PAR_IDX += int(UNIQ)
                COUNTED_QUAD_PARAMS += int(UNIQ)

            print ("")

            ADD_LINES += 2

    # Remaining tidbids

    mapsfile=open(args.map,"r").readlines()

    print ("")

    for i in range(0,len(mapsfile)):
        print (mapsfile[i].rstrip('\n'))

    print ("")

    total_params = effective_pair * SNUM_2B + COUNTED_TRIP_PARAMS + COUNTED_QUAD_PARAMS + COUNTED_COUL_PARAMS 
    
    N_ENER_OFFSETS = int(hf[9].split()[2]) - len(EXCL_1B)
    N_ATOM_TYPES   = int(hf[9].split()[2])




## Parameter count could be off by natom_types, if energies are included in the fit
    if (total_params != len(x)) and (len(x) != (total_params+N_ENER_OFFSETS)) :
        sys.stderr.write( "Error in counting parameters\n") 
        sys.stderr.write("total_params " + str(total_params) + "\n")
        sys.stderr.write("N_ENER_OFFSETS " + str(N_ENER_OFFSETS) + "\n")
        sys.stderr.write("total_params+N_ENER_OFFSETS " + str(total_params+N_ENER_OFFSETS) + "\n")	
        sys.stderr.write("len(x) " + str(len(x)) + "\n") 
        sys.stderr.write("TOTAL_PAIRS " + str(TOTAL_PAIRS) + "\n") 
        sys.stderr.write("SNUM_2B " + str(SNUM_2B) + "\n") 
        sys.stderr.write("COUNTED_TRIP_PARAMS " + str(COUNTED_TRIP_PARAMS) + "\n") 
        sys.stderr.write("COUNTED_QUAD_PARAMS " + str(COUNTED_QUAD_PARAMS) + "\n")
        sys.stderr.write("COUNTED_COUL_PARAMS " + str(COUNTED_COUL_PARAMS) + "\n")
        sys.stderr.write("============= ")
        for i in range(10):
                sys.stderr.write(str(i) + "	" + hf[i])
        exit(1)


    if len(x) == (total_params+N_ENER_OFFSETS):
        print ("NO ENERGY OFFSETS: ", N_ATOM_TYPES)
	
        eff_idx = 0
    
        for i in range(N_ATOM_TYPES):
	
            if i in EXCL_1B:
                print ("ENERGY OFFSET %d %21.13e" % (i+1,0.0))	    
            else:
                print ("ENERGY OFFSET %d %21.13e" % (i+1,x[total_params+eff_idx]))
                eff_idx += 1

    if args.test_suite:
        test_suite_params=open("test_suite_params.txt","w")		
        for i in range(0,len(x)):
            phrase = "%5d %21.13e\n" % (i,x[i])
            test_suite_params.write(phrase)
        test_suite_params.close()

    print ("ENDFILE")
    return 0

#############################################
#############################################
# Small helper functions
#############################################
#############################################

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
                       

def count_nonzero_vars(x):
    np = 0
    for i in range(0, len(x)):
        if ( abs(x[i]) > 1.0e-05 ):
            np = np + 1
    return np


#############################################
#############################################
# DLARS wrapper
#############################################
#############################################

def fit_dlars(dlasso_dlars_path, nodes, cores, alpha, split_files, algorithm, read_output, weights, normalize, A , b, restart_dlasso_dlars, mpistyle):

    # Use the Distributed LARS/LASSO fitting algorithm.  Returns both the solution x and
    # the estimated force vector A * x, which is read from Ax.txt.    
    
    if dlasso_dlars_path == '':
        print ("ERROR: DLARS/DLASSO  path not provided.")
        print ("Please run again with --dlasso_dlars_path </absolute/path/to/dlars/dlasso/src/>")
        exit(0)

    if algorithm == 'dlasso' :
        print ('! DLARS code for LASSO used')
    elif algorithm == 'dlars' :
        print ('! DLARS code for LARS used')
    else:
        print ("Bad algorithm in fit_dlars:" + algorithm)
        exit(1)
    print ('! DLARS alpha = %10.4e' % alpha)

    if not read_output:
    
        dlars_file = dlasso_dlars_path + "dlars"
        
        if os.path.exists(dlars_file):
	
            exepath = ""

            if mpistyle == "srun": 
                exepath = "srun -N " + str(nodes) + " -n " + str(cores) + " " + dlars_file
            elif mpistyle == "ibrun":
                exepath = "ibrun" + " " + dlars_file  
            else:
                print("Unrecognized mpistyle:",args.mpistyle,". Recognized options are srun or ibrun")
           
            command = None

            command = exepath + " " + A  + " " + b + " dim.txt --lambda=" + str(alpha)

            #else:
            #    command = ("{0} A.txt b.txt dim.txt --lambda={1}".format(exepath, alpha))

            if ( split_files ) :
                command = command + " --split_files"
            if ( algorithm == 'dlars' ):
                command = command + " --algorithm=lars"
            elif ( algorithm == 'dlasso' ):
                command = command + " --algorithm=lasso"

            if ( weights != 'None' ):
                command = command + " --weights=" + weights

            if ( normalize ):
                command = command + " --normalize=y" 
            else:
                command = command + " --normalize=n" 

            if restart_dlasso_dlars != "":
                print ("Will run a dlars/dlasso restart job with file:", restart_dlasso_dlars)

                command = command + " --restart=" + restart_dlasso_dlars
                
            command = command +  " >& dlars.log"

            print("! DLARS run: " + command + "\n")

            if ( os.system(command) != 0 ) :
                print(command + " failed")
                sys.exit(1)
        else:
            print (dlars_file + " does not exist")
            sys.exit(1)
    else:
        print ("! Reading output from prior DLARS calculation")

    x = numpy.genfromtxt("x.txt", dtype='float')
    y = numpy.genfromtxt("Ax.txt", dtype='float') 
    return x,y



# Python magic to allow having a main function definition.    
if __name__ == "__main__":
    main()
    














