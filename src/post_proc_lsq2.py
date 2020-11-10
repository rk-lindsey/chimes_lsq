import sys

# Parameter file post-processing .. remove zeroed-out parameters
#
# Usage: python <this script> <parameter file> 
#
# To do:
#
# 1. allow zeroing of 2B parameters
# 2. Account for the case where every interaction is set to zero


class CLU_TYPE:
    
    def __init__(self,TYPE):
        
        self.NUM_PAR_COL = 0
        self.UNQ_PAR_COL = 0
        self.PAR_VAL_COL = 0
        
        if TYPE == "PAIRTYPE":
            print "ERROR: This class shouldn't encounter \"PAIRTYPE\". Exiting."
            exit()            

        elif TYPE == "TRIPLETTYPE":
            
            self.NUM_PAR_COL = 7
            self.UNQ_PAR_COL = 5
            self.PAR_VAL_COL = 6

        elif TYPE == "QUADRUPLETYPE":
            
            self.NUM_PAR_COL = 10
            self.UNQ_PAR_COL = 8
            self.PAR_VAL_COL = 9

        else: # Undefined type.. probably a custer of greater than 4 atoms
            
            print "ERROR: Cluster type " + TYPE + " undefined. Exiting."
            exit()

INFILE  = open(sys.argv[1],'r')  # A parameter file produced by lsq2.py script
OUTFILE = open(sys.argv[1]+".reduced",'w') 

# Count up number of TOTAL remaining parameters after reduction, update the parameter list accordingly

FINISHED = False

while not FINISHED:
    
    PARSE = INFILE.readline()
    
    if len(PARSE.split())>0 and PARSE.split()[0] == "ENDFILE":
        
        OUTFILE.write(PARSE)        
        FINISHED = True
        break
        
    if len(PARSE.split())==2 and PARSE.split()[1] == "PARAMS:" and PARSE.split()[0] != "PAIRTYPE":
        
        CLUSTER = CLU_TYPE(PARSE.split()[0]) # Is this a pair, triplet, quad, etc
        
        OUTFILE.write(PARSE)        # "TRIPLETTYPE PARAMS:"
        
        PARSE = INFILE.readline() 
        
        OUTFILE.write(PARSE)        # "INDEX: 0 ATOMS: Al Al Al"
        
        NO_PARAM_LINE  = INFILE.readline().split()
        NUM_PARAMS = int(NO_PARAM_LINE[CLUSTER.NUM_PAR_COL])
        
        HOLD_1 = INFILE.readline()
        HOLD_2 = INFILE.readline()
        
        TMP_PARAMS = []
        
        # Read triple type parameters, decide which will remain
        
        for i in xrange(NUM_PARAMS): 
            
            PARSE = INFILE.readline()
            
            if abs(float(PARSE.split()[CLUSTER.PAR_VAL_COL])) > 0.000001: # Parameter NOT ignored
                
                TMP_PARAMS.append(PARSE)
                
        # Re-index triple params. Note: equiv index and param index aren't actually used by the MD code,
        # so they don't need to be re-numbered
        
        for i in xrange(len(TMP_PARAMS)):
            
            TMP_PARAMS[i] = TMP_PARAMS[i].split()
            TMP_PARAMS[i][0] = `i`
            
            ' '.join(TMP_PARAMS[i])
        
        # Update parameter counts
        
        NO_PARAM_LINE[CLUSTER.UNQ_PAR_COL] = "-1"
        NO_PARAM_LINE[CLUSTER.NUM_PAR_COL] = `len(TMP_PARAMS)`
        
        ' '.join(NO_PARAM_LINE)
        
        # Now print out the updated parameter file
        
        OUTFILE.write(' '.join(NO_PARAM_LINE)+'\n')
        OUTFILE.write(HOLD_1)
        OUTFILE.write(HOLD_2)
        
        for i in xrange(len(TMP_PARAMS)):
            OUTFILE.write(' '.join(TMP_PARAMS[i])+'\n')

    else:
        OUTFILE.write(PARSE)    

INFILE.close()
OUTFILE.close()
