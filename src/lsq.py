#! /usr/bin/python

import sys
import numpy
from numpy import *
##sys.path.append('/g/g90/koziol3/codes/mol/')
##import matrix
from numpy.linalg import lstsq
from datetime import *

TEST_SUITE_RUN = False

if(len(sys.argv) < 4):
    print "The following args were supplied:"
    print sys.argv
    
    print "usage is:"
    print "./lsq.py A.txt b.txt params.header"
    sys.exit()
elif( (len(sys.argv) == 5) and (sys.argv[4]=="TEST_SUITE_RUN")):
    # Then this run is being used for the test suite... 
    # print out parameters without any fanciness so tolerances can be checked  
    TEST_SUITE_RUN = "do"


#################################

###################################



af=open(sys.argv[1],"r").readlines()
bf=open(sys.argv[2],"r").readlines()
hf=open(sys.argv[3],"r").readlines() ;

nlines=len(af)
nlines2=len(bf)
if ( nlines != nlines2 ):
    print "Error: the number of lines in the input files do not match\n"
    exit(1) 

np=len(af[0].split())
print "# Date ", date.today() ;
print "# Number of columns = ", np
print "# Number of rows    = ", nlines

A=zeros((nlines,np))
b=zeros((nlines))

for i in range(0,nlines):
    af[i]=af[i].split()
    for n in range(0,np):
#    	print `len(A)` + " " + `len(A[i])`+ `len(af)` + " " + `len(af[i])` + " " + `i` + " " + `n` 
        A[i][n]=float(af[i][n])
    b[i]=float(bf[i])

####
#print "finished reading."
#print "Number of configurations", (nlines-2)/3
#for i in range(0,(nlines-2)/3):
#    k = 3*i
#    scale = sqrt(b[k]*b[k] + b[k+1]*b[k+1] + b[k+2]*b[k+2])+10.0
#    b[k]=b[k]/scale
#    b[k+1]=b[k+1]/scale
#    b[k+2]=b[k+2]/scale
#    for n in range(0,np):
#        A[k][n]=A[k][n]/scale
#        A[k+1][n]=A[k+1][n]/scale
#        A[k+2][n]=A[k+2][n]/scale


####SVD
U,D,VT=numpy.linalg.svd(A)


Dmat=array((transpose(A)))



dmax = 0.0
for i in range(0,len(Dmat)):
    if ( abs(D[i]) > dmax ) :
        dmax = abs(D[i])
    for j in range(0,len(Dmat[i])):
        Dmat[i][j]=0.0

# Cut off singular values based on fraction of maximum value as per 
# numerical recipes (LEF).
eps=1.0e-5 * dmax

for i in range(0,len(D)):
    if(abs(D[i])>eps):
        Dmat[i][i]=1.0/D[i]


x=dot(transpose(VT),Dmat)
x=dot(x,dot(transpose(U),b))


#print "SHAPE OF Dmat: ***** " + `shape(Dmat)`
#print "SHAPE OF VT  : ***** " + `shape(VT)`
#print "SHAPE OF U   : ***** " + `shape(U)`
#print "SHAPE OF b   : ***** " + `shape(b)`

y=dot(A,x)
Z=0.0

# Put calculated forces in force.txt
yfile = open("force.txt", "w")
for a in range(0,len(b)):
    Z=Z+(y[a]-b[a])**2.0
    yfile.write("%13.6e\n"% y[a]) 

print "# RMS force error = " , sqrt(Z/float(nlines))

print "# max variable = ",  max(x)

#########################################################################################################################################################################
# Now, instead of just printing the header file, we need to parse the info so it can be printed in a format consistent with the old version of the md cod


#for i in range(0, len(hf)):
#    sys.stdout.write(hf[i])

# VARS WE NEED TO POPULATE:

pair_type         = hf[6].split()
pair_type         = pair_type[1]


# 2. Figure out how many coeffs each atom type will have

CASE_BY_CASE = False

SNUM_2B = 0
SNUM_3B = 0

if pair_type == "CHEBYSHEV" or pair_type == "DFTBPOLY":
	TMP = hf[6].split()
	
	if len(TMP) == 4:
		SNUM_3B = int(TMP[3])
	
	SNUM_2B = int(TMP[2])
	
elif pair_type == "INVRSE_R":
	TMP = hf[6].split()
	SNUM_2B = int(TMP[2])
else:
	CASE_BY_CASE = True # We'll need to do it per atom pair type.
	
coulomb = hf[0].split()
coulomb = coulomb[1]
	
fit_coulomb = hf[1].split()
fit_coulomb = fit_coulomb[1]	

fit_pover = hf[3].split()
fit_pover = fit_pover[1]

nover=5

overcoord = hf[2].split()
overcoord = overcoord[1]

ATOM_TYPES_LINE=8

TOTAL_ATOM_TYPES = hf[ATOM_TYPES_LINE].split()

TOTAL_ATOM_TYPES = int(TOTAL_ATOM_TYPES[2])

ATOM_PAIRS_LINE=10+TOTAL_ATOM_TYPES+2

npair =  hf[ATOM_PAIRS_LINE].split()
npair = int(npair[2])

threeb_cheby = hf[4].split()
threeb_cheby = threeb_cheby[1]	

pair_type_params  = []
overcoord_params  = []

for i in range(0,npair):
	
	MIN = hf[ATOM_PAIRS_LINE+2+i+1].split()
	MAX = float(MIN[4])
	DEL = float(MIN[5])
	MIN = float(MIN[3])	
	
	if CASE_BY_CASE:  
		
		# Figure out individual SNUM_2B
		
		SNUM_2B = int((2+m.floor((MAX - MIN)/DEL))*2)
		
		PAR_LINE = `MIN` + " "  + `MAX` + " " + `DEL` + " " + `SNUM_2B`
	elif pair_type == "CHEBYSHEV" > 0:
		PAR_LINE = `MIN` + " "  + `MAX` + " " + `DEL` + " " + `SNUM_2B` + " " + `SNUM_3B`
		pair_type_params.append(PAR_LINE)
	else: 
		PAR_LINE = `MIN` + " "  + `MAX` + " " + `DEL` + " " + `SNUM_2B`
	pair_type_params.append(PAR_LINE)

if(overcoord == "true"):
	for i in range(0,npair):
		POVER = hf[ATOM_PAIRS_LINE+npair+2+2+i+1].split()
		if(POVER[3] == "1"):
			overcoord_params.append(POVER[5])
			overcoord_params.append(POVER[6])
			overcoord_params.append(POVER[7])
			overcoord_params.append(POVER[8])
			overcoord_params.append(POVER[9])
			break
		


if pair_type == "SPLINE":
	pair_type = spline
elif pair_type == "CHEBYSHEV":
	pair_type = chebyshev
elif pair_type == "DFTBPOLY":
	pair_type = dftbpoly
elif pair_type == "INVRSE_R":
	pair_type = inverse_r


print "npair " + `npair`
print "3b_cheby " + threeb_cheby
print "pair_type " + pair_type.lower()
for i in range(0,npair):
	print pair_type_params[i]
print "coulomb " + coulomb
print "fit_coulomb " + fit_coulomb
print "overcoord " + overcoord
if(overcoord == "true"):
	print "nover 5"
	print overcoord_params[0]
	print overcoord_params[1]
	print overcoord_params[2]
	print overcoord_params[3]
	print overcoord_params[4]

print "fit_pover " + fit_pover
print "least squares parameters"	


for i in range(0,len(x)):
    print i,x[i]

if TEST_SUITE_RUN == "do":
    test_suite_params=open("test_suite_params.txt","w")		
    for i in range(0,len(x)):
        phrase = `i` + " " + `x[i]` + '\n'
        test_suite_params.write(phrase)
    test_suite_params.close()
