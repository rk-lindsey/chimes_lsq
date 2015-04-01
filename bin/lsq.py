#! /usr/bin/python

import sys
import numpy
from numpy import *
##sys.path.append('/g/g90/koziol3/codes/mol/')
##import matrix
from numpy.linalg import lstsq
from datetime import *

if(len(sys.argv) != 4):
    print "./lsq.py A.txt b.txt params.header"
    sys.exit()


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
for i in range(0,len(Dmat)):
    for j in range(0,len(Dmat[i])):
        Dmat[i][j]=0.0

eps=1.0e-5
for i in range(0,len(D)):
    if(abs(D[i])>eps):
        Dmat[i][i]=1.0/D[i]


x=dot(transpose(VT),Dmat)
x=dot(x,dot(transpose(U),b))

y=dot(A,x)
Z=0.0

# Put calculated forces in force.txt
yfile = open("force.txt", "w")
for a in range(0,len(b)):
    Z=Z+(y[a]-b[a])**2.0
    yfile.write("%13.6e\n"% y[a]) 

print "# RMS force error = " , sqrt(Z/float(nlines))

print "# max variable = ",  max(x)

for i in range(0, len(hf)):
    sys.stdout.write(hf[i])

for i in range(0,len(x)):
    print i,x[i]

