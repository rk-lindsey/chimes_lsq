import sys
import math as m
import matplotlib.pyplot as plt

DATA = open(sys.argv[1],"r").readlines()

MIN = -150
MAX =  150
BINS = 150

HIST = [0]*BINS

WIDTH = (MAX-MIN)/BINS
SUM   = 0

for i in xrange(len(DATA)):
	DATA[i] = float(DATA[i])

for i in xrange(len(DATA)):

	BIN = int(   (m.floor((DATA[i]-MIN)/WIDTH)) )

	if BIN>=0 and BIN<len(HIST):
		HIST[BIN] += 1.0
	SUM       += 1.0

OFSTREAM = open("generated_histogram.dat",'w')

X = [0]*BINS
Y = [0]*BINS
	
for i in xrange(len(HIST)):
	X[i] = i*WIDTH+0.5*WIDTH+MIN
	Y[i] = HIST[i]/SUM/WIDTH
	OFSTREAM.write(`X[i]` + " " + `Y[i]` + "\n")
	
OFSTREAM.close()

#plt.plot(X, Y)
#plt.show()
