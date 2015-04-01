CC=gcc
CFLAGS=-g -Wall

all : nrfit

nrfit : nrfit.o nrutil.o dfpmin.o amebsa.o powell.o linmin.o ran1.o mnbrak.o f1dim.o brent.o amoeba.o amotry.o
	$(CC) nrfit.o nrutil.o dfpmin.o amebsa.o powell.o linmin.o ran1.o \
brent.o mnbrak.o f1dim.o amoeba.o amotry.o -lm -o nrfit
	cp nrfit ../../bin

.c.o :
	$(CC) -c  $(CFLAGS) $(DEFINES) $(INCLUDES) $*.c 

clean : 
	rm -f *.o *~ *.out

realclean :
	rm -f nrfit

## All available optimization options are tested.  In some cases
## Powell's method follows a global optimization method.
test : nrfit
	rm -f *coeff.out
	./nrfit -p < fit.in >& fit.out
	./checkparam.pl pcoeff.out
	./nrfit -a asa.in < fit.in >& fit.out
	./checkparam.pl asacoeff.out
	./nrfit -s sa.in < fit.in >& fit.out
	./checkparam.pl scoeff.out
	./nrfit -b < fit.in >& fit.out
	./checkparam.pl bcoeff.out
	rm -f pcoeff.out
	./nrfit -r record.in -p < fit.in >& fit.out
	./checkparam.pl pcoeff.out
	rm -f pcoeff.out
	./nrfit -g gda.in -p < fit.in >& fit.out
	./checkparam.pl pcoeff.out
	./nrfit -d dfp.in < fit.in >& fit.out
	./checkparam.pl dcoeff.out

