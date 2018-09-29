NP=36
RUN_JOB=srun -n $(NP)
PYTHON=python
COMPARE=perl ../../contrib/compare/compare.pl

all: A.txt svd lassolars ridge

A.txt:
	$(RUN_JOB) ../chimes_lsq < fm_setup.in > fm_setup.out

%.cmp : %.txt
	$(COMPARE) current_output/$*.txt correct_output/$*.txt >& $*.cmp

svd: params.svd.1.0e-04.cmp params.svd.1.0e-05.cmp

lassolars: params.lassolars.1.0e-04.cmp 

ridge: params.ridge.1.0e-03.cmp

ridgecv: params.ridgecv.cmp

params.ridgecv.txt : A.txt
	$(PYTHON) ../../src/lsq2.py --algorithm=ridgecv > params.ridgecv.txt	
	mv params.ridgecv.txt current_output/

params.ridge.1.0e-03.txt: A.txt
	$(PYTHON) ../../src/lsq2.py --algorithm=ridge --alpha=1.0e-03 > params.ridge.1.0e-03.txt
	mv params.ridge.1.0e-03.txt current_output/

params.lassolars.1.0e-04.txt: A.txt
	$(PYTHON) ../../src/lsq2.py --algorithm=lassolars --alpha=1.0e-04 > params.lassolars.1.0e-04.txt
	mv params.lassolars.1.0e-04.txt current_output/

params.svd.1.0e-04.txt: A.txt
	$(PYTHON) ../../src/lsq2.py --algorithm=svd --eps=1.0e-04 > params.svd.1.0e-04.txt
	mv params.svd.1.0e-04.txt current_output/

params.svd.1.0e-05.txt: A.txt
	$(PYTHON) ../../src/lsq2.py --algorithm=svd --eps=1.0e-05 > params.svd.1.0e-05.txt
	mv params.svd.1.0e-05.txt current_output/

clean:
	rm *.txt *.cmp
