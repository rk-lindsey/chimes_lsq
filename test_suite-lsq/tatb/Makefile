## These variables are defaults.
NP=36
NODES=1
RUN_JOB=srun -n $(NP) -N $(NODES)
PYTHON=python # /usr/tce/bin/python
COMPARE=perl ../../contrib/compare/compare.pl
DLARS_PATH=../../../contrib/dlars/src/

all: fm_setup.out params.dlars.txt

fm_setup.out:
	cp fm_setup.in current_output/
	cp SCAN_all.xyzf current_output/
	cd current_output ; $(RUN_JOB) ../../../build/chimes_lsq fm_setup.in > fm_setup.out
	-$(COMPARE) current_output/fm_setup.out correct_output/fm_setup.out
	-$(COMPARE) current_output/b-labeled.txt correct_output/b-labeled.txt

generate: all
	cp current_output/fm_setup.out correct_output/
	cp current_output/params.dlars.txt correct_output/
	cp current_output/b-labeled.txt correct_output/

params.dlars.txt:
	cd current_output ; $(PYTHON) ../../../src/chimes_lsq.py --algorithm=dlars --dlasso_dlars_path=${DLARS_PATH} --alpha=1.0e-04 --cores=$(NP) --nodes=$(NODES) --split_files=true --normalize=false > params.dlars.txt
	-$(COMPARE) current_output/params.dlars.txt correct_output/params.dlars.txt

clean:
	rm -f current_output/* 
