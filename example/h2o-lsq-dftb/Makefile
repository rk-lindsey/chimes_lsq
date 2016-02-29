dftb :
	cd dftb_compute_erep_forces ; dftb_fmatch_input.pl input_openmx_1.5gcc.xyzf O H > dftb.out ;  \
	mv input_dftb.xyzf ../lsq/input.xyzf
	cd lsq ; ../../../bin/splines_ls dftb_ls.in 63 > dftb_ls.out ; \
	../../../bin/lsq.py A.txt b.txt params.header > params.dftb.txt ; \
	head -n 5 params.dftb.txt ; \
	./get_new_skf_h2o.pl params.dftb.txt ; \
	./compute_erep.pl params.dftb.txt 63

clean :
	rm -f A.txt b.txt params.header params.dftb.txt force.txt


