all : cheby_fit cheby_force cheby_md

cheby_fit :
	../../bin/splines_ls cheby_ls_3b.in 49 > cheby_ls_3b.out
	../../bin/lsq.py A.txt b.txt params.header > params.3b.txt
	../../bin/cheby_print.pl < params.3b.txt
	-diff chebypot.0.out output/chebypot.0.out
	-diff chebypot.1.out output/chebypot.1.out
	-diff chebypot.2.out output/chebypot.2.out
	-diff cheby.0.out output/cheby.0.out
	-diff cheby.1.out output/cheby.1.out
	-diff cheby.2.out output/cheby.2.out
	touch cheby_fit

cheby_md : cheby_fit
	../../bin/splines_md cheby_md_3b.in > cheby_md_3b.out
	-diff cheby_md_3b.out output/cheby_md_3b.out

cheby_force : cheby_fit
	../../bin/splines_md cheby_3b_force.in > cheby_3b_force.out
	-diff cheby_3b_force.out output/cheby_3b_force.out

# Special check of code.
check_sorting :
	echo 'splines_md must be compiled with TESTING on in Cheby3B.C'
	echo 'Need a lot of disk space'
	../../bin/splines_md cheby_3b_force.in > cheby_3b_force.out
	../../bin/check_sorting.pl < cheby_3b_force.out > check_sorting.out

# Check of scrambling atom order.  Should be same results as cheby_md
check_scramble:
	mv input.xyz input.save
	cp scramble.xyz input.xyz
	-../../bin/splines_md cheby_md_3b.in > scramble_md_3b.out	
	cp input.save input.xyz

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt *.out cheby_fit


