all : cheby_fit cheby_md

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

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt *.out cheby_fit


