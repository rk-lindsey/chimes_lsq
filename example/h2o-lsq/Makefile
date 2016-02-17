all : cheby spline

cheby :
	../../bin/splines_ls cheby_ls.in 49 > cheby_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.cheby.txt
	../../bin/cheby_print.pl < params.cheby.txt
	diff chebypot.0.out output/chebypot.0.out
	diff chebypot.1.out output/chebypot.1.out
	diff chebypot.2.out output/chebypot.2.out
	diff cheby.0.out output/cheby.0.out
	diff cheby.1.out output/cheby.1.out
	diff cheby.2.out output/cheby.2.out


spline :
	../../bin/splines_ls splines_ls.in 49 > splines_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.spline.txt
	../../bin/spline_vals.pl < params.spline.txt
	diff splinepot.0.out output/splinepot.0.out
	diff splinepot.1.out output/splinepot.1.out
	diff splinepot.2.out output/splinepot.2.out
	diff spline.0.out output/spline.0.out
	diff spline.1.out output/spline.1.out
	diff spline.2.out output/spline.2.out

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt


