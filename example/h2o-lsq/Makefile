all : cheby spline

cheby :
	../../bin/splines_ls cheby_ls.in 49 > cheby_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.cheby.txt
	../../bin/cheby_print.pl < params.cheby.txt

spline :
	../../bin/splines_ls splines_ls.in 49 > splines_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.spline.txt
	../../bin/spline_vals.pl < params.spline.txt

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt


