all :
	../../bin/splines_ls splines_ls.in 50
	../../bin/lsq.py A.txt b.txt params.header > params.cheby.txt
	../../bin/cheby_print.pl < params.cheby.txt

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt


