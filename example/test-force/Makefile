# This tests the consistency of forces produced by lsq.py with those
# calculated by splines_md.C.
#
all : cheby spline

cheby :
	../../bin/splines_ls cheby_ls.in 49 > cheby_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.cheby.txt
	../../bin/cheby_print.pl < params.cheby.txt
	../../bin/splines_md cheby_md.in > cheby_md.out
	-diff cheby_md.out output/cheby_md.out

spline :
	../../bin/splines_ls splines_ls.in 49 > splines_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.spline.txt
	../../bin/splines_md splines_md.in > splines_md.out
	-diff splines_md.out output/splines_md.out

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt *.out


