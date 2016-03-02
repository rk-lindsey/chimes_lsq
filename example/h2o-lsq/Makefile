all : cheby cheby2 spline spline2

# tests spline2 and cheby2 use separate fitting parameters for each interaction.

cheby :
	../../bin/splines_ls cheby_ls.in 49 > cheby_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.cheby.txt
	../../bin/cheby_print.pl < params.cheby.txt
	-diff chebypot.0.out output/chebypot.0.out
	-diff chebypot.1.out output/chebypot.1.out
	-diff chebypot.2.out output/chebypot.2.out
	-diff cheby.0.out output/cheby.0.out
	-diff cheby.1.out output/cheby.1.out
	-diff cheby.2.out output/cheby.2.out


cheby2 :
	../../bin/splines_ls cheby_ls2.in 49 > cheby_ls2.out
	../../bin/lsq.py A.txt b.txt params.header > params.cheby2.txt
	../../bin/cheby_print.pl < params.cheby2.txt
	mv chebypot.0.out chebypot2.0.out
	mv chebypot.1.out chebypot2.1.out
	mv chebypot.2.out chebypot2.2.out
	mv cheby.0.out cheby2.0.out
	mv cheby.1.out cheby2.1.out
	mv cheby.2.out cheby2.2.out
	-diff chebypot2.0.out output/chebypot2.0.out
	-diff chebypot2.1.out output/chebypot2.1.out
	-diff chebypot2.2.out output/chebypot2.2.out
	-diff cheby2.0.out output/cheby2.0.out
	-diff cheby2.1.out output/cheby2.1.out
	-diff cheby2.2.out output/cheby2.2.out

spline :
	../../bin/splines_ls splines_ls.in 49 > splines_ls.out
	../../bin/lsq.py A.txt b.txt params.header > params.spline.txt
	../../bin/spline_vals.pl < params.spline.txt
	-diff splinepot.0.out output/splinepot.0.out
	-diff splinepot.1.out output/splinepot.1.out
	-diff splinepot.2.out output/splinepot.2.out
	-diff spline.0.out output/spline.0.out
	-diff spline.1.out output/spline.1.out
	-diff spline.2.out output/spline.2.out

spline2 :
	../../bin/splines_ls splines_ls2.in 49 > splines_ls2.out
	../../bin/lsq.py A.txt b.txt params.header > params.spline2.txt
	../../bin/spline_vals.pl < params.spline2.txt
	-mv splinepot.0.out splinepot2.0.out
	-mv splinepot.1.out splinepot2.1.out
	-mv splinepot.2.out splinepot2.2.out
	-mv spline.0.out spline2.0.out
	-mv spline.1.out spline2.1.out
	-mv spline.2.out spline2.2.out
	-diff splinepot2.0.out output/splinepot2.0.out
	-diff splinepot2.1.out output/splinepot2.1.out
	-diff splinepot2.2.out output/splinepot2.2.out
	-diff spline2.0.out output/spline2.0.out
	-diff spline2.1.out output/spline2.1.out
	-diff spline2.2.out output/spline2.2.out

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt


