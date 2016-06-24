all : cheby_fit cheby_force cheby_md check_scramble

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

# Check of scrambling atom order.  The forces are compared for 2 xyz
# files.  One that is unscrambled (input.xyz) and one that is scrambled
# (scramble.xyz)
# The "unscrambled" forces are places in forceout.scramble.order.txt for
# comparison to forceout.orig.txt
check_scramble:
	../../bin/splines_md scramble_3b_force.in > scramble_3b_force.out	
	mv forceout.txt forceout.scramble.txt
	sort forceout.scramble.txt > scramble.sort.txt
	../../bin/splines_md orig_3b_force.in > orig_3b_force.out	
	mv forceout.txt forceout.orig.txt
	../../bin/match_xyz.pl input.xyz scramble.xyz forceout.scramble.txt > forceout.scramble.order.txt
	sort forceout.orig.txt > orig.sort.txt
	diff scramble.sort.txt orig.sort.txt

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt *.out cheby_fit forceout*.txt *.sort.txt


