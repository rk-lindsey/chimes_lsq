all : cheby_fit cheby_force cheby_md check_scramble cheby_fit2 cheby_md2 cheby_force2

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
	-grep -v Date params.3b.txt > params.3b.nodate.txt
	-diff params.3b.nodate.txt output/params.3b.nodate.txt
	touch cheby_fit

cheby_fit2 :
	../../bin/splines_ls cheby_ls_3b2.in 49 > cheby_ls_3b2.out
	../../bin/lsq.py A.txt b.txt params.header > params.3b2.txt
	../../bin/cheby_print.pl < params.3b2.txt
	-diff chebypot.0.out output/chebypot2.0.out
	-diff chebypot.1.out output/chebypot2.1.out
	-diff chebypot.2.out output/chebypot2.2.out
	-diff cheby.0.out output/cheby2.0.out
	-diff cheby.1.out output/cheby2.1.out
	-diff cheby.2.out output/cheby2.2.out
	-grep -v Date params.3b2.txt > params.3b2.nodate.txt
	-diff params.3b2.nodate.txt output/params.3b2.nodate.txt
	touch cheby_fit2

cheby_md : cheby_fit
	../../bin/splines_md cheby_md_3b.in > cheby_md_3b.out
	-diff cheby_md_3b.out output/cheby_md_3b.out

cheby_md2 : cheby_fit2
	../../bin/splines_md cheby_md_3b2.in > cheby_md_3b2.out
	-diff cheby_md_3b2.out output/cheby_md_3b2.out

cheby_force : cheby_fit
	../../bin/splines_md cheby_3b_force.in > cheby_3b_force.out
	-diff cheby_3b_force.out output/cheby_3b_force.out

cheby_force2 : cheby_fit2
	../../bin/splines_md cheby_3b_force2.in > cheby_3b_force2.out
	-diff cheby_3b_force2.out output/cheby_3b_force2.out

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

check_scramble2:
	../../bin/splines_md scramble_3b_force2.in > scramble_3b_force2.out	
	mv forceout.txt forceout.scramble2.txt
	sort forceout.scramble2.txt > scramble.sort2.txt
	../../bin/splines_md orig_3b_force2.in > orig_3b_force2.out	
	mv forceout.txt forceout.orig2.txt
	../../bin/match_xyz.pl input.xyz scramble.xyz forceout.scramble2.txt > forceout.scramble2.order.txt
	sort forceout.orig2.txt > orig.sort2.txt
	diff scramble.sort2.txt orig.sort2.txt

clean :
	rm -f A.txt b.txt params.header params.spline.txt force.txt *.out cheby_fit cheby_fit2\
 forceout*.txt *.sort.txt *.sort2.txt *.nodate.txt traj.gen params.*.txt
