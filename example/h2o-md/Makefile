all : spline cheby

spline : 
	../../bin/splines_md splines_md.in > splines_md.out
	-diff splines_md.out output/splines_md.out

cheby : 
	../../bin/splines_md cheby_md.in > cheby_md.out
	-diff cheby_md.out output/cheby_md.out
clean :
	rm -f splines_md.out cheby_md.out traj.gen output.xyz

realclean : clean
