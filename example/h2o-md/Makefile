all : spline

spline : 
	../../bin/splines_md splines_md.in > splines_md.out
	-diff splines_md.out output/splines_md.out
clean :
	rm -f splines_md.out traj.gen output.xyz

realclean : clean
