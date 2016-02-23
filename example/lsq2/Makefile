all :
	cat omx1.xyzf omx2.xyzf > input.xyzf
	../../bin/splines_ls splines_ls.in 60
	../../bin/lsq.py A.txt b.txt params.header > params.spline2.txt
	-diff params.spline2.txt output/params.spline2.txt

clean :
	rm -f A.txt b.txt params.header params.spline2.txt force.txt


