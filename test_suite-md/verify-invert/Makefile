# Check that forces are not dependent on translation of all atoms.
# The forces are compared for 2 xyz
# files.  

NP ?= 1

check_invert: 
	./invert.pl < input.xyz > invert.xyz
	$(RUN_JOB) ../../build/chimes_md invert_md.in > invert_3b_force.out	
	mv forceout.txt forceout.invert.txt
	$(RUN_JOB) ../../build/chimes_md run_md.in > orig_3b_force.out	
	mv forceout.txt forceout.orig.txt
	./compare-invert.pl forceout.orig.txt forceout.invert.txt

clean :
	rm -f force.txt *.out forceout*.txt *.sort.txt *.sort2.txt output.* traj.gen *bad.xyz
