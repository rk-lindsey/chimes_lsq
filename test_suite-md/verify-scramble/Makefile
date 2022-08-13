# Check that forces are not dependent on atom order.  
# The forces are compared for 2 xyz
# files.  One that is unscrambled (input.xyz) and one that is scrambled
# (scramble.xyz)

NP ?= 1

check_scramble: 
	$(RUN_JOB) ../../build/chimes_md scramble_md.in > scramble_4b_force.out	
	mv forceout.txt forceout.scramble.txt
	sort --general-numeric forceout.scramble.txt > scramble.sort.txt
	$(RUN_JOB) ../../build/chimes_md run_md.in > orig_4b_force.out	
	mv forceout.txt forceout.orig.txt
	sort --general-numeric forceout.orig.txt > orig.sort.txt
	perl ../../contrib/compare/compare.pl scramble.sort.txt orig.sort.txt

clean :
	rm -f force.txt *.out forceout*.txt *.sort.txt *.sort2.txt output.* traj.gen
