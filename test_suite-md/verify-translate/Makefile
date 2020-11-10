# Check that forces are not dependent on translation of all atoms.
# The forces are compared for 2 xyz
# files.  

NP ?= 1

check_translate: 
	./translate.pl < input.xyz > translate.xyz
	srun -n $(NP)	../chimes_md translate_md.in > translate_4b_force.out	
	mv forceout.txt forceout.translate.txt
	srun -n $(NP)	../chimes_md run_md.in > orig_4b_force.out	
	mv forceout.txt forceout.orig.txt
	perl ../../contrib/compare/compare.pl forceout.orig.txt forceout.translate.txt

clean :
	rm -f force.txt *.out forceout*.txt *.sort.txt *.sort2.txt output.* traj.gen
