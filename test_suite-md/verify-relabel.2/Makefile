## Test to see whether 1 element and 2 element systems are identical
## if the 1 element system is relabeled as a 2 element system.
##

all:
	../chimes_md run_md.oh.in > run_md.oh.out
	mv forceout.txt forceout.oh.txt
	../chimes_md run_md.o.in > run_md.o.out
	mv forceout.txt forceout.o.txt
	diff forceout.oh.txt forceout.o.txt

clean:
	rm -f forceout*.txt *.out *.gen *.bak traj_bad.xyz

