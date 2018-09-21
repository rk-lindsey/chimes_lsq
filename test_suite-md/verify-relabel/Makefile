## Test to see whether 1 element and 2 element systems are identical
## if the 1 element system is relabeled as a 2 element system.
##
all: 2body 3body 4body

2body:
	../chimes_md < run_md.ch.2b.in > run_md.ch.2b.out
	mv forceout.txt forceout.ch.2b.txt
	../chimes_md < run_md.c.2b.in > run_md.c.2b.out
	mv forceout.txt forceout.c.2b.txt
	diff forceout.ch.2b.txt forceout.c.2b.txt


3body:
	../chimes_md < run_md.ch.3b.in > run_md.ch.3b.out
	mv forceout.txt forceout.ch.3b.txt
	../chimes_md < run_md.c.3b.in > run_md.c.3b.out
	mv forceout.txt forceout.c.3b.txt
	diff forceout.ch.3b.txt forceout.c.3b.txt

4body:
	../chimes_md < run_md.ch.4b.in > run_md.ch.4b.out
	mv forceout.txt forceout.ch.4b.txt
	../chimes_md < run_md.c.4b.in > run_md.c.4b.out
	mv forceout.txt forceout.c.4b.txt
	diff forceout.ch.4b.txt forceout.c.4b.txt

clean:
	rm -f forceout*.txt *.out *.gen *.bak traj_bad.xyz

