# Check that forces are consistent with the energy expression by
# evaluating forces numerically.
#

NP ?= 16

check_forces:
	$(RUN_JOB) ../../build/chimes_md run_md.in > run_md.out
	grep -i 'fail' run_md.out ; if [ $$? -eq 0 ] ; then echo 'run failed' ; exit 1 ; fi

clean:
	rm -f force.txt *.out forceout*.txt output.* traj.gen traj_bad.xyz *~

