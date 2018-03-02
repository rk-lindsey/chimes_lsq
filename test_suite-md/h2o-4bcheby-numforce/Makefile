# Check that forces are consistent with the energy expression by
# evaluating forces numerically.
#
check_forces:
	srun -N 1 -n 16 ../../src/chimes_md < run_md.in > run_md.out
	grep -i 'fail' run_md.out ; if [ $$? -eq 0 ] ; then echo 'run failed' ; exit 1 ; fi

clean:
	rm -f force.txt *.out forceout*.txt output.* traj.gen traj_bad.xyz *~

