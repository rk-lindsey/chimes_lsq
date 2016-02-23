all:
	../../bin/gen_to_omx_inputs.pl traj.gen openmx.header openmx.footer msub.header h2o
	-diff h2o.001.dat output/h2o.001.dat

clean:
	rm -f *.dat h2o.msub

realclean : clean
	rm -f *~

