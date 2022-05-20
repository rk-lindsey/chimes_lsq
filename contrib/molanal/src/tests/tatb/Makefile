tatbtest : clean
	export OMP_NUM_THREADS=4 ; ../../molanal.new tatb.gen > molanal.out
	../../skiphead.pl molanal.out > molanal.tail
	../../findmolecules.pl molanal.out > findmolecules.out
	../../skiphead.pl findmolecules.out > findmolecules.tail
	diff molanal.tail molanal.save
	diff molanal.xyz molanal.xyz.save
	diff -r molecules molecules.save --exclude=.svn
	diff findmolecules.tail findmolecules.save

clean :
	rm -rf molecules *.out *.xyz *.com *.tail *.log


