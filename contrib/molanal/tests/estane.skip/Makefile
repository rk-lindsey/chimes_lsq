estanetest : clean
	export OMP_NUM_THREADS=4 ; ../../molanal.new estane.gen > estane.out
	../../skiphead.pl estane.out > estane.tail
	../../findmolecules.pl estane.out > findmolecules.out
	diff estane.tail estane.save
	diff molanal.xyz molanal.save
	diff -r molecules molecules.save --exclude=.svn
	../../skiphead.pl findmolecules.out > findmolecules.tail
	diff findmolecules.tail findmolecules.save

clean:
	rm -rf molanal.xyz molecules *.out *.com *.tail *.log
