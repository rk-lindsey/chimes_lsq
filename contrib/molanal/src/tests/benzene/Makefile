benzenetest : clean
	../../molanal.new Benzene.md.gen > benzene.out
	../../skiphead.pl benzene.out > benzene.tail
	../../findmolecules.pl benzene.out > findmolecules.out
	../../skiphead.pl findmolecules.out > findmolecules.tail
	-diff benzene.tail benzene.save
	-diff findmolecules.tail findmolecules.save

clean:
	rm -rf molanal.xyz molecules *.out *.com *.tail *.log

