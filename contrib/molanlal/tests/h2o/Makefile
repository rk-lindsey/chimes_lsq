h2otest : clean
	../../molanal.new ice_mix_test.gen > ice_mix_test.out
	../../skiphead.pl ice_mix_test.out > ice_mix_test.tail
	../../findmolecules.pl ice_mix_test.out > ice_mix_test.rxn
	../../skiphead.pl ice_mix_test.rxn > ice_mix_test2.rxn
	-diff ice_mix_test.tail ice_mix_test.save
	-diff molanal.xyz molanal.save
	-diff ice_mix_test2.rxn ice_mix_test.rxn.save
	-diff -r molecules molecules.save --exclude=.svn

clean :
	rm -rf molecules *.out *.xyz *.rxn *.com *.tail *.log
