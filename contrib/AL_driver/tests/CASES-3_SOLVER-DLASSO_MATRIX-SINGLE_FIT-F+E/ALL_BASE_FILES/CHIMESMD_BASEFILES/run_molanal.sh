base="traj"

# Molanal

/g/g17/rlindsey/CURR_TRACKED-GEN/contrib/molanlal/molanal.new ${base}.gen | tee ${base}.gen-molanal.out
/g/g17/rlindsey/CURR_TRACKED-GEN/contrib/molanlal/findmolecules.pl ${base}.gen-molanal.out | tee ${base}.gen-find_molecs.out

nedit *-find_molecs.out
