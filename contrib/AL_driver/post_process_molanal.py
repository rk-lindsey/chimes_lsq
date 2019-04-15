import sys

import helpers

# Post-process (extract relevant speciation)


args_species = sys.argv[1:]
nspecies = len(args_species)

for i in xrange(nspecies):

	ofname   = "speciation.species-" + `i`
	ofstream = open(ofname,'w')
	
	ofstream.write("# " + args_species[i]    + '\n')
	ofstream.write("# <mol_frac> <lifetime>" + '\n')		
	
	ifname   = "traj.gen-find_molecs.out"
	ifstream = open(ifname,'r')
	ifile    = ifstream.readlines()

	for j in xrange(len(ifile)):
	
		if args_species[i] in ifile[j]:
			
			print "found"
		
			tmp = ifile[j].split()
			
			ofstream.write(tmp[len(tmp)-5] + " " + tmp[len(tmp)-2] + '\n')
			
			break
	ofstream.close()
	ifstream.close()

	
