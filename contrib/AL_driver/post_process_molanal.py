# Global (python) modules

import sys
import os

# Localmodules

import helpers

def post_process_molanal(args):

	""" 

	Extract relevant speciation from an MD run.

	Usage: python post_process_molanal.py "C1 O1 1(O-C)" "C1 O2 2(O-C)" "C2 O2 1(C-C) 2(O-C)" "C3 O2 2(C-C) 2(O-C)"

	Notes: Expects to find traj.gen-find_molecs.out in the working directory.

	       	
	"""
		
	# Post-process (extract relevant speciation)


	args_species = args
	nspecies     = len(args_species)

	for i in xrange(nspecies):

		ofname   = "speciation.species-" + `i`
		ofstream = open(ofname,'w')
		
		ofstream.write("# " + args_species[i]    + '\n')
		ofstream.write("# <mol_frac> <lifetime>" + '\n')	
		
		if not os.path.isfile("traj.gen-find_molecs.out"):	
		
			print "ERROR: Cannot find file traj.gen-find_molecs.out"
			os._exit(1)
		
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

	
if __name__=='__main__':

	""" 
	
	Allows commandline calls to post_process_molanal().
	
	      	
	"""	
	
	post_process_molanal(sys.argv[1:])
	
	

	
