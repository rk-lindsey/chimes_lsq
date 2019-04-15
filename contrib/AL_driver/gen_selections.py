# Global (python) modules

import glob
import os
import sys
import argparse
import numpy as np
import random
import matplotlib 			# This and the next command prevents matplotlib from requiring an x-server (useful when screen is used)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from cycler import cycler
import math as m

# Local modules

import helpers


def populate_repo(my_ALC):

	if not os.path.isdir("../CENTRAL_REPO"):
		helpers.run_bash_cmnd("mkdir ../CENTRAL_REPO")

	# Create the list of selected species for the current ALC
	
	helpers.run_bash_cmnd("rm -f ../CENTRAL_REPO/ALC-" + `my_ALC` + ".all_selections.xyzlist")
	
	ifstream = open("all.xyzlist.dat",'r')
	all_xyz  = ifstream.readlines()
	ifstream .close()
	
	ifstream = open("all.selection.dat",'r')
	all_sel  = ifstream.readlines()
	ifstream .close()
	
	ofstream = open("../CENTRAL_REPO/ALC-" + `my_ALC` + ".all_selections.xyzlist", 'w')
	
	
	for j in xrange(len(all_sel)):
		for i in xrange(len(all_xyz)):

			if i == int(all_sel[j].rstrip()):
			
				line = all_xyz[i].split()

				line = line[0] + " " + line[1] + " " + line[2] + " " + "../ALC-" + `my_ALC` + "/" + line[3]

				ofstream.write(line + '\n')
				continue
	ofstream.close()
	
	
	# Recompile the central repo list of files

	helpers.run_bash_cmnd("rm -f ../CENTRAL_REPO/full_repo.xyzlist")
		
	helpers.cat_pattern("../CENTRAL_REPO//full_repo.xyzlist", ' '.join(glob.glob("../CENTRAL_REPO/*.all_selections.xyzlist")))
	
	
	

def GET_BIN(val,bins):

	for i in xrange(len(bins)-1):
	
		if (val >= bins[i]) and (val < bins[i+1]):
			return i
		elif val == bins[i+1]:
			return i
			
	print "PROBLEM: No bin was found for value ", val
	print bins[0]
	print bins[len(bins)-1]
	exit()
			
def GEN_NORM_HIST(idx_list, ener_list, central_repo_enerlist, minval, maxval, nbins):

	ener_vals = []

	for i in xrange(len(idx_list)):
		ener_vals.append(ener_list[idx_list[i]])
	
	ener_vals += list(central_repo_enerlist)

	prob, bins = np.histogram(ener_vals, bins=nbins, range=(minval,maxval), density=False)

	if(float(np.sum(prob)) == 0.0):
		print "Problem, found a zero sum:"
		print prob
		print 
		exit()
	
	prob = prob/float(np.sum(prob))

	return prob, bins	
		
	
def SET_CONDITIONS(nsweeps):
	
	arr = []

	for i in xrange(5):
		arr.append(m.exp(float(i)*2))
		
	for i in xrange(5):
		arr[i] = m.ceil((arr[i] - arr[0])/(arr[4] - arr[0])*(nsweeps-1))
		
	print "Stats will be printed for sweeps: ",arr
	return arr


def gen_subset(**kwargs): # time python gen_subset.py  all.energies_normed $SELECTIONS $SWEEPS 0 # Last 2 args: # to select, # sweeps, (optional:) E-cutoff


	# Notes:
	#
	# Expects to be run from the ALC-X directory
	# Expects "all.energies_normed" in the ALC-X directory 
	# 
	# 
	# 
	
	
	################################
	# 0. Set up an argument parser
	################################	
	
	default_keys   = [""]*6
	default_values = [""]*6        
	
	# Cluster specific controls
	
	default_keys[0 ] = "energies"	 ; default_values[0 ] = 'all.energies_normed' # List of energies to select from  
	default_keys[1 ] = "nsel"	 ; default_values[1 ] = '400'		      # Number of selections to make	 
	default_keys[2 ] = "nsweep"	 ; default_values[2 ] = '200000'	      # Number of MC sqeeps		 
	default_keys[3 ] = "nbins"	 ; default_values[3 ] = '20'		      # Number of histogram bins	 
	default_keys[4 ] = "ecut"	 ; default_values[4 ] = '1.0E10'	      # Maximum energy to consider	 
	default_keys[5 ] = "repo"	 ; default_values[5 ] = ''		      # Location of central repo energies

		
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	ENER    =       args["energies"]
	NSELECT =   int(args["nsel"    ])
	NSWEEPS =   int(args["nsweep"  ]) 
	NBINS   =   int(args["nbins"   ]) 
	ECUTOFF = float(args["ecut"    ])
	REPENER =       args["repo"    ]	
	
	print "Parsed arguments: "
	print "energies: ",ENER   
	print "nsel:     ",NSELECT
	print "ncyc:     ",NSWEEPS
	print "nbins:    ",NBINS
	print "ecut:     ",ECUTOFF
	print "repo:     ",REPENER	
	

	################################
	# 1. Read the energy file(s), Identify max/min values
	################################

	# Energies to select from:
	
	ENER    = np.loadtxt(ENER)

	if NSELECT > ENER.shape[0]:
		print "ERROR: NSELECT is larger than the number of available energies: "
		print "NSELECT:   ", NSELECT
		print "NENERGIES: ", ENER.shape[0]
		exit()
		
	# Find the max and min energy values
	
	MAX_VAL = ENER[0]
	MIN_VAL = ENER[0]
	
	MIN_FROM_MAIN = True # Assume max/min values are coming from ENER and not REPENER
	MAX_FROM_MAIN = True
	
	MAX_IDX = 0
	MIN_IDX = 0

	PAST_CUT_VAL = []
	PAST_CUT_IDX = []

	for i in xrange(len(ENER)):
	
		if abs(ENER[i]) >= ECUTOFF:
			PAST_CUT_VAL.append(ENER[i])
			PAST_CUT_IDX.append(i)	
		else:	

			if MAX_VAL < ENER[i]:
				MAX_VAL = ENER[i]
				MAX_IDX = i

			if MIN_VAL > ENER[i]:
				MIN_VAL = ENER[i]
				MIN_IDX = i
			
	
	# Central repository energies
	
	IGNORE_REPO_GTMAX = 0
	IGNORE_REPO_LTMIN = 0
	
	if REPENER != '':
	
		REPENER = np.loadtxt(REPENER)	
		
		REPO_POP_LIST = []
		
		for i in xrange(len(REPENER)):
			if abs(REPENER[i]) >= ECUTOFF:
				print "WARNING: Repo energy is outside ECUTOFF!"
				REPO_POP_LIST.append(i)
				
			if REPENER[i] > MAX_VAL:
				#print "WARNING: Found central repo energy larger than current iteration\'s ... ignoring."
				IGNORE_REPO_GTMAX += 1
				REPO_POP_LIST.append(i)
			elif REPENER[i] < MIN_VAL:
				#print "WARNING: Found central repo energy smaller than current iteration\'s ... ignoring."
				IGNORE_REPO_LTMIN += 1
				REPO_POP_LIST.append(i)
				
		if len(	REPO_POP_LIST ) >  0:
		
			TMP = list(REPENER)
			
			REPO_POP_LIST.sort(reverse=True)
			
			for i in xrange(len(REPO_POP_LIST)):
			
				TMP.pop(REPO_POP_LIST[i])
			
			REPENER = np.array(TMP)
			
		
		for i in xrange(len(REPENER)):

			if MAX_VAL < REPENER[i]:
				MAX_VAL = REPENER[i]
				MAX_IDX = i
				MAX_FROM_MAIN = False

			if MIN_VAL > REPENER[i]:
				MIN_VAL = REPENER[i]
				MIN_IDX = i
				MIN_FROM_MAIN = False
	else:
		REPENER = []
		
	#PAST_CUT_VAL.sort(reverse=True)	
	#PAST_CUT_IDX.sort(reverse=True)
				
	print "FOUND MIN/MAX VALS (from ENER?): ", MIN_VAL, MAX_VAL, MIN_FROM_MAIN, MAX_FROM_MAIN 
	print "FOUND MIN/MAX IDXS             : ", MIN_IDX, MAX_IDX	
	print "Ignored repo vals (<min)       : ", IGNORE_REPO_LTMIN
	print "Ignored repo vals (>max)       : ", IGNORE_REPO_GTMAX
				
	########################
	# Generate the probability distribution for the energy list
	# (a normalized numpy histogram; NBINS bins; NOT including REPENER
	# values except max/min)
	########################
	
	# Generate

	PROB, BINS = np.histogram(ENER, bins=NBINS, range=(MIN_VAL,MAX_VAL), density=False)
	
	if(float(np.sum(PROB)) == 0.0):
		print "Problem, found a zero sum:"
		print PROB
		print 
		exit()
		
	print "Total number of configurations: ",sum(PROB)
	print "Number of configurations in each bin: ",PROB
	
	NSWEEPS *= len(ENER)
		
	PROB = PROB/float(np.sum(PROB))
	


	# Plot Original random selection results ( and set color scheme)
	
	colors = plt.cm.OrRd(np.linspace(0, 1, 10)[::-1]) # See https://matplotlib.org/tutorials/colors/colormaps.html for options
	plt.rcParams['axes.prop_cycle'] = cycler('color', colors)	

	MID = (BINS[:-1] + BINS[1:]) / 2.0	# Plot (use bin midpoints)
	plt.plot(MID, PROB, marker='o', fillstyle='none', label='all current iteration configs')


	########################
	# Generate the initial sub-selection
	########################
	
	random.seed(1)
	
	# Step 0: Combine the new (ENER) and old (REPENER) energies into REPO. If
	#         min/max are from ENER, remove (pop off) the from REPO b/c they 
	#         will always be in the selected subset then we can always pick from 
	#         the first len(ENER)-NPOPPED configuations
	
	NENER = ENER.shape[0]
	
	#COMBINED_ENERS = np.concatenate((ENER, REPENER)) # ENER + REPENER # First len(ENER) entries are ALWAYS from ENER	
	
	REPO = list(xrange(NENER))

	NPOPPED = 0
	NCUTOFF = 0
	
	POP_LIST = list(PAST_CUT_IDX)

	NCUTOFF = len(PAST_CUT_IDX)
	
        if MAX_FROM_MAIN:
		POP_LIST.append(MAX_IDX)
                NPOPPED += 1

        if MIN_FROM_MAIN:
		POP_LIST.append(MIN_IDX)
                NPOPPED += 1
		
	POP_LIST.sort(reverse=True)
	
	for i in xrange(len(POP_LIST)):
		REPO.pop(REPO.index(POP_LIST[i]))
	
	print "Values above energy cutoff: ", NCUTOFF
		
	# Step 1: Select NSELECT-NPOPPED random energies only from SELE elemnts from ENER! (repetition not allowed)
	
	#print "Selecting random values on: ",0, NENER-NPOPPED
	#print "Selecting nvalues:          ",NSELECT-NPOPPED
	

	SELE = random.sample(REPO[0:(NENER-NPOPPED-NCUTOFF)],NSELECT-NPOPPED) 	# random sample returns *VALUES* of REPO to SELE... so SELE[i] = energy list index

	# Need to figure out the corresponding REPO index for the values stored in SELE, before popping them off

	for i in xrange(len(SELE)):
		REPO.pop(REPO.index(SELE[i]))
		
	# Determine the number of elements in SELE arising from ENER, and in REPO arising from ENER
	
	N_ENER_SELE = len(SELE)
	N_ENER_REPO = NENER - NPOPPED - N_ENER_SELE - NCUTOFF
	
	
	########################
	# Do a MC sweep to update selections
	########################
	#
	# Our scheme: accept if P_OLD - P_NEW > rand(0,1)
	# This will bias toward new configurations that were 
	# lower in probability, and should flatten our hist
	
	
	# Set up conditions for printing current histogram
	
	CONDITIONS = SET_CONDITIONS(NSWEEPS)

	# Generate the probability histogram for the MC run... histogram is over selected new energies (ENER) and all old energies (REPENER)

	HIST_SELE = SELE[:]

	if MIN_FROM_MAIN:
		HIST_SELE += [MIN_IDX]
	if MAX_FROM_MAIN:
		HIST_SELE += [MAX_IDX]
		
	REMAINING =  REPO[N_ENER_REPO:] # These are values from the central repository
		
	HIST_SELE += REMAINING

	SELE_PROB, SELE_BINS = GEN_NORM_HIST(HIST_SELE, ENER, REPENER, MIN_VAL, MAX_VAL, NBINS)

	# Do MC sweeps

	HIST_SELE = []
	SSQR_LIST = []
	SSQR      = "not yet calculated"

	if NSELECT == ENER.shape[0]:
		NSWEEPS = 1
		HIST_SELE = SELE[:]
		
		if MIN_FROM_MAIN:
			HIST_SELE += [MIN_IDX]
		if MAX_FROM_MAIN:
			HIST_SELE += [MAX_IDX]
			
		HIST_SELE += REMAINING

		SELE_PROB, SELE_BINS = GEN_NORM_HIST(HIST_SELE, ENER, REPENER, MIN_VAL, MAX_VAL)

	
	for sweep in xrange(NSWEEPS):

		if sweep in CONDITIONS: #if ((NSWEEPS/5>0) and (sweep+1)%(NSWEEPS/5)==0):
			print "Running MC sweep " + `sweep` + " of " + `NSWEEPS` + " ... SSQR: " +  `SSQR` 
		
		NONE_ACC = True

		for i in xrange(NSELECT-NPOPPED):
		
			# print "	Running MC step " + `i` + " of " + `NSELECT-2`

			if NSELECT < ENER.shape[0]:		

				# Select an existing and new energy (OLD and NEW = index)

				
				OLD = random.randint(0,N_ENER_SELE-1) # provides a SELE index		   # One from the current subset
				NEW = random.randint(0,N_ENER_REPO-1) # Provides a REPO index		   # One from the possible repo of configs


				if (OLD >len(ENER)) or (NEW >len(ENER)):
					print "++++"				
					print "ERROR: ", OLD, NEW
					
					print len(ENER)
					ENER.sort(reverse=True)
					print ENER[0]
					print "++++"
					
					exit()

				# Apply the acceptance criteria

				P_OLD = SELE_PROB[GET_BIN(ENER[SELE[OLD]],SELE_BINS)]
				P_NEW = SELE_PROB[GET_BIN(ENER[REPO[NEW]],SELE_BINS)]
						
				
				RAND  = random.random()
				CRIT  = 0.5*(P_NEW-P_OLD) + 1.0
				#CRIT  = (P_NEW-P_OLD) #  + 1.0	
				#CRIT = 1.0 + (P_OLD - P_NEW)*0.5
				
				if False:			
				
					print ""
					print " 	      old:  ", P_OLD
					print " 	      new:  ", P_NEW
					print " 	      rand: ", RAND
					print " 	      crit: ", CRIT
			
				if ( CRIT > RAND): # This is a high prob cfg... we want to bias against it
					continue
	
				NONE_ACC = False
					
				# We've accepted the move... update the selected and stored respositories and re-compute the probabilities

				SELE_VAL = SELE.pop(OLD)
				REPO_VAL = REPO.pop(NEW) # .pop returns REPO[NEW], and removes element [NEW] from REPO
			
				SELE.append(REPO_VAL)
				REPO.append(SELE_VAL)

			
			
			if any( x>len(ENER) for x in SELE):
				
				SELE.sort(reverse=True)
			
				print "Warning: Selected value too large: "
				print len(SELE)
				print SELE[0]
				print len(ENER)
				print ENER			
				exit()
			
			
			HIST_SELE = SELE[:]
			
			if MIN_FROM_MAIN:
				HIST_SELE += [MIN_IDX]
			if MAX_FROM_MAIN:
				HIST_SELE += [MAX_IDX]
				
			HIST_SELE += REMAINING

			SELE_PROB, SELE_BINS = GEN_NORM_HIST(HIST_SELE, ENER, REPENER, MIN_VAL, MAX_VAL, NBINS)
			

			if NSELECT == ENER.shape[0]:
				break
			
			
		# Compute sum of squared residuals (our "equilibration" criteria)	
			
		if NONE_ACC:
			if len(SSQR_LIST) == 0:
				SSQR_LIST.append(-1)
				
			SSQR = SSQR_LIST[len(SSQR_LIST)-1]
			SSQR_LIST.append(SSQR)
		else:
			TARGET = 1/float(NBINS) # 1/nbins
			SSQR   = 0.0
			
			for i in xrange(len(SELE_PROB)):
		
				#SSQR += (TARGET - ENER[HIST_SELE[i]])**2.0
				SSQR += (TARGET - SELE_PROB[i])**2.0
				
			try: 
				SSQR = m.sqrt( SSQR / len(HIST_SELE)) 
			except: 
				print SELE
				print HIST_SELE
		
			SSQR_LIST.append(SSQR)	
		
		
		# Plot current selection results 

		SELE_MID = (SELE_BINS[:-1] + SELE_BINS[1:]) / 2.0
		LABEL    = "sweep " + `sweep+1`

		if sweep in CONDITIONS: # if ((NSWEEPS/5>0) and (sweep+1)%(NSWEEPS/5)==0):
			plt.plot(SELE_MID, SELE_PROB, marker='x', label=LABEL)

	
	TMP_SELE = SELE[:]

	if MIN_FROM_MAIN:
		TMP_SELE += [MIN_IDX]
	if MAX_FROM_MAIN:
		TMP_SELE += [MAX_IDX]	
		
	# Save resulting selection

	np.savetxt("selection.dat",np.array(TMP_SELE).astype(int),fmt='%5d')

	#print "Probabilites for final selection: "
	#print SELE_PROB

	# Plot results


	#plt.legend(loc='upper left')
	plt.legend_ = None
	plt.savefig('energy_hist.pdf')
	plt.clf()
	plt.cla()
	plt.close()

	plt.plot(list(xrange(1,len(SSQR_LIST)+1)), SSQR_LIST)
	plt.xscale('log')
	plt.savefig('residuals.pdf')

	helpers.run_bash_cmnd("mv selection.dat all.selection.dat")










	
	
	
	
