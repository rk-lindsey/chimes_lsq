# Global (python) modules

import glob # Warning: glob is unserted... set my_list = sorted(glob.glob(<str>)) if sorting needed
import helpers
import sys
import os

# See definition of main at the bottom of this file

def get_repo_energies(*argv, **kwargs):

	""" 
	
	Gets ChIMES energies for clusters in the CFG_REPO directory.
	
	Usage: get_repo_energies(<arguements>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Runs out of the ALC folder.
	       Expects a "params.txt" file to be there.
	       Calls .sh scripts in the utilities folder.	
	      	
	"""

	################################
	# 0. Set up an argument parser
	################################	
	
	default_keys   = [""]*12
	default_values = [""]*12	
	
	# Cluster specific controls
	
	default_keys[0 ] = "calc_central"      ; default_values[0 ] = False					# Recalculate central repository energies?
	default_keys[1 ] = "base_runfile"      ; default_values[1 ] = ""	        			# Directory containing run_md.base, etc.
	default_keys[2 ] = "driver_dir"	       ; default_values[2 ] = ""	        			# Directory containing these driver source files
	
	# Overall job controls
	
	default_keys[3 ] = "job_name"	       ; default_values[3 ] = 	 "dmb_enr"	        		# Name for ChIMES lsq job
	default_keys[4 ] = "job_nodes"         ; default_values[4 ] = 	 "4"			  	        # Number of nodes for ChIMES lsq job
	default_keys[5 ] = "job_ppn"	       ; default_values[5 ] = 	 "36"  		  	       	        # Number of processors per node for ChIMES lsq job
	default_keys[6 ] = "job_walltime"      ; default_values[6 ] = 	 "1"			  	        # Walltime in hours for ChIMES lsq job
	default_keys[7 ] = "job_queue"         ; default_values[7 ] = 	 "pdebug"		  	        # Queue for ChIMES lsq job
	default_keys[8 ] = "job_account"       ; default_values[8 ] = 	 "pbronze"		  	        # Account for ChIMES lsq job
	default_keys[9 ] = "job_executable"    ; default_values[9 ] = 	 ""			  	        # Full path to executable for ChIMES lsq job
	default_keys[10] = "job_system"        ; default_values[10] = 	 "slurm"	       	       	        # slurm or torque 
	default_keys[11] = "job_email"         ; default_values[11] =     True					# Send slurm emails?
		 
		
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)

	################################
	# 1. Prepare the energy jobs
	################################

	# Grab the submit scripts from the utitilities folder

	helpers.run_bash_cmnd("cp " + args["driver_dir"] + "/utilities/new-get_dumb_ener.sh	   .")
	helpers.run_bash_cmnd("cp " + args["driver_dir"] + "/utilities/new-get_dumb_ener_subjob.sh .")

	# Grab the run_md.base file, prepare for this job
	
	ifstream = open(args["base_runfile"],'r')
	runfile  = ifstream.readlines()

	ofstream = open("tmp",'w')	
	
	found1 = False
		
	for i in xrange(len(runfile)):
	
		if found1:
			ofstream.write('\t' + "false" + '\n')
			found1 = False		
		else:
	
			ofstream.write(runfile[i])
	
			if "PRNTBAD" in runfile[i]:
				found1 = True
	ofstream.close()
	helpers.run_bash_cmnd("mv tmp run_md.base")	
	
	
	################################
	# 2. Submit the local energy job
	################################	
	
	curr_dir = helpers.run_bash_cmnd("pwd").rstrip()

	run_cmnd  = "msub -l nodes=" + args["job_nodes"] + ":ppn=" + args["job_ppn"] + " -l walltime=" + args["job_walltime"] + ":00:00 -N dmb_enr " + " -q " +  args["job_queue"] + " " + " -A " +  args["job_account"] + " "
	if args["job_email"]:
		run_cmnd += " -m abe "
	run_cmnd += " new-get_dumb_ener.sh GEN_FF/params.txt.reduced "  + args["job_executable"] + " " + curr_dir + "/run_md.base " + args["driver_dir"]
	
	# Submit and monitor the jobs

	active_jobs = []
	active_jobs.append(int(helpers.run_bash_cmnd(run_cmnd).split()[0]))

	
	################################
	# 3. Submit the central energy job
	################################
	
	if args["calc_central"]:
	
		print "Launching central repo energy calculations as well"
	
		# Grab the submit scripts from the utitilities folder

		currdir = helpers.run_bash_cmnd("pwd").rstrip()
		
		helpers.run_bash_cmnd("cp " + args["driver_dir"]   + "/utilities/new-get_dumb_ener.sh	       ../CENTRAL_REPO")
		helpers.run_bash_cmnd("cp " + args["driver_dir"]   + "/utilities/new-get_dumb_ener_subjob.sh   ../CENTRAL_REPO")
		helpers.run_bash_cmnd("cp " + currdir + "/GEN_FF/params.txt.reduced			       ../CENTRAL_REPO")
		helpers.run_bash_cmnd("cp " + args["base_runfile"] + "  				       ../CENTRAL_REPO")   
		
		os.chdir("../CENTRAL_REPO")	
			
		
		# Update the base run_md file
		
		ifstream = open(args["base_runfile"],'r')
		runfile  = ifstream.readlines()
	
		ofstream = open("tmp",'w')			
		
		found1 = False

		for i in xrange(len(runfile)):
	
			if found1:
				ofstream.write('\t' + "params.txt.reduced" + '\n')
				found1 = False		
			else:
				ofstream.write(runfile[i])

				if "PRMFILE" in runfile[i]:
					found1 = True
				
		ofstream.close()
		helpers.run_bash_cmnd("mv tmp " + args["base_runfile"])		
		
		
		# Set up the job	
		
		helpers.run_bash_cmnd("cp full_repo.xyzlist xyzlist.dat")
	
		run_cmnd  = "msub -l nodes=" + args["job_nodes"] + ":ppn=" + args["job_ppn"] + " -l walltime=" + args["job_walltime"] + ":00:00 -N dmb_enr " + " -q " +  args["job_queue"] + " " + " -A " +  args["job_account"] + " "
		run_cmnd += " new-get_dumb_ener.sh params.txt.reduced "  + args["job_executable"] + " " + args["base_runfile"] + " " + args["driver_dir"]
		run_cmnd += " REPO"
	
		active_jobs.append(int(helpers.run_bash_cmnd(run_cmnd).split()[0]))

		################################
		# 4. Wait for repo job to end, post-process
		################################
	
		helpers.wait_for_job(active_jobs[1], job_system = args["job_system"], verbose = True, job_name = "get_repo_energies-central")	
	
		helpers.run_bash_cmnd("mv all.xyzlist.dat     full_repo.xyzlist.dat")
		helpers.run_bash_cmnd("mv all.energies        full_repo.energies")
		helpers.run_bash_cmnd("mv all.energies_normed full_repo.energies_normed")
	
		os.chdir(currdir)
		
		print "Central repo jobs finished"
	
	return active_jobs
	


def list_clusters(CFG_REPO, *argv): # This is where we need to start caring about the number of atoms types and thier order!

	""" 
	
	Generates a list of wrapped clusters in the repo directory.
	
	Usage: get_repo_energies(<arguements>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Assumes repo is located at current_ALC/<repo>.
	       Does NOT include "species" that contain fewer than 2 atoms.
	      	
	"""
	
	atm_types = argv[0] # This is a pointer!

	helpers.run_bash_cmnd("rm -f xyzlist.dat ts_xyzlist.dat")
	
	################################
	# X. Get a list of all tight species
	################################
	
	#ofstream = open("xyzlist.dat",'w')
	ofstream = open("tmp",'w')	# Sorting is required for reproducible runs
	
	repo_files = sorted(glob.glob(CFG_REPO + "/tight*.wrap*xyz"))
	
	for i in xrange(len(repo_files)):

		no_atoms = int(helpers.head(repo_files[i],1)[0].rstrip())
		no_types = [0]*len(atm_types)

		with open(repo_files[i], "r") as ifstream:
			
			for line in ifstream:
			
				for j in xrange(len(atm_types)):
				
					if atm_types[j] in line:
					
						no_types[j] += 1
		if sum(no_types) < 2:
			continue
				
		ofstream.write(str(no_atoms) + " " + ' '.join(map(str,no_types)) + " " + repo_files[i] + '\n')

	ofstream.close()
	
	helpers.run_bash_cmnd_to_file("xyzlist.dat", "sort tmp")
	helpers.run_bash_cmnd("rm -f tmp")
	

	
	################################
	# X. Get a list of all loose species
	################################
	
	#ofstream = open("ts_xyzlist.dat",'w')
	ofstream = open("tmp",'w')

	repo_files = sorted(glob.glob(CFG_REPO + "/ts*.wrap*xyz"))
	
	for i in xrange(len(repo_files)):
	
		no_atoms = int(helpers.head(repo_files[i],1)[0].rstrip())
		no_types = [0]*len(atm_types)
		
		for j in xrange(len(atm_types)):

			with open(repo_files[i], "r") as ifstream:
			
				for line in ifstream:
				
					if atm_types[j] in line:
					
						no_types[j] += 1
		if sum(no_types) < 2:
			continue

		ofstream.write(str(no_atoms) + " " + ' '.join(map(str,no_types)) + " " + repo_files[i] + '\n')

	ofstream.close()
	
	helpers.run_bash_cmnd_to_file("ts_xyzlist.dat", "sort tmp")
	helpers.run_bash_cmnd("rm -f tmp")	


def generate_clusters(**kwargs):

	""" 
	
	Extracts tight- and loose-clusters from a trajectory file.
	
	Usage: generate_clusters(<arguements>)
	
	Notes: See function definition in helpers.py for a full list of options. 

	WARNING: The cluster is hard coded for two atom types only, at the moment.
	      	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	default_keys   = [""]*6
	default_values = [""]*6
	
	# Clustering-specific controls
	
	default_keys[0 ] = "my_dir"	 ; default_values[0 ] = ""								# Current MD simulation
	default_keys[1 ] = "traj_file"   ; default_values[1 ] = "traj+box_250F.xyz"						# Trajectory file on which to run clustering
	default_keys[2 ] = "tight_crit"  ; default_values[2 ] = "../../../../tight_bond_crit.dat"				# File with tight bonding criteria for clustering
	default_keys[3 ] = "loose_crit"  ; default_values[3 ] = "../../../../loose_bond_crit.dat"				# File with loose bonding criteria for clustering
	default_keys[4 ] = "clu_code"	 ; default_values[4 ] = "/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/new_ts_clu.cpp"	# Clustering code
	default_keys[5 ] = "compilation" ; default_values[5 ] = "g++ -std=c++11 -O3"						# Command to compile clustering code

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	
	################################
	# 0. Set up the extraction
	################################
	
	if args["my_dir"]:
		os.chdir(args["my_dir"])

	# Compile the cluster code

	print helpers.run_bash_cmnd(args["compilation"] + " -o new_ts_clu " + args["clu_code"])

	# Dermine trajectory properties
	
	atoms = int(helpers.head(args["traj_file"],1)[0].rstrip())
		
	nframe = int(helpers.wc_l(args["traj_file"])) / (atoms+2)
	
	################################
	# 1. Run the extraction, move results to CFG_REPO
	################################
	
	print "Running the extraction code..."

	print helpers.run_bash_cmnd("./new_ts_clu " + args["traj_file"] + " " + `nframe` + " " + args["tight_crit"] + " " + args["loose_crit"])
	
	print "Code run."
	
	helpers.run_bash_cmnd("rm -rf CFG_REPO")
	helpers.run_bash_cmnd("mkdir CFG_REPO")
	
	# Save
	
	items = ' '.join(glob.glob("*.wrapped.*xyz"))
	
	helpers.run_bash_cmnd("mv " + items + " CFG_REPO")
	
	# Cleanup
	
	items =  ' '.join(glob.glob("*wrap*xyz")) + " " + ' '.join(glob.glob("*wrap*lammpstrj")) + " " + ' '.join(glob.glob("*cluster*stats"))

	helpers.run_bash_cmnd("rm -f " + items)
	

	if args["my_dir"]:
		os.chdir("..")


def get_pared_trajs(do_cluster):

	""" 
	
	Pares a trajectory down to a smaller number of frames.
	
	Usage: get_pared_trajs(True)
	
	Notes: This is run from inside of the CASE/INDEP folder.
	      	
	"""

	do_cluster = bool(do_cluster)
	
	################################
	# 0. Convert traj.gen to traj.xyz
	################################	

	#helpers.run_bash_cmnd("cd INDEP-" + args["my_indep"])
	
	frames = helpers.run_bash_cmnd("grep Step traj.gen").count("Step")

	helpers.dftbgen_to_xyz(frames, "traj.gen")
	
	# Replace comment line with box lengths
	
	box = helpers.run_bash_cmnd("head -n 1 traj.box")
	
	ifstream = open("traj.xyz"    ,'r')	
	ofstream = open("traj+box.xyz",'w')
	
	
	with open("traj.xyz") as ifstream:
		for line in ifstream:
		
			if not "Frame" in line:
				ofstream.write(line)
			else:
				ofstream.write(box)
	ofstream.close()
	
	
	################################
	# 1. Generate 20 evenly spaced frames
	################################
	
	print "Generating 20 evenly spaced frames"
	
	skip = int(frames / 20)
	
	helpers.break_apart_xyz(frames, "traj+box.xyz", skip, True)
	
	# Cat the resulting files, add penalty frames
	
	helpers.cat_pattern("traj_20F.xyz", "traj+box_#*")
	
	print "Adding any available traj_bad_r.lt.rin+dp.xyz frames"
	
	if os.path.isfile("traj_bad_r.lt.rin+dp.xyz"):
	
		helpers.cat_specific("tmp", ["traj_20F.xyz", "traj_bad_r.lt.rin+dp.xyz"])
	
		helpers.run_bash_cmnd("mv tmp traj_20F.xyz")
		
	# Clean up
	
	helpers.run_bash_cmnd("rm -f " + ' '.join(glob.glob("traj+box_#*")))
	helpers.run_bash_cmnd("rm -f " + ' '.join(glob.glob("*FORCES*")))

	################################
	# 2. If requested, generate 250 evenly spaced frames
	################################	
	
	if do_cluster:
		
		print "Generating 250 evenly spaced frames"

		skip = int(frames / 250)

		helpers.break_apart_xyz(frames, "traj+box.xyz", skip, True)

		# Cat the resulting files, add penalty frames

		helpers.cat_pattern("traj_250F.xyz", "traj+box_#*")
		
		print "Adding any available traj_bad_r.lt.rin+dp.xyz frames"
		
		if os.path.isfile("traj_bad_r.lt.rin+dp.xyz"):
		
			helpers.cat_specific("tmp", ["traj_250F.xyz", "traj_bad_r.lt.rin+dp.xyz"])

			helpers.run_bash_cmnd("mv tmp traj_250F.xyz")
			
		# Clean up

		helpers.run_bash_cmnd("rm -f " + ' '.join(glob.glob("traj+box_#*")))
		helpers.run_bash_cmnd("rm -f " + ' '.join(glob.glob("*FORCES*")))
		
	#helpers.run_bash_cmnd("cd ..")
	

if __name__=='__main__':

	""" 
	
	Allows commandline calls to generate_clusters() and get_pared_trajs().
	
	      	
	"""

	if   sys.argv[1] == "generate_clusters":
	
		generate_clusters(
			traj_file   = sys.argv[2],
			tight_crit  = sys.argv[3],
			loose_crit  = sys.argv[4],
			clu_code    = sys.argv[5],
			compilation = str(' '.join(sys.argv[6:])))
		
	elif sys.argv[1] == "get_pared_trajs":
	
		get_pared_trajs(bool(sys.argv[2]))
	
	else:
		print "ERROR: Unknown option in call to cluster.py: ", sys.argv[1]
		exit()











