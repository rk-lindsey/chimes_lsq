# Global (python) modules

import glob
import os

# Local modules

import helpers

def run_md(my_ALC, my_case, my_indep, *argv, **kwargs):

	""" 
	
	Launches a ChIMES md simulation and post-md analysis.
	
	Usage: run_md(1, 0, 0, <arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Requrires config.CHIMES_MD.
	       Requrires config.CHIMES_MOLANAL (should contain molanal.new and findmolecules.pl)
	       Expects to be called from ALC-my_ALC's base folder.
	       Assumes job is being launched from ALC-X.
	       Supports ~parallel learning~ via file structure:
	              ALC-X/CASE-1_INDEP-1/<md simulation/stuff>
	       Expects input files named like:
	              case-1.indep-1.input.xyz and case-1.indep-1.run_md.in
	       will also run molanal on finished simulation.
	       Will also post-process the molanal output.
	       Will generate clusters.
	       Will save clusters to CASE-X_INDEP-0/CFG_REPO.
	       Returns a job_id for the submitted job.
	       	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	### ...argv
	
	args_species = argv # This is a pointer!
	
	### ...kwargs
	
	default_keys   = [""]*20
	default_values = [""]*20


	# MD specific controls
	
	default_keys[0 ] = "basefile_dir"  ; default_values[0 ] = "../CHIMES-MD_BASEFILES/"	# Directory containing run_md.base, etc.
	default_keys[1 ] = "driver_dir"    ; default_values[1 ] = "" 				# Post_proc_lsq*py file... should also include the python command
	default_keys[2 ] = "penalty_pref"  ; default_values[2 ] = 1.0E6				# Penalty function pre-factor
	default_keys[3 ] = "penalty_dist"  ; default_values[3 ] = 0.02				# Pentaly function kick-in distance
	default_keys[4 ] = "molanal_dir"   ; default_values[4 ] = ""				# Path to the molanal directory
	default_keys[5 ] = "local_python"  ; default_values[5 ] = ""				# Local computer's python executable

	# Cluster specific controls
		
	default_keys[6 ] = "do_cluster " ; default_values[6 ] = True			# Local computer's python executable
	default_keys[7 ] = "tight_crit"  ; default_values[7 ] = "../../../../tight_bond_crit.dat"				# File with tight bonding criteria for clustering
	default_keys[8 ] = "loose_crit"  ; default_values[8 ] = "../../../../loose_bond_crit.dat"				# File with loose bonding criteria for clustering
	default_keys[9 ] = "clu_code"	 ; default_values[9 ] = "/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/new_ts_clu.cpp"	# Clustering code
	default_keys[10] = "compilation" ; default_values[10] = "g++ -std=c++11 -O3"						# Command to compile clustering code	

	# Overall job controls	
	
	default_keys[11] = "job_name"	   ; default_values[11] = "ALC-"+ `my_ALC`+"-md"	# Name for ChIMES md job
	default_keys[12] = "job_nodes"     ; default_values[12] = "2"				# Number of nodes for ChIMES md job
	default_keys[13] = "job_ppn"	   ; default_values[13] = "36"  			# Number of processors per node for ChIMES md job
	default_keys[14] = "job_walltime"  ; default_values[14] = "1"				# Walltime in hours for ChIMES md job
	default_keys[15] = "job_queue"     ; default_values[15] = "pdebug"			# Queue for ChIMES md job
	default_keys[16] = "job_account"   ; default_values[16] = "pbronze"			# Account for ChIMES md job
	default_keys[17] = "job_executable"; default_values[17] = ""				# Full path to executable for ChIMES md job
	default_keys[18] = "job_system"    ; default_values[18] = "slurm"			# slurm or torque       
	default_keys[19] = "job_file"	   ; default_values[19] = "run.cmd"			# Name of the resulting submit script   
	

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)	
	
	################################
	# 1. Create the MD directory for this specific case and independent simulation, grab the base files
	################################
	
	my_md_path = "CASE-" + str(my_case) + "_INDEP_" + str(my_indep) + "/"

	helpers.run_bash_cmnd("rm -rf " + my_md_path)
	helpers.run_bash_cmnd("mkdir -p " + my_md_path)
	
	helpers.run_bash_cmnd("cp "+ ' '.join(glob.glob(args["basefile_dir"] + "/*"  )) + " " + my_md_path)
	helpers.run_bash_cmnd("cp GEN_FF/params.txt.reduced " + my_md_path)
	

	################################
	# 2. Post-process the parameter file
	################################
	
	os.chdir(my_md_path)


	ifstream   = open("params.txt.reduced",'r')
	paramsfile = ifstream.readlines()

	ofstream = open("tmp",'w')
	
	found = False
		
	for i in xrange(len(paramsfile)):
	
		if found:
			ofstream.write(paramsfile[i])
			ofstream.write("PAIR CHEBYSHEV PENALTY DIST:    " + str(args["penalty_dist"]) + '\n')
			ofstream.write("PAIR CHEBYSHEV PENALTY SCALING: " + str(args["penalty_pref"]) + '\n\n')

			found = False
		else:
	
			ofstream.write(paramsfile[i])
	
			if "FCUT TYPE" in paramsfile[i]:
				found  = True
	ofstream.close()
	helpers.run_bash_cmnd("mv tmp params.txt.reduced")
	
	
	################################
	# 3. Post-process the run_md.in file
	################################
	
	md_infile  = "case-" + str(my_case) + ".indep-" + str(my_indep) + ".run_md.in"
	md_xyzfile = "case-" + str(my_case) + ".indep-" + str(my_indep) + ".input.xyz"	
	
	ifstream = open(md_infile,'r')
	runfile  = ifstream.readlines()
	
	ofstream = open("tmp",'w')	
	
	found1 = False
	found2 = False
	found3 = False
		
	for i in xrange(len(runfile)):
	
		if found1:
			ofstream.write('\t' + str(1+my_indep) + str(1+my_indep) + str(1+my_indep) + str(1+my_indep) + '\n')
			found1 = False
		elif found2:
			ofstream.write('\t' + md_xyzfile + '\n')
			found2 = False	
		elif found3:
			ofstream.write('\t' + "params.txt.reduced" + '\n')
			found3 = False		
		else:
	
			ofstream.write(runfile[i])
	
			if "RNDSEED" in runfile[i]:
				found1 = True
				
			if "CRDFILE" in runfile[i]:
				found2 = True
				
			if "PRMFILE" in runfile[i]:
				found3 = True
				
	ofstream.close()
	helpers.run_bash_cmnd("mv tmp " + md_infile)	
	
	
	################################
	# 4. Launch the md simulation
	# ... will also run molanal on finished simulation
	# ... Will also post-process the molanal output
	# ... Will generate clusters
	################################
	
	
	# Create the task string
	
	job_task  = "-n " + `int(args["job_nodes"])*int(args["job_ppn"])` + " " +  args["job_executable"] + " " + md_infile + " > run_md.out"
	job_task += '\n'
	job_task += args["molanal_dir"] + "/molanal.new traj.gen > traj.gen-molanal.out"
	job_task += '\n'		
	job_task += args["molanal_dir"] + "/findmolecules.pl traj.gen-molanal.out > traj.gen-find_molecs.out"
	job_task += '\n'
	job_task += "rm -rf molecules molanal*"	
	job_task += '\n'	
	job_task += args["local_python"] + " " + args["driver_dir"] + "/post_process_molanal.py " + " \"" + '\" \"'.join(args_species) + "\""
	job_task += '\n'		
	job_task += args["local_python"] + " " + args["driver_dir"] + "/cluster.py get_pared_trajs "      + `args["do_cluster"]`
	
	print "compilation: ", args["compilation"]
	
	if args["do_cluster"]:
		job_task += '\n'	
		job_task += args["local_python"] + " " + args["driver_dir"] + "/cluster.py generate_clusters " + " traj_250F.xyz" + " " + args["tight_crit" ] + " " + args["loose_crit" ] + " " + args["clu_code"   ] + " " + args["compilation"]
	

	if args["job_system"] == "slurm":
		job_task = "srun "   + job_task
	else:
		job_task = "mpirun " + job_task	
	
	md_jobid = helpers.create_and_launch_job(
		job_name       =     args["job_name"	] ,
		job_nodes      = str(args["job_nodes"   ]),
		job_ppn        = str(args["job_ppn"     ]),
		job_walltime   = str(args["job_walltime"]),
		job_queue      =     args["job_queue"	] ,
		job_account    =     args["job_account" ] ,
		job_executable =     job_task,
		job_system     =     args["job_system"  ] ,
		job_file       =     "run_chimesmd.cmd")
		
	
	os.chdir("..")
	
	return md_jobid	
	
