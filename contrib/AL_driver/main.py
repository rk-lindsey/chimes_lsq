# Global (python) modules

import os
import sys

# Local modules

import config  # User-specified "global" vars
import helpers
import gen_ff
import run_md
import cluster
import gen_selections
import vasp_driver


""" 

Active Learning Driver Main.

Usage: unbuffer python ../PYTHON_DRIVER/main.py 0 1 2 | driver.log 

Notes: Second argument is a list of cycles to run (i.e. ALC-0, ALC-, ALC-2)
       Still working on support for parallel learning
       Still working on support for DLARS
       Still working on restart functionality
       
NOTE TO SELF: (4/15/19, 12:00 pm): The code has been set up for parallel edits, 
              but none of these new changes have been tested
	
"""	       
	


################################
# Read commandline args
################################

#THIS_ALC   = int(sys.argv[1])
THIS_CASE  = 0
THIS_INDEP = 0


################################
# Pre-process user specified variables
################################

config.ATOM_TYPES.sort()

config.CHIMES_SOLVER  = config.HPC_PYTHON + " " + config.CHIMES_SOLVER
config.CHIMES_POSTPRC = config.HPC_PYTHON + " " + config.CHIMES_POSTPRC

config.VASP_POSTPRC   = config.HPC_PYTHON + " " + config.VASP_POSTPRC


################################
################################
# Begin Active Learning
################################
################################

print "Will run for the following active learning cycles:"

for THIS_ALC in sys.argv[1:]:
	print THIS_ALC
	
helpers.run_bash_cmnd("rm -f restart.dat")

for THIS_ALC in sys.argv[1:]:

	THIS_ALC = int(THIS_ALC)
	
	# Prepare the restart file

	
	restream = open("restart.dat",'a',0)	
	restream.write("ALC: " + str(THIS_ALC) + '\n')
		
	print "Working on ALC:", THIS_ALC

	os.chdir(config.WORKING_DIR)

	# Begins in the working directory (WORKING_DIR)

	if THIS_ALC == 0: # Then this is the first ALC, so we need to do things a bit differently ... "FULLY" TESTED

		# Set up/move into the ALC directory
		
		helpers.run_bash_cmnd("rm -rf ALC-" + str(THIS_ALC))
		helpers.run_bash_cmnd("mkdir  ALC-" + str(THIS_ALC))
		
		os.chdir("ALC-" + str(THIS_ALC))
		
		################################
		# Generate the force field	
		################################
		
		active_job = gen_ff.build_amat(THIS_ALC,
					prev_gen_path      = config.ALC0_FILES,
					job_ppn            = str(config.HPC_PPN),
					job_walltime       = "1",			
					job_account        = config.HPC_ACCOUNT, 
					job_system         = config.HPC_SYSTEM,
					job_executable     = config.CHIMES_LSQ)
					
		helpers.wait_for_job(active_job, job_system = config.HPC_SYSTEM, verbose = True, job_name = "build_amat")
		
		restream.write("BUILD_AMAT: COMPLETE" + '\n')	
		
		active_job = gen_ff.solve_amat(THIS_ALC, 
					weights_force  = config.WEIGHTS_FORCE,
					weights_energy = config.WEIGHTS_ENER,
					regression_alg = config.REGRESS_ALG,
					regression_var = config.REGRESS_VAR,							
					job_ppn        = str(config.HPC_PPN),
					job_walltime   = "1",					
					job_account    = config.HPC_ACCOUNT, 
					job_system     = config.HPC_SYSTEM,
					job_executable = config.CHIMES_SOLVER)	
					
		helpers.wait_for_job(active_job, job_system = config.HPC_SYSTEM, verbose = True, job_name = "solve_amat")
		
		restream.write("SOLVE_AMAT: COMPLETE" + '\n')	
		
		helpers.run_bash_cmnd(config.CHIMES_POSTPRC + " GEN_FF/params.txt")
		
		
		################################				
		# Extract/process/select clusters
		################################
		
		# Get a list of files from which to extract clusters
		
		traj_files = helpers.cat_to_var("GEN_FF/traj_list.dat")[1:]
		
		
		# Extract clusters from each file, save into own repo, list
		
		cat_xyzlist_cmnd    = ""
		cat_ts_xyzlist_cmnd = ""
		
		for i in xrange(len(traj_files)):
		
			# Pre-process name
			
			traj_files[i] = traj_files[i].split()[1]
			
			print "Extracting clusters from file: ", traj_files[i]
			
			# Extract
		
			cluster.generate_clusters(
						traj_file   = "GEN_FF/" + traj_files[i].split()[0],
						tight_crit  = config.TIGHT_CRIT,
						loose_crit  = config.LOOSE_CRIT,
						clu_code    = config.CLU_CODE,
						compilation = "g++ -std=c++11 -O3")
			
			repo = "CFG_REPO-" + traj_files[i].split()[0]
			
			helpers.run_bash_cmnd("mv CFG_REPO " + repo)
			
			# list
		
			cluster.list_clusters(repo, 
						config.ATOM_TYPES)
						
			helpers.run_bash_cmnd("mv xyzlist.dat    " + traj_files[i].split()[0] + ".xyzlist.dat"   )
			helpers.run_bash_cmnd("mv ts_xyzlist.dat " + traj_files[i].split()[0] + ".ts_xyzlist.dat")
			
			cat_xyzlist_cmnd    += traj_files[i].split()[0] + ".xyzlist.dat "
			cat_ts_xyzlist_cmnd += traj_files[i].split()[0] + ".ts_xyzlist.dat "
			
		helpers.cat_specific("xyzlist.dat"   , cat_xyzlist_cmnd   .split())
		helpers.cat_specific("ts_xyzlist.dat", cat_ts_xyzlist_cmnd.split())
		
		helpers.run_bash_cmnd("rm -f " + cat_xyzlist_cmnd   )
		helpers.run_bash_cmnd("rm -f " + cat_ts_xyzlist_cmnd)
		
		restream.write("CLUSTER_EXTRACTION: COMPLETE" + '\n')
		
		
		# Compute cluster energies
		
		active_jobs = cluster.get_repo_energies(
					base_runfile   = config.ALC0_FILES + "/run_md.base",
					driver_dir     = config.DRIVER_DIR,
					job_ppn        = str(config.HPC_PPN),
					job_walltime   = "1",					
					job_account    = config.HPC_ACCOUNT, 
					job_system     = config.HPC_SYSTEM,
					job_executable = config.CHIMES_MD)	
					
		helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "get_repo_energies")
		
		print helpers.run_bash_cmnd("pwd")
		print helpers.run_bash_cmnd("ls -lrt")	
		
		restream.write("CLUENER_CALC: COMPLETE" + '\n')
		
		
		# Generate cluster sub-selection and store in central repository
		
		gen_selections.cleanup_repo(THIS_ALC)
		
		gen_selections.gen_subset(
					nsel     = config.MEM_NSEL, # Number of selections to make    
					nsweep   = config.MEM_CYCL, # Number of MC sqeeps	      
					nbins    = config.MEM_BINS, # Number of histogram bins  	
					ecut     = config.MEM_ECUT) # Maximum energy to consider	
					
		restream.write("CLU_SELECTION: COMPLETE" + '\n')
		
		gen_selections.populate_repo(THIS_ALC)

		repo = "CASE-" + str(THIS_CASE) + "_INDEP_" + str(THIS_INDEP) + "/CFG_REPO/"

		
		################################
		# Launch VASP
		################################
		
		vasp_driver.cleanup_and_setup(["all"], ".")
		
		restream.write("CLEANSETUP_VASP: COMPLETE" + '\n')		
		
		active_jobs = vasp_driver.setup_vasp(THIS_ALC,
						["all"], 
						config.ATOM_TYPES, 
						basefile_dir   = config.VASP_FILES,
						traj_list      = config.ALC0_FILES + "/traj_list.dat", # Has a temperature for each file ... expected as integer
						vasp_exe       = config.VASP_EXE,
						job_nodes      = config.VASP_NODES,
						job_ppn        = config.HPC_PPN,
						job_walltime   = config.VASP_WALLT,
						job_queue      = config.VASP_QUEUE,
						job_account    = "pbronze",
						job_system     = config.HPC_SYSTEM)
		
		helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "setup_vasp")
		
		restream.write("INIT_VASPJOB: COMPLETE" + '\n')
		
		# Check that the job was complete
		
		while True:

			active_jobs = vasp_driver.continue_job(
						["all"], 
						job_system     = config.HPC_SYSTEM)
						
			print "active jobs: ", active_jobs			
						
			if len(active_jobs) > 0:
			
				print "waiting for restarted vasp job."
			
				helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "setup_vasp - restarts")
			else:
				print "All jobs are complete"
				break		
		
		
		restream.write("ALL_VASPJOBS: COMPLETE" + '\n')

		
		# Post-process the vasp jobs
		
		print "post-processing..."	
		
		vasp_driver.post_process(["all"], "ENERGY",
					vasp_postproc = config.VASP_POSTPRC)

		os.chdir("..")
		
		print "ALC-", THIS_ALC, "is complete"	
		
		restream.write("THIS_ALC: COMPLETE" + '\n')
		restream.close()
		
		
	else:

		# Set up/move into the ALC directory
		
		helpers.run_bash_cmnd("rm -rf ALC-" + str(THIS_ALC))
		helpers.run_bash_cmnd("mkdir  ALC-" + str(THIS_ALC))
		
		os.chdir("ALC-" + str(THIS_ALC))
		
		vasp_all_path = config.WORKING_DIR + "/ALC-" + `THIS_ALC-1` + "/VASP-all/"
		vasp_20F_path = ""
		
		if THIS_ALC > 1:
			vasp_20F_path = config.WORKING_DIR + "/ALC-" + `THIS_ALC-1` + "/VASP-20/"
		
		active_job = gen_ff.build_amat(THIS_ALC, 
				prev_vasp_all_path = vasp_all_path,
				prev_vasp_20_path  = vasp_20F_path,
				job_ppn        = str(config.HPC_PPN),
				job_walltime   = "1",			
				job_account    = config.HPC_ACCOUNT, 
				job_system     = config.HPC_SYSTEM,
				job_executable = config.CHIMES_LSQ)
		
		helpers.wait_for_job(active_job, job_system = config.HPC_SYSTEM, verbose = True, job_name = "build_amat")
		
		restream.write("BUILD_AMAT: COMPLETE" + '\n')	
		
		active_job = gen_ff.solve_amat(THIS_ALC, 
				weights_force  = config.WEIGHTS_FORCE,
				weights_energy = config.WEIGHTS_ENER,
				regression_alg = config.REGRESS_ALG,
				regression_var = config.REGRESS_VAR,							
				job_ppn        = str(config.HPC_PPN),
				job_walltime   = "1",					
				job_account    = config.HPC_ACCOUNT, 
				job_system     = config.HPC_SYSTEM,
				job_executable = config.CHIMES_SOLVER)	
				
		helpers.wait_for_job(active_job, job_system = config.HPC_SYSTEM, verbose = True, job_name = "solve_amat")	
		
		restream.write("SOLVE_AMAT: COMPLETE" + '\n')
		
		helpers.run_bash_cmnd(config.CHIMES_POSTPRC + " GEN_FF/params.txt")
		
		################################				
		# Run MD
		################################
		
		
		# ... May want to consider making speciation optional ... can add another key word that allows the user to set up different 
		#    types of post-processing jobs
		
				
		# Run the MD/cluster jobs
		
		active_jobs = []
		
		for THIS_CASE in xrange(config.NO_CASES):		
		
		
			active_job = run_md.run_md(THIS_ALC, THIS_CASE, THIS_INDEP,
					"C1 O1 1(O-C)",
					"C1 O2 2(O-C)",
					"C2 O2 1(C-C) 2(O-C)",
					"C3 O2 2(C-C) 2(O-C)",
					basefile_dir   = config.CHIMES_MDFILES, 
					driver_dir     = config.DRIVER_DIR,
					penalty_pref   = 1.0E6,		
					penalty_dist   = 0.02, 		
					molanal_dir    = config.CHIMES_MOLANAL, 
					local_python   = config.HPC_PYTHON, 	
					do_cluster     = config.DO_CLUSTER,	
					tight_crit     = config.TIGHT_CRIT,	
					loose_crit     = config.LOOSE_CRIT,	
					clu_code       = config.CLU_CODE,  	
					compilation    = "g++ -std=c++11 -O3",
					job_name       = "ALC-"+ str(THIS_ALC) +"-md-c" + str(THIS_CASE) +"-i" + str(THIS_INDEP),
					job_nodes      = 8,	   	 
					job_ppn        = 36,	   	 
					job_walltime   = 1,	   	 
					job_queue      = "pdebug", 	 
					job_account    = "pbronze",	 
					job_executable = config.CHIMES_MD,	 
					job_system     = "slurm",  	 
					job_file       = "run.cmd")	

		active_jobs.append(active_job.split()[0])										
		helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "run_md")

		restream.write("RUN_MD: COMPLETE" + '\n')
		
		
		# list ... remember, we only do clustering/active learning on a single indep (0)
		
		cat_xyzlist_cmnd    = ""
		cat_ts_xyzlist_cmnd = ""		
		
		for THIS_CASE in xrange(config.NO_CASES):
		        
		        repo = "CASE-" + str(THIS_CASE) + "_INDEP_" + str(THIS_INDEP) + "/CFG_REPO/"
		        
		        cluster.list_clusters(repo, 
		        		config.ATOM_TYPES)	
		        		
		        helpers.run_bash_cmnd("mv xyzlist.dat	 " + "CASE-" + str(THIS_CASE) + ".xyzlist.dat"   )
		        helpers.run_bash_cmnd("mv ts_xyzlist.dat " + "CASE-" + str(THIS_CASE) + ".ts_xyzlist.dat")
		        
		        cat_xyzlist_cmnd    += "CASE-" + str(THIS_CASE) + ".xyzlist.dat "
		        cat_ts_xyzlist_cmnd += "CASE-" + str(THIS_CASE) + ".ts_xyzlist.dat "
		        
		helpers.cat_specific("xyzlist.dat"   , cat_xyzlist_cmnd   .split())
		helpers.cat_specific("ts_xyzlist.dat", cat_ts_xyzlist_cmnd.split())
		
		helpers.run_bash_cmnd("rm -f " + cat_xyzlist_cmnd   )
		helpers.run_bash_cmnd("rm -f " + cat_ts_xyzlist_cmnd)
		
		
		restream.write("CLUSTER_EXTRACTION: COMPLETE" + '\n')					
		
		# Compute cluster energies
		
		gen_selections.cleanup_repo(THIS_ALC)	
		
		active_jobs = cluster.get_repo_energies(
				calc_central   = True,
				base_runfile   = config.CHIMES_MDFILES + "/" + "run_md.base",
				driver_dir     = config.DRIVER_DIR,
				job_ppn        = str(config.HPC_PPN),
				job_walltime   = "1",					
				job_account    = config.HPC_ACCOUNT, 
				job_system     = config.HPC_SYSTEM,
				job_executable = config.CHIMES_MD)	
				
		helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "get_repo_energies")

		restream.write("CLUENER_CALC: COMPLETE" + '\n')	


		# Generate cluster sub-selection and store in central repository

		gen_selections.gen_subset(
				repo     = "../CENTRAL_REPO/full_repo.energies_normed",
				nsel     = config.MEM_NSEL, # Number of selections to make    
				nsweep   = config.MEM_CYCL, # Number of MC sqeeps	      
				nbins    = config.MEM_BINS, # Number of histogram bins  	
				ecut     = config.MEM_ECUT) # Maximum energy to consider	
				
		restream.write("CLU_SELECTION: COMPLETE" + '\n')
						


		gen_selections.populate_repo(THIS_ALC)	


		################################
		# Launch VASP
		################################
		
		vasp_driver.cleanup_and_setup(["20", "all"], "..")
		
		restream.write("CLEANSETUP_VASP: COMPLETE" + '\n')
		
		active_jobs = vasp_driver.setup_vasp(THIS_ALC,
						["20", "all"], 
						config.ATOM_TYPES,
						THIS_CASE, 
						config.THIS_SMEAR,
						basefile_dir   = config.VASP_FILES,
						build_dir      = "..", # Put the VASP-x directories in the ALC-X folder
						vasp_exe       = config.VASP_EXE,
						job_nodes      = config.VASP_NODES,
						job_ppn        = config.HPC_PPN,
						job_walltime   = config.VASP_WALLT,
						job_queue      = config.VASP_QUEUE,
						job_account    = "pbronze",
						job_system     = config.HPC_SYSTEM)
		
		helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "setup_vasp")
		
		restream.write("INIT_VASPJOB: COMPLETE" + '\n')
		
		# Check that the job was complete
		
		while True:

			active_jobs = vasp_driver.continue_job(
						["all","20"], 
						job_system     = config.HPC_SYSTEM)
						
			print "active jobs: ", active_jobs			
						
			if len(active_jobs) > 0:
			
				print "waiting for restarted vasp job."
			
				helpers.wait_for_jobs(active_jobs, job_system = config.HPC_SYSTEM, verbose = True, job_name = "setup_vasp - restarts")
			else:
				print "All jobs are complete"
				break	
				
		restream.write("ALL_VASPJOBS: COMPLETE" + '\n')	
		
		
		# Post-process the vasp jobs
		
		print "post-processing..."	
		
		vasp_driver.post_process(["all","20"], "ENERGY",
					vasp_postproc = config.VASP_POSTPRC)
					
		os.chdir("..")
		
		print "ALC-", THIS_ALC, "is complete"	
		
		restream.write("THIS_ALC: COMPLETE" + '\n')
		restream.close()					
