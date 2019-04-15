import os
import sys
import helpers
import gen_ff
import run_md
import cluster
import gen_selections
import vasp_driver


# !!!! TO DO: Add cleanup for central repository (ALC-X.all_selections.xyzlist removed if not yet run, full_repo.xyzlist rebuilt)


################################
# Set user specified variables ... should this be read from an input file?
################################

##### General variables

ATOM_TYPES    = ["C", "O"]
NO_CASES      = 1

WORKING_DIR   = "/p/lustre1/rlindsey/RC4B_RAG/11-12-18/TEST_PYTHON_DRIVER-A/"	# Directory from which all ALCs will be run
DRIVER_DIR    = "/p/lustre1/rlindsey/RC4B_RAG/11-12-18/PYTHON_DRIVER/"		# This driver's src location
CHIMES_SRCDIR = "/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/SVN_SRC-11-12-18/"	# Location of ChIMES src files ... expects also post_proc_lsq2.py there


##### General HPC paths

HPC_PPN       = 36
HPC_ACCOUNT   = "pbronze"
HPC_SYSTEM    = "slurm"
HPC_PYTHON    = "/usr/tce/bin/python"


##### ChIMES LSQ and MD paths

ALC0_FILES    = WORKING_DIR + "ALC-0_BASEFILES/"

CHIMES_LSQ    = CHIMES_SRCDIR + "chimes_lsq"
CHIMES_SOLVER = CHIMES_SRCDIR + "lsq2.py"
CHIMES_POSTPRC= CHIMES_SRCDIR + "post_proc_lsq2.py"

WEIGHTS_FORCE = 1.0 
WEIGHTS_ENER  = 5.0

REGRESS_ALG   = "lassolars"
REGRESS_VAR   = "1.0E-4"

CHIMES_MD     = CHIMES_SRCDIR + "chimes_md"
CHIMES_MOLANAL= "/g/g17/rlindsey/CURR_TRACKED-GEN/contrib/molanlal/"  
CHIMES_MDFILES= WORKING_DIR + "CHIMESMD_BASEFILES/"

##### Cluster specific paths/variables

DO_CLUSTER    = True
TIGHT_CRIT    = WORKING_DIR + "tight_bond_crit.dat"
LOOSE_CRIT    = WORKING_DIR + "loose_bond_crit.dat"
CLU_CODE      = "/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/new_ts_clu.cpp"


##### ALC-specific variables ... Note: Hardwired for partial memory mode

MEM_BINS = 40
MEM_CYCL = MEM_BINS/10
MEM_NSEL = 400
MEM_ECUT = 100.0


##### VASP Specific variables

VASP_FILES   = WORKING_DIR + "VASP_BASEFILES"
VASP_POSTPRC = CHIMES_SRCDIR + "vasp2xyzf.py"
VASP_NODES   = 6
VASP_WALLT   = 1
VASP_QUEUE   = "pdebug"
VASP_EXE     = "/usr/gapps/emc-vasp/vasp.5.4.1/build/std/vasp"
VASP_MODULES = "mkl"


################################
# Pre-process user specified variables
################################

ATOM_TYPES.sort()

CHIMES_SOLVER  = HPC_PYTHON + " " + CHIMES_SOLVER
CHIMES_POSTPRC = HPC_PYTHON + " " + CHIMES_POSTPRC

VASP_POSTPRC   = HPC_PYTHON + " " + VASP_POSTPRC


################################
# Read commandline args
################################

#THIS_ALC   = int(sys.argv[1])
THIS_CASE  = 0
THIS_SMEAR = 2400
THIS_INDEP = 0


################################
################################
# Begin Active Learning
################################
################################

print "Will run for the following active learning cycles:"

for THIS_ALC in sys.argv[1:]:
	print THIS_ALC
	
################################
# Set up the next ALC
################################


for THIS_ALC in sys.argv[1:]:

	THIS_ALC = int(THIS_ALC)
	
	print "Working on ALC:", THIS_ALC
	print "Returning to base working directory:", WORKING_DIR
	
	os.chdir(WORKING_DIR)
	
	print helpers.run_bash_cmnd("pwd")

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
					prev_gen_path      = ALC0_FILES,
					job_ppn            = str(HPC_PPN),
					job_walltime       = "1",			
					job_account        = HPC_ACCOUNT, 
					job_system         = HPC_SYSTEM,
					job_executable     = CHIMES_LSQ)
					
		helpers.wait_for_job(active_job, job_system = HPC_SYSTEM, verbose = True, job_name = "build_amat")	
		
		active_job = gen_ff.solve_amat(THIS_ALC, 
					weights_force  = WEIGHTS_FORCE,
					weights_energy = WEIGHTS_ENER,
					regression_alg = REGRESS_ALG,
					regression_var = REGRESS_VAR,							
					job_ppn        = str(HPC_PPN),
					job_walltime   = "1",					
					job_account    = HPC_ACCOUNT, 
					job_system     = HPC_SYSTEM,
					job_executable = CHIMES_SOLVER)	
					
		helpers.wait_for_job(active_job, job_system = HPC_SYSTEM, verbose = True, job_name = "solve_amat")	
		
		helpers.run_bash_cmnd(CHIMES_POSTPRC + " GEN_FF/params.txt")
		
		
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
						tight_crit  = TIGHT_CRIT,
						loose_crit  = LOOSE_CRIT,
						clu_code    = CLU_CODE,
						compilation = "g++ -std=c++11 -O3")
			
			repo = "CFG_REPO-" + traj_files[i].split()[0]
			
			helpers.run_bash_cmnd("mv CFG_REPO " + repo)
			
			# list
		
			cluster.list_clusters(repo, 
						ATOM_TYPES)
						
			helpers.run_bash_cmnd("mv xyzlist.dat    " + traj_files[i].split()[0] + ".xyzlist.dat"   )
			helpers.run_bash_cmnd("mv ts_xyzlist.dat " + traj_files[i].split()[0] + ".ts_xyzlist.dat")
			
			cat_xyzlist_cmnd    += traj_files[i].split()[0] + ".xyzlist.dat "
			cat_ts_xyzlist_cmnd += traj_files[i].split()[0] + ".ts_xyzlist.dat "
			
		helpers.cat_specific("xyzlist.dat"   , cat_xyzlist_cmnd   .split())
		helpers.cat_specific("ts_xyzlist.dat", cat_ts_xyzlist_cmnd.split())
		
		helpers.run_bash_cmnd("rm -f " + cat_xyzlist_cmnd   )
		helpers.run_bash_cmnd("rm -f " + cat_ts_xyzlist_cmnd)
		
		
		# Compute cluster energies
		
		active_jobs = cluster.get_repo_energies(
					base_runfile   = ALC0_FILES + "/run_md.base",
					driver_dir     = DRIVER_DIR,
					job_ppn        = str(HPC_PPN),
					job_walltime   = "1",					
					job_account    = HPC_ACCOUNT, 
					job_system     = HPC_SYSTEM,
					job_executable = CHIMES_MD)	
					
		helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "get_repo_energies")
		
		print helpers.run_bash_cmnd("pwd")
		print helpers.run_bash_cmnd("ls -lrt")	
		
		
		# Generate cluster sub-selection and store in central repository
		
		gen_selections.gen_subset(
					nsel     = MEM_NSEL, # Number of selections to make    
					nsweep   = MEM_CYCL, # Number of MC sqeeps	      
					nbins    = MEM_BINS, # Number of histogram bins  	
					ecut     = MEM_ECUT) # Maximum energy to consider	
		
		gen_selections.populate_repo(THIS_ALC)
		
		

		repo = "CASE-" + str(THIS_CASE) + "_INDEP_" + str(THIS_INDEP) + "/CFG_REPO/"

		
		################################
		# Launch VASP
		################################
		
		active_jobs = vasp_driver.setup_vasp(THIS_ALC,
						["all"], 
						ATOM_TYPES, 
						basefile_dir   = VASP_FILES,
						traj_list      = ALC0_FILES + "/traj_list.dat", # Has a temperature for each file ... expected as integer
						vasp_exe       = VASP_EXE,
						job_nodes      = VASP_NODES,
						job_ppn        = HPC_PPN,
						job_walltime   = VASP_WALLT,
						job_queue      = VASP_QUEUE,
						job_account    = "pbronze",
						job_system     = HPC_SYSTEM)
		
		helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "setup_vasp")
		
		# Check that the job was complete
		
		while True:

			active_jobs = vasp_driver.continue_job(
						["all"], 
						job_system     = HPC_SYSTEM)
						
			print "active jobs: ", active_jobs			
						
			if len(active_jobs) > 0:
			
				print "waiting for restarted vasp job."
			
				helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "setup_vasp - restarts")
			else:
				print "All jobs are complete"
				break		
		
		


		
		# Post-process the vasp jobs
		
		print "post-processing..."	
		
		vasp_driver.post_process(["all"], "ENERGY",
					vasp_postproc = VASP_POSTPRC)

		os.chdir("..")
		
		print "ALC-", THIS_ALC, "is complete"	
		
		
	else:
			
		# Set up/move into the ALC directory
		
		helpers.run_bash_cmnd("rm -rf ALC-" + str(THIS_ALC))
		helpers.run_bash_cmnd("mkdir  ALC-" + str(THIS_ALC))
		
		os.chdir("ALC-" + str(THIS_ALC))
		
		vasp_all_path = WORKING_DIR + "/ALC-" + `THIS_ALC-1` + "/VASP-all/"
		vasp_20F_path = ""
		
		if THIS_ALC > 1:
			vasp_20F_path = WORKING_DIR + "/ALC-" + `THIS_ALC-1` + "/VASP-20/"
		
		active_job = gen_ff.build_amat(THIS_ALC, 
				prev_vasp_all_path = vasp_all_path,
				prev_vasp_20_path  = vasp_20F_path,
				job_ppn        = str(HPC_PPN),
				job_walltime   = "1",			
				job_account    = HPC_ACCOUNT, 
				job_system     = HPC_SYSTEM,
				job_executable = CHIMES_LSQ)
		
		helpers.wait_for_job(active_job, job_system = HPC_SYSTEM, verbose = True, job_name = "build_amat")	
		
		active_job = gen_ff.solve_amat(THIS_ALC, 
				weights_force  = WEIGHTS_FORCE,
				weights_energy = WEIGHTS_ENER,
				regression_alg = REGRESS_ALG,
				regression_var = REGRESS_VAR,							
				job_ppn        = str(HPC_PPN),
				job_walltime   = "1",					
				job_account    = HPC_ACCOUNT, 
				job_system     = HPC_SYSTEM,
				job_executable = CHIMES_SOLVER)	
				
		helpers.wait_for_job(active_job, job_system = HPC_SYSTEM, verbose = True, job_name = "solve_amat")	
		
		
		
		helpers.run_bash_cmnd(CHIMES_POSTPRC + " GEN_FF/params.txt")
		
		################################				
		# Run MD
		################################

		# ...time to start working on the md part, with considerations for parallel learning
		# ...Should expect something like: (among also  bonds.dat)
		#	CHIMES_MDFILES= "/p/lustre1/rlindsey/RC4B_RAG/11-12-18/2400K_TRAJ-REFINED-ALC-40BINS/CHIMES-MD_BASEFILES"
		#	...where:
		#		.../CHIMES-MD_BASEFILES/case-1.indep-1.run_md.in # differ by seed and input file name
		#		.
		#		.
		#		.
		#		.../CHIMES-MD_BASEFILES/case-n.indep-m.run_md.in
		#	...and similarly:
		#		.../CHIMES-MD_BASEFILES/case-n.indep-m.input.xyz
		#
		#
		# ... May want to consider making speciation optional ... can add another key word that allows the user to set up different 
		#    types of post-processing jobs

		#active_jobs = []
		#for c in xrange(ncases)):
		#	
		#	for i in xrange(nindep)):

		active_job = run_md.run_md(THIS_ALC, THIS_CASE, THIS_INDEP,
					"C1 O1 1(O-C)",
					"C1 O2 2(O-C)",
					"C2 O2 1(C-C) 2(O-C)",
					"C3 O2 2(C-C) 2(O-C)",
					basefile_dir   = CHIMES_MDFILES, 
					driver_dir     = DRIVER_DIR,
					penalty_pref   = 1.0E6,		
					penalty_dist   = 0.02, 		
					molanal_dir    = CHIMES_MOLANAL, 
					local_python   = HPC_PYTHON, 	
					do_cluster     = DO_CLUSTER,	
					tight_crit     = TIGHT_CRIT,	
					loose_crit     = LOOSE_CRIT,	
					clu_code       = CLU_CODE,  	
					compilation    = "g++ -std=c++11 -O3",
					job_name       = "ALC-"+ str(THIS_ALC) +"-md-c" + str(THIS_CASE) +"-i" + str(THIS_INDEP),
					job_nodes      = 8,	   	 
					job_ppn        = 36,	   	 
					job_walltime   = 1,	   	 
					job_queue      = "pdebug", 	 
					job_account    = "pbronze",	 
					job_executable = CHIMES_MD,	 
					job_system     = "slurm",  	 
					job_file       = "run.cmd")	

		#active_jobs.append(str(active_job).split()[0])										
		#helpers.wait_for_jobs(active_job, job_system = HPC_SYSTEM) # This is untested!					

		helpers.wait_for_job(active_job, job_system = HPC_SYSTEM, verbose = True, job_name = "run_md")	



		repo = "CASE-" + str(THIS_CASE) + "_INDEP_" + str(THIS_INDEP) + "/CFG_REPO/"
		base = "case-" + str(THIS_CASE) + ".indep-" + str(THIS_INDEP) + ".run_md.in"

		# list

		cluster.list_clusters(repo, 
					ATOM_TYPES)
					


		# Compute cluster energies

		active_jobs = cluster.get_repo_energies(
				calc_central   = True,
				base_runfile   = CHIMES_MDFILES + "/" + "run_md.base",
				driver_dir     = DRIVER_DIR,
				job_ppn        = str(HPC_PPN),
				job_walltime   = "1",					
				job_account    = HPC_ACCOUNT, 
				job_system     = HPC_SYSTEM,
				job_executable = CHIMES_MD)	
				
		helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "get_repo_energies")	
		
		# Generate cluster sub-selection and store in central repository
		
		gen_selections.gen_subset(
				repo     = "../CENTRAL_REPO/full_repo.energies_normed",
				nsel     = MEM_NSEL, # Number of selections to make    
				nsweep   = MEM_CYCL, # Number of MC sqeeps	      
				nbins    = MEM_BINS, # Number of histogram bins  	
				ecut     = MEM_ECUT) # Maximum energy to consider	

		gen_selections.populate_repo(THIS_ALC)	




		################################
		# Launch VASP
		################################
		
		active_jobs = vasp_driver.setup_vasp(THIS_ALC,
						["20", "all"], 
						ATOM_TYPES,
						THIS_CASE, 
						THIS_SMEAR,
						basefile_dir   = VASP_FILES,
						traj_list      = ALC0_FILES + "/traj_list.dat", # Has a temperature for each file ... expected as integer
						vasp_exe       = VASP_EXE,
						job_nodes      = VASP_NODES,
						job_ppn        = HPC_PPN,
						job_walltime   = VASP_WALLT,
						job_queue      = VASP_QUEUE,
						job_account    = "pbronze",
						job_system     = HPC_SYSTEM)
		
		helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "setup_vasp")
		
		# Check that the job was complete
		
		while True:

			active_jobs = vasp_driver.continue_job(
						["all","20"], 
						job_system     = HPC_SYSTEM)
						
			print "active jobs: ", active_jobs			
						
			if len(active_jobs) > 0:
			
				print "waiting for restarted vasp job."
			
				helpers.wait_for_jobs(active_jobs, job_system = HPC_SYSTEM, verbose = True, job_name = "setup_vasp - restarts")
			else:
				print "All jobs are complete"
				break		
		
		
		# Post-process the vasp jobs
		
		print "post-processing..."	
		
		vasp_driver.post_process(["all"], "ENERGY",
					vasp_postproc = VASP_POSTPRC)
					
		os.chdir("..")
		
		print "ALC-", THIS_ALC, "is complete"						
