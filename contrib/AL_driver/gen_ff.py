# Global (python) modules

import os.path
import os
import glob

# Local modules

import helpers


def build_amat(my_ALC, **kwargs):  

	""" 
	
	Generates the A- and b-matrices for an input trajectory.
	
	Usage: build_amat(1, <arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Requrires config.CHIMES_LSQ.
	       Expects to be called from ALC-my_ALC's base folder.
	       Returns a job_id for the submitted job.
	       	
	"""

	################################
	# 0. Set up an argument parser
	################################
	
	default_keys   = [""]*11
	default_values = [""]*11
	
	# Paths
	
	default_keys[0 ] = "prev_gen_path"     ; default_values[0 ] = 	 "../ALC-" + `my_ALC-1` + "/GEN_FF/"    # Path to previous ALCs GEN_FF folder -- absolute is best
	default_keys[1 ] = "prev_vasp_all_path"; default_values[1 ] = 	 ""	       	       	       	        # Path to previous ALCs VASP-all folder -- absolute is best
	default_keys[2 ] = "prev_vasp_20_path" ; default_values[2 ] = 	 ""	       	       	       	        # Path to previous ALCs VASP-20 folder --absolute is best
	
	# Job controls
	
	default_keys[3 ] = "job_name"	       ; default_values[3 ] = 	 "ALC-"+ `my_ALC`+"-lsq-1"	        # Name for ChIMES lsq job
	default_keys[4 ] = "job_nodes"         ; default_values[4 ] = 	 "2"			  	        # Number of nodes for ChIMES lsq job
	default_keys[5 ] = "job_ppn"	       ; default_values[5 ] = 	 "36"  		  	       	        # Number of processors per node for ChIMES lsq job
	default_keys[6 ] = "job_walltime"      ; default_values[6 ] = 	 "1"			  	        # Walltime in hours for ChIMES lsq job
	default_keys[7 ] = "job_queue"         ; default_values[7 ] = 	 "pdebug"		  	        # Queue for ChIMES lsq job
	default_keys[8 ] = "job_account"       ; default_values[8 ] = 	 "pbronze"		  	        # Account for ChIMES lsq job
	default_keys[9 ] = "job_executable"    ; default_values[9 ] = 	 ""			  	        # Full path to executable for ChIMES lsq job
	default_keys[10] = "job_system"        ; default_values[10] = 	 "slurm"	       	       	        # slurm or torque 	 
	
	
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	################################
	# 1. Create the GEN_FF directory
	################################
	
	helpers.run_bash_cmnd("rm -rf GEN_FF")
	helpers.run_bash_cmnd("mkdir  GEN_FF")
	
	################################
	# 2. grab the fm.in and trajlist.in previous iteration, update the contents
	################################
	
	nfiles      = 0
	nframes_all = 0
	nframes_20  = 0	
	
	if my_ALC == 0:
		
		helpers.run_bash_cmnd("cp " + args["prev_gen_path"] + "/fm_setup.in"   + " GEN_FF/fm_setup.in")
		helpers.run_bash_cmnd("cp " + args["prev_gen_path"] + "/traj_list.dat" + " GEN_FF/traj_list.dat")
		helpers.run_bash_cmnd("cp " + ' '.join(glob.glob(args["prev_gen_path"] + "/*xyzf"  )) + " GEN_FF/")
		
		
		nfiles = int(helpers.head("GEN_FF/traj_list.dat",1)[0])
	else:

		helpers.run_bash_cmnd("cp " + ' '.join(glob.glob(args["prev_gen_path"] + "/*fm_setup.in"  )) + " GEN_FF/fm_setup.in")
		helpers.run_bash_cmnd("cp " + ' '.join(glob.glob(args["prev_gen_path"] + "/*traj_list.dat")) + " GEN_FF/traj_list.dat")

		# Get the number of files and number of frames in each file

		if args["prev_vasp_all_path"]:
	
			nfiles     += 1
			nframes_all = helpers.count_xyzframes_general(args["prev_vasp_all_path"] + "/OUTCAR.xyzf")

		if args["prev_vasp_20_path"]:
	
			nfiles    += 1	
			nframes_20 = helpers.count_xyzframes_general(args["prev_vasp_20_path"] + "/OUTCAR.xyzf")
			
		# Update the fm_setup.in file
	
	
		ifstream = open("GEN_FF/fm_setup.in",'r')
		runfile  = ifstream.readlines()
	
		ofstream = open("tmp",'w')	
	
		found1 = False
		found2 = False
	
		for i in xrange(len(runfile)):
	
			if found1:
				ofstream.write('\t' + "MULTI traj_list.dat" + '\n')
				found1 = False
				
			elif found2:
				ofstream.write('\t' + str(nframes_20+nframes_all) + '\n')
				found2 = False		
			else:
	
				ofstream.write(runfile[i])
	
				if "TRJFILE" in runfile[i]:
					found1 = True
					
				if "NFRAMES" in runfile[i]:
					found2 = True				
				
		helpers.run_bash_cmnd("mv tmp GEN_FF/fm_setup.in")				
		ofstream.close()
		ifstream.close()
		

		# Update the traj_list file
	
		ifstream = open("GEN_FF/traj_list.dat",'w')
		ifstream.write(`nfiles` + '\n')
	
		if args["prev_vasp_all_path"]:
			ifstream.write(`nframes_all` + " " + args["prev_vasp_all_path"] + "/OUTCAR.xyzf\n")
			
		if args["prev_vasp_20_path"]:
			ifstream.write(`nframes_20`  + " " + args["prev_vasp_20_path"]  + "/OUTCAR.xyzf\n")	
			
		ifstream.close()

	################################
	# 3. Set up and submit the .cmd file for the job
	################################
	
	# Create the task string
	
	job_task = "-n " + `int(args["job_nodes"])*int(args["job_ppn"])` + " " + args["job_executable"] + " fm_setup.in | tee fm_setup.log"
	
	if args["job_system"] == "slurm":
		job_task = "srun "   + job_task
	else:
		job_task = "mpirun " + job_task	

	# Launch the job
	
	os.chdir("GEN_FF")
	
	lsq_jobid_1 = helpers.create_and_launch_job(
		job_name       =      args["job_name"	 ] ,
		job_nodes      =  str(args["job_nodes"   ]),
		job_ppn        =  str(args["job_ppn"	 ]),
		job_walltime   =  str(args["job_walltime"]),
		job_queue      =      args["job_queue"   ] ,
		job_account    =      args["job_account" ] ,
		job_executable =      job_task,
		job_system     =      args["job_system"  ] ,
		job_file       = "run_chimeslsq.cmd")

	os.chdir("..")

	return lsq_jobid_1.split()[0]
	

def solve_amat(my_ALC, **kwargs):  


	""" 
	
	Generates parameters based on generated A- and b-matrices.
	
	Usage: solve_amat(1, <arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Requrires config.CHIMES_SOLVER.
	       Currently only supports lassolars and svd.
	       Returns a job_id for the submitted job.
	       Assumes last ALC's GEN_FF folder can be accessed from current 
	       ALC's base folder via ../ALC-(n-1)/GEN_FF.
	       
	WARNING: Only coded to work with forces and energies.
	
	WARNING: Weight generation function will need to be updated for stresses.
	       
	WARNING: Concatenation is untested and may be problematic for large A-matrices.
	       	
	"""

	# WARNING: Only coded to work with forces and energies ... weights writing function
	#          will need to be updated if stresses are included
	#
	# WARNING: Havent tested concatenation method ... it may be problematic for large A-matrices
	#

	# Notes:
	#
	# Assumes last ALC's GEN_FF folder can be accessed from current ALC's base folder via ../ALC-(n-1)/GEN_FF

	################################
	# 0. Set up an argument parser
	################################
	
	default_keys   = [""]*12
	default_values = [""]*12
	
	# Weights
	
	default_keys[0 ] = "weights_force"     ; default_values[0 ] = 	 "1.0" # Weights to be added to per-atom forces
	default_keys[1 ] = "weights_energy"    ; default_values[1 ] = 	 "5.0" # Weights to be added to per-frame and per cluster energies
	
	# LSQ controls
	
	default_keys[2 ] = "regression_alg"    ; default_values[2 ] = 	 "lassolars" # Regression algorithm to be used in lsq2
	default_keys[3 ] = "regression_var"    ; default_values[3 ] = 	 "1.0E-4"    # SVD eps or Lasso alpha
	
	# Overall job controls
	
	default_keys[4 ] = "job_name"	       ; default_values[4 ] =	 "ALC-"+ `my_ALC`+"-lsq-2"		# Name for ChIMES lsq job
	default_keys[5 ] = "job_nodes"         ; default_values[5 ] =	 "1"					# Number of nodes for ChIMES lsq job
	default_keys[6 ] = "job_ppn"	       ; default_values[6 ] =	 "36"					# Number of processors per node for ChIMES lsq job
	default_keys[7 ] = "job_walltime"      ; default_values[7 ] =	 "1"					# Walltime in hours for ChIMES lsq job
	default_keys[8 ] = "job_queue"         ; default_values[8 ] =	 "pdebug"				# Queue for ChIMES lsq job
	default_keys[9 ] = "job_account"       ; default_values[9 ] =	 "pbronze"				# Account for ChIMES lsq job
	default_keys[10] = "job_executable"    ; default_values[10] =	 ""					# Full path to executable for ChIMES lsq job
	default_keys[11] = "job_system"        ; default_values[11] =	 "slurm"				# slurm or torque	
	

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)

	################################
	# 1. Generate weights for the current ALC's trajectory
	################################
	
	os.chdir("GEN_FF")
	helpers.run_bash_cmnd("rm -f weights.dat")
	
	ifstream = open("b-labeled.txt",'r')
	weightfi = open("weights.dat",'w')
	
	contents = ifstream.readlines()
	
	for i in xrange(len(contents)):
	
		if contents[i].split()[0] == "+1":
			weightfi.write(str(args["weights_energy"])+'\n')
		else:
			weightfi.write(str(args["weights_force"])+'\n')
			
	ifstream.close()
	weightfi.close()
	
	os.chdir("..")

	
	################################
	# 2. Combine current weights, A.txt, and b.txt with previous ALC's
	#    ... Only needed for ALC >= 1
	################################
	
	print helpers.run_bash_cmnd("pwd")
	#print helpers.run_bash_cmnd(
	#print helpers.run_bash_cmnd(
	#print helpers.run_bash_cmnd(
	
	if my_ALC > 0: # "*_comb" files should always exist, because we create them for ALC-0 too
	
		# A-files
	
		prevfile    = "../ALC-" + `my_ALC-1` + "/GEN_FF/A_comb.txt"
	
		helpers.cat_specific("GEN_FF/A_comb.txt",       ["../ALC-" + `my_ALC-1` + "/GEN_FF/A_comb.txt",       "GEN_FF/A.txt"]      )
	
		helpers.cat_specific("GEN_FF/b_comb.txt",       ["../ALC-" + `my_ALC-1` + "/GEN_FF/b_comb.txt",       "GEN_FF/b.txt"]      )
	
		helpers.cat_specific("GEN_FF/weights_comb.dat", ["../ALC-" + `my_ALC-1` + "/GEN_FF/weights_comb.dat", "GEN_FF/weights.dat"])

		os.chdir("GEN_FF")
			
	else:
		os.chdir("GEN_FF")
		
		helpers.run_bash_cmnd("mv A.txt       A_comb.txt"      )
		helpers.run_bash_cmnd("mv b.txt       b_comb.txt"      ) 
		helpers.run_bash_cmnd("mv weights.dat weights_comb.dat")
			
	
	# Sanity checks
	
	print "A-mat entries:  ",helpers.run_bash_cmnd("wc -l A_comb.txt").split()[0]
	print "b-mat entries:  ",helpers.run_bash_cmnd("wc -l b_comb.txt").split()[0]
	print "weight entries: ",helpers.run_bash_cmnd("wc -l weights_comb.dat").split()[0]
	

	################################
	# 3. Run the actual fit
	################################
	
	# Create the task string
			
	job_task = args["job_executable"] + " --A A_comb.txt --b b_comb.txt --weights weights_comb.dat --algorithm " + args["regression_alg"]

	if "svd" in args["regression_alg"]:
		job_task += " --eps " + args["regression_var"]
	elif "lasso" in args["regression_alg"]:
		job_task += " --alpha " + args["regression_var"]
	else:
		print "ERROR: unknown regression algorithm: ", args["regression_alg"]
		exit()

	# Launch the job
	 
	job_task += " | tee params.txt "

	run_py_jobid = helpers.create_and_launch_job(
		job_name       =     args["job_name"	] ,
		job_nodes      = str(args["job_nodes"	]),
		job_ppn        = str(args["job_ppn"	]),
		job_walltime   = str(args["job_walltime"]),
		job_queue      =     args["job_queue"	] ,
		job_account    =     args["job_account" ] ,
		job_executable =     job_task,
		job_system     =     args["job_system"  ] ,
		job_file       =     "run_lsqpy.cmd")

	os.chdir("..")
	
	return run_py_jobid.split()[0]




