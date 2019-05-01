# Global (python) modules

import glob
import os

# Localmodules

import helpers

def cleanup_and_setup(*argv, **kwargs):

	""" 
	
	Removes VASP-X folders, if they exist in the build_dir.
	
	Usage: cleanup_and_setup(<arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Expects to be run from ALC-X folder
	       Should only be run once per ALC
	       	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	### ...argv
	
	args_targets   = argv[0] # This is a pointer!
	args_this_case = -1
	
	if len(argv) == 2:
		args_this_case = argv[1]

	### ...kwargs
	
	default_keys   = [""]*1
	default_values = [""]*1

	# VASP specific controls

	default_keys[0 ] = "build_dir" 	   ; default_values[0 ] = "." 				# Post_proc_lsq*py file... should also include the python command

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)	
	
	print helpers.run_bash_cmnd("pwd")
	
	for i in args_targets: # 20 all
	
		if args_this_case == -1:
	
			for i in args_targets: # 20 all

				vasp_dir = args["build_dir"] + "/VASP-" + i
			
				if os.path.isdir(vasp_dir):
			
					helpers.run_bash_cmnd("rm -rf " + vasp_dir)
					helpers.run_bash_cmnd("mkdir  " + vasp_dir)
								
				else:
					helpers.run_bash_cmnd("mkdir  " + vasp_dir)
				
			return

	for i in args_targets: # 20 all

		vasp_dir = args["build_dir"] + "/VASP-" + i
		
		if os.path.isdir(vasp_dir):
		
			helpers.run_bash_cmnd("rm -f " + vasp_dir + "/*")

			vasp_dir = args["build_dir"] + "/VASP-" + i + "/CASE-" + `args_this_case`

			if os.path.isdir(vasp_dir):
				helpers.run_bash_cmnd("rm -f " + vasp_dir + "/*")
			else:
				helpers.run_bash_cmnd("mkdir  -p " + vasp_dir)
			
		else:
			helpers.run_bash_cmnd("mkdir -p  " + vasp_dir + "/CASE-" + `args_this_case`)


def continue_job(*argv, **kwargs):

	""" 
	
	Checks whether all VASP ingle point calculations ran, resubmits if needed.
	
	Usage: continue_job(<arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Returns a SLURM jobid list
	       	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	### ...argv
	
	args_targets = argv[0] # This is a pointer!
	args_case    = argv[1]
	
	### ...kwargs
	
	default_keys   = [""]*1
	default_values = [""]*1
	
	default_keys[0 ] = "job_system"    ; default_values[0 ] = "slurm"		      # slurm or torque       
	
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)	
	
	
	################################
	# 1. Check on job status, resubmit if needed
	################################	
	
	job_list = []
	
	for i in xrange(len(args_targets)): # 20 all
	
		if os.path.isdir("VASP-" + args_targets[i] + "/CASE-" + `args_case`):
		
			os.chdir("VASP-" + args_targets[i] + "/CASE-" + `args_case`)

			# Count the number of possible jobs
			
			count_POSCAR = len(helpers.run_bash_cmnd("ls " + ' '.join(glob.glob("*.POSCAR"))).split())
			
			# Count the number of completed jobs
			
			count_OUTCAR = len(helpers.run_bash_cmnd("ls " + ' '.join(glob.glob("*.OUTCAR"))).split())

			# If all jobs haven't completed, resubmit
			
			if count_POSCAR > count_OUTCAR:

				if args["job_system"] == "slurm":
					job_list.append(helpers.run_bash_cmnd("msub run_vasp.cmd").replace('\n', ''))
				else:	
					job_list.append(helpers.run_bash_cmnd("msub run_vasp.cmd").replace('\n', ''))
					
			os.chdir("../..")

	return job_list	


def generate_POSCAR(inxyz, *argv):

	""" 
	
	Converts a .xyz file to a VASP POSCAR file.
	
	Usage: continue_job(my_file.xyz, ["C","O"], [1,2])
	
	Notes: Final two arguments should be two lists:
	              the first giving possible atom types
		      the second giving number of atoms of each type 
		      
	WARNING: Assumes an orthorhombic box
	       	
	"""

	atm_types = argv[0] # pointer!
	smearing  = int(argv[1])
	ifstream = open(inxyz,'r')
	ofstream = open(inxyz + ".POSCAR", 'w')
	
	# Set up the header portion
	
	ofstream.write(' '.join(atm_types) + " ( " + `smearing` + " eV)" + '\n')
	ofstream.write("1.0" + '\n')
	
	box = ifstream.readline()		# No. atoms
	box = ifstream.readline().split()	# Box lengths
	
	ofstream.write(box[0] + " 0.000       0.000" + '\n')
	ofstream.write("0.000 " + box[1]  + " 0.000" + '\n')
	ofstream.write("0.000     0.000 " +   box[2] + '\n')
	
	
	# Count up the number of atoms of each type
	
	contents = ifstream.readlines()
	contents.sort() # Vasp expects sorted coordinates
	
	natm_types = [0]*len(atm_types)
	
	for i in xrange(len(contents)):
	
		for j in xrange(len(atm_types)):
			
			if atm_types[j] in contents[i]:
				
				natm_types[j] += 1
				
	
	# Remove any elements that don't appear in the coordinate file
	# But first, make sure everything is sorted 
	
	tmp = zip(atm_types,natm_types)
	tmp.sort() # Sorts by the first entry (atom types) ... alphabetical
	
	for i in xrange(len(tmp)-1,-1,-1): # Goes through list backwards
		
		if tmp[i][1] == 0: 
		
			tmp.pop(i)
	
	atm_types  = list(zip(*tmp)[0])
	natm_types = map(str, list(zip(*tmp)[1]))
	
	ofstream.write( ' '.join(atm_types ) +'\n')
	ofstream.write( ' '.join(natm_types) +'\n')
	ofstream.write("Cartesian" + '\n')
	
	for i in xrange(len(contents)):
		
		line = contents[i].split()
		line = ' '.join(line[1:]) + '\n'
		ofstream.write(line)
		
	ofstream.close()


def post_process(*argv, **kwargs):

	""" 
	
	Converts a converts a VASP OUTCAR file to .xyzf file
	
	Usage: post_process(<arguements>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	### ...argv
	
	args_targets    = argv[0] # ... all ... 20
	args_properties = argv[1] # "ENERGY STRESS" ...etc for the post process script
	args_cases      = argv[2]
	
	
	### ...kwargs
	
	default_keys   = [""]*1
	default_values = [""]*1


	# VASP specific controls
	
	default_keys[0 ] = "vasp_postproc"  ; default_values[0 ] = ""		# POTCAR, KPOINTS, and INCAR	

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)	
	
	################################

	for i in xrange(len(args_targets)): # 20 all

		if not os.path.isdir("VASP-" + args_targets[i]):

			continue
			
		print "Working on:","VASP-" + args_targets[i]
	
		os.chdir("VASP-" + args_targets[i])
	
		helpers.run_bash_cmnd("rm -f OUTCAR.xyzf")
		
		outcar_list = []
		
		for j in xrange(args_cases):
			
			outcar_list += glob.glob("CASE-" + `j` + "/*.OUTCAR")
			
		for j in xrange(len(outcar_list)):

			print helpers.run_bash_cmnd(args["vasp_postproc"] + " " + outcar_list[j] + " 1 " + args_properties + " | grep ERROR ")
			
			if os.path.isfile(outcar_list[j] + ".xyzf"):
			
				if os.path.isfile("OUTCAR.xyzf"):
			
					helpers.cat_specific("tmp.dat", ["OUTCAR.xyzf", outcar_list[j] + ".xyzf"])
				else:
					helpers.cat_specific("tmp.dat", [outcar_list[j] + ".xyzf"])
				
				helpers.run_bash_cmnd("mv tmp.dat OUTCAR.xyzf")
		
		os.chdir("..")

def setup_vasp(my_ALC, *argv, **kwargs):

	""" 
	
	Sets up and launches VASP single point calculations
	
	Usage: setup_vasp(1, <arguments>)
	
	Notes: See function definition in helpers.py for a full list of options. 
	       Expects to be run from teh ALC-X folder.
	       Expects the "all" file in the ALC-X folder.
	       Expects the "20" file in the ALC-X/INDEP_X folder.
	       Requries a list of atom types.
	       Expects a potcar file for each possible atom type, named like:
	              X.POSCAR
	       All other input files (INCAR, KPOINTS, etc) are taken from config.VASP_FILES.
	       Returns a SLURM jobid
	       
	WARNING: For ALC-0, smearing is taken to be the corresponding temeprature
	         in the traj_list.dat file.
	
	WARNING: For ALC-1+, smearing is a specified parameter.	
	
	WARNING: Largely only intended for liquids, so far.	
	      	
	"""
	
	################################
	# 0. Set up an argument parser
	################################
	
	### ...argv
	
	args_targets = argv[0] # This is a pointer!
	atm_types    = argv[1]
	my_case      = ""
	
	if len(argv) >= 3:
		my_case = str(argv[2])
		
	if len(argv) >= 4:
		my_smear = str(argv[3])
	else:
		my_smear = 0
		
	
	### ...kwargs
	
	default_keys   = [""]*12
	default_values = [""]*12


	# VASP specific controls
	
	default_keys[0 ] = "basefile_dir"  ; default_values[0 ] = "../VASP_BASEFILES/"		# POTCAR, KPOINTS, and INCAR
	default_keys[1 ] = "traj_list" 	   ; default_values[1 ] = "" 				# Post_proc_lsq*py file... should also include the python command
	default_keys[2 ] = "modules" 	   ; default_values[2 ] = "mkl" 			# Post_proc_lsq*py file... should also include the python command
	default_keys[3 ] = "build_dir" 	   ; default_values[3 ] = "." 				# Post_proc_lsq*py file... should also include the python command


	# Overall job controls	
	
	default_keys[4 ] = "job_nodes"     ; default_values[4 ] = "2"			      # Number of nodes for ChIMES md job
	default_keys[5 ] = "job_ppn"	   ; default_values[5 ] = "36"  		      # Number of processors per node for ChIMES md job
	default_keys[6 ] = "job_walltime"  ; default_values[6 ] = "1"			      # Walltime in hours for ChIMES md job
	default_keys[7 ] = "job_queue"     ; default_values[7 ] = "pdebug"		      # Queue for ChIMES md job
	default_keys[8 ] = "job_account"   ; default_values[8 ] = "pbronze"		      # Account for ChIMES md job
	default_keys[9 ] = "job_executable"; default_values[9 ] = ""			      # Full path to executable for ChIMES md job
	default_keys[10] = "job_system"    ; default_values[10] = "slurm"		      # slurm or torque       
	default_keys[11] = "job_file"	   ; default_values[11] = "run.cmd"		      # Name of the resulting submit script   
	

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)	
	
	
	run_vasp_jobid = []
	
	################################
	# 1. Set up and launch the vasp single point calculations
	################################
		
	for i in xrange(len(args_targets)): # 20 all
	
		curr_dir = helpers.run_bash_cmnd("pwd").rstrip() # This should be either CASE-X... (ALC>=1) or ...? (ALC==0)
	
		vasp_dir = "VASP-" + args_targets[i] + "/CASE-" + my_case
		
		helpers.run_bash_cmnd("mkdir -p " + vasp_dir)		
			
		# Set up an launch the job
		
		os.chdir(vasp_dir)
		
		if args_targets[i] == "20":
		
			my_md_path = curr_dir + "/CASE-" + str(my_case) + "_INDEP_0/" # ONLY EVER COMES FROM FIRST INDEP " + str(my_indep) + "/"

			my_file = "case_" + str(my_case) + ".indep_0.traj_" + args_targets[i] + "F.xyz"

			helpers.run_bash_cmnd("cp " + my_md_path + "/traj_20F.xyz " + my_file)
			
			atoms  = helpers.head(my_file,1)[0].rstrip()
			frames = helpers.wc_l(my_file  )
			frames = frames / (2+int(atoms))

			helpers.break_apart_xyz(frames, my_file)
			
			helpers.run_bash_cmnd("rm -f " + ' '.join(glob.glob("*FORCES*")))
			

			# Generate the POSCAR files

			target_files = ' '.join(glob.glob("case_" + str(my_case) + ".indep_0.traj_" + args_targets[i] + "F_#*")).split()

			for j in target_files:

				generate_POSCAR(j, atm_types, my_smear)	
				
		else:

			ifstream     = open(curr_dir + "/" + args_targets[i] + ".selection.dat", 'r') # all.selection.dat ... a list of indices
			target_files = ifstream.readlines()
			ifstream     .close()
			
			ifstream     = open(curr_dir + "/" + args_targets[i] + ".xyzlist.dat",   'r') # all.xyzlist.dat ... a list of file names
			contents     = ifstream.readlines()
			ifstream     .close()
			
			temp = my_smear
			
			if args["traj_list"]:
			
				# Need to figure out the possible temperatures ... 
			
				ifstream = open(args["traj_list"], 'r')
			
				temps    = ifstream.readlines()[1:]
			
				for j in xrange(len(temps)):
					temps[j] = temps[j].split()
	
				ifstream.close()
			
			for j in xrange(len(target_files)):
			
				sel_file = contents[int(target_files[j].rstrip())].split()[3]
				
				case_check = "CASE-" + my_case
				
				if case_check not in sel_file:
					continue

				if args["traj_list"]:
				
					for k in xrange(len(temps)):
					
						if temps[k][1] in sel_file:
					
							temp = temps[k][2]
							break
	
				generate_POSCAR(curr_dir + "/" + sel_file, atm_types, temp)
				
				helpers.run_bash_cmnd("mv " + curr_dir + "/" + sel_file + ".POSCAR tmp.POSCAR")
				
	
				sel_POSCAR = sel_file + ".POSCAR"
				sel_POSCAR = sel_POSCAR.replace('/','.')
				
				helpers.run_bash_cmnd("mv  tmp.POSCAR " + sel_POSCAR)


		################################
		# 2. Launch the actual job 
		################################
	
		# Grab the necessary files
	
		helpers.run_bash_cmnd("cp " + ' '.join(glob.glob(args["basefile_dir"] + "/*")) + " .")
	
		# Create the task string
				
		job_task = []
		job_task.append("module load " + args["modules"] + '\n')
	
	
		for i in xrange(len(atm_types)):
	
			job_task.append("ATOMS[" + `i` + "]=" + atm_types[i] + '\n')
	
		job_task.append("for j in $(ls *.POSCAR)	")	
		job_task.append("do				")
		job_task.append("	TAG=${j%*.POSCAR}	")	
		job_task.append("	cp ${TAG}.POSCAR POSCAR	")
		job_task.append("	CHECK=${TAG}.OUTCAR	")
		job_task.append("	if [ -e ${CHECK} ] ; then ")	
		job_task.append("		continue	")	
		job_task.append("	fi			")
		job_task.append("	TEMP=`awk '{print $4; exit}' POSCAR`")
		job_task.append("	cp ${TEMP}.INCAR INCAR		")	
		job_task.append("	rm -f POTCAR		")
		job_task.append("	for k in ${ATOMS[@]}	")
		job_task.append("	do			")
		job_task.append("		NA=`awk -v atm=\"$k\" \'{if(NR==6){for(i=1;i<=NF;i++){ if($i==atm){getline;print $i;exit}} print \"0\"}}\' POSCAR` ")
		job_task.append("		if [ $NA -gt 0 ] ; then ")
		job_task.append("			cat ${ATOMS[$k]}.POTCAR >> POTCAR ")
		job_task.append("		fi ")	
		job_task.append("	done	")	
		job_task.append("	srun -N " + `args["job_nodes" ]` + " -n " + `int(args["job_nodes"])*int(args["job_ppn"])` + " " + args["vasp_exe"] + " > ${TAG}.out  ")
		job_task.append("	cp OUTCAR  ${TAG}.OUTCAR	")
		job_task.append("	cp OSZICAR ${TAG}.OSZICAR	")
		job_task.append("	rm -f OUTCAR CHG DOSCAR XDATCAR POSCAR CHGCAR EIGENVAL PCDAT XDATCAR CONTCAR IBZKPT OSZICAR WAVECAR  ")	
		job_task.append("done	")
	
		this_jobid = helpers.create_and_launch_job(job_task,
			job_name       =          "vasp_spcalcs"  ,
			job_nodes      = str(args["job_nodes"	]),
			job_ppn        = str(args["job_ppn"	]),
			job_walltime   = str(args["job_walltime"]),
			job_queue      =     args["job_queue"	] ,
			job_account    =     args["job_account" ] ,
			job_system     =     args["job_system"  ] ,
			job_file       =     "run_vasp.cmd")
			
		run_vasp_jobid.append(this_jobid.split()[0])	
	
		os.chdir(curr_dir)
	
	return run_vasp_jobid		
