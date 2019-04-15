from subprocess import check_output 
from subprocess import CalledProcessError
import time
import glob
import sys
import math

def run_bash_cmnd(cmnd_str):

	# Notes:
	#
	# Runs a shell command (expects bash), and captures any returned output

	msg = ""
	
	try:
		msg = check_output(cmnd_str.split())
	except CalledProcessError as err_msg:
		msg = err_msg

	return msg
	
def run_bash_cmnd_to_file(outfile, cmnd_str):
	
	ofstream = open(outfile,'r')
	ofstream .write(run_bash_cmnd(cmnd_str))
	ofstream .close()
	
	
def run_bash_cmnd_suppress(cmnd_str):

	# Notes:
	#
	# Runs a shell command (expects bash), ... output is ignored

	call(cmnd_str.split())	
	return
	
def cat_to_var(*argv):

	files_to_cat = argv # This is a pointer!
	
	contents = []
	
	for f in files_to_cat:
		
		ifstream = open(f,'r')
			
		contents += ifstream.readlines()	
			
		ifstream.close()
		
	return contents
	
	
def cat_specific(outfilename, *argv):

	files_to_cat = argv[0] # This is a pointer!

	with open(outfilename, "w") as ofstream:
	
		for f in files_to_cat:
			with open(f, "r") as ifstream:
				ofstream.write(ifstream.read())

def cat_pattern(outfilename, pattern):

	files_to_cat = glob.glob(pattern)
	
	with open(outfilename, "wb") as ofstream:
		for f in files_to_cat:
			with open(f, "rb") as ifstream:
				ofstream.write(ifstream.read())
				
def head(*argv): # Syntax is: head(file_to_head) or head(file_to_head,nlines_requested)
	
	nlines = 10
	
	if len(argv) == 2:
		nlines = int(argv[1])
	elif len(argv) > 2:
		print "ERROR: Unrecognized head command: ", argv
		exit()
		
	ifstream = open(argv[0],'r')
	
	contents = []
	
	for i in xrange(nlines):
		contents.append(ifstream.readline())
	
	ifstream.close()
	return contents
	
def wc_l(infile):

	nlines = 0
	
	with open(infile, "r") as ifstream:
		for line in ifstream:
			nlines += 1
	return nlines		
	
def count_xyzframes_general(infile):

	nframes = 0
	
	with open(infile, "r") as ifstream:
		for line in ifstream:
			if len(line.split()) == 1:
				nframes += 1
	return nframes	
	
def create_and_launch_job(*argv, **kwargs):

	# Notes:
	#
	# "job_executable" is the string to be executed in the submit script
	# if "job_executable" is empty, will use whatever list of commands are specified in *argv

	################################
	# 0. Set up an argument parser
	################################
	
	default_keys   = [""]*9
	default_values = [""]*9

	# Overall job controls
	
	default_keys[0 ] = "job_name"	       ; default_values[0 ] =	 "ALC-x-lsq-1"		# Name for ChIMES lsq job
	default_keys[1 ] = "job_nodes"         ; default_values[1 ] =	 "2"			# Number of nodes for ChIMES lsq job
	default_keys[2 ] = "job_ppn"	       ; default_values[2 ] =	 "36"			# Number of processors per node for ChIMES lsq job
	default_keys[3 ] = "job_walltime"      ; default_values[3 ] =	 "1"			# Walltime in hours for ChIMES lsq job
	default_keys[4 ] = "job_queue"         ; default_values[4 ] =	 "pdebug"		# Queue for ChIMES lsq job
	default_keys[5 ] = "job_account"       ; default_values[5 ] =	 "pbronze"		# Account for ChIMES lsq job
	default_keys[6 ] = "job_executable"    ; default_values[6 ] =	 ""			# Full path to executable for ChIMES lsq job
	default_keys[7 ] = "job_system"        ; default_values[7 ] =	 "slurm"		# slurm or torque	
	default_keys[8 ] = "job_file"          ; default_values[8 ] =	 "run.cmd"		# Name of the resulting submit script	
	

	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	################################
	# 1. Create the job file
	################################	
	
	
	run_bash_cmnd("rm -f " + args["job_file"])
	
	JOB = []
	JOB.append(" -N " + args["job_name"])
	JOB.append(" -l " + "nodes=" + args["job_nodes"] + ":ppn=" + args["job_ppn"])
	JOB.append(" -l " + "walltime=" + args["job_walltime"] + ":00:00")	  	   
	JOB.append(" -q " + args["job_queue"])
	JOB.append(" -m " + "abe")   
	JOB.append(" -A " + args["job_account"])  
	JOB.append(" -V " )
	JOB.append(" -o " + "stdoutmsg")
	
	ofstream = open(args["job_file"],'w')
	ofstream.write("#!/bin/bash\n")
	
	for i in xrange(len(JOB)):
	
		if args["job_system"] == "slurm":
			JOB[i] = "#MSUB" + JOB[i]
		elif args["job_system"] == "torque":
			JOB[i] = "#PBS"  + JOB[i]
		else:
			print "ERROR: Unknown job_system: ", args["job_system"]
			exit()
			
		ofstream.write(JOB[i] + '\n')
	
	if args["job_executable"]:
	
		ofstream.write(args["job_executable"] + '\n')
	else:
		job_list = argv[0]
		
		for i in xrange(len(job_list)):
			
			ofstream.write(job_list[i] + '\n')
	
	ofstream.close()
	
	################################
	# 2. Launch the job file
	################################

	jobid = None
	
	if args["job_system"] == "slurm":
		jobid = run_bash_cmnd("msub " + args["job_file"])
	else:	
		jobid = run_bash_cmnd("qsub " + args["job_file"])

	return jobid	
	

def wait_for_job(active_job, **kwargs):

	# Notes:
	#
	# Accepts a jobid and queries the queueing system to determine
	# whether the job is still active ... does not return until job is complete
	# will likely need to be modified for ~parallel learning~

	################################
	# 0. Set up an argument parser
	################################

	default_keys   = [""]*3
	default_values = [""]*3
	
	default_keys  [0] = "job_system" ; default_values[0] = "slurm"
	default_keys  [1] = "verbose"    ; default_values[1] = False 
	default_keys  [2] = "job_name"   ; default_values[2] = "unspecified" 
	
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	active_job = str(active_job).split()[0]
	
	
	################################
	# 1. Determine job status, hold until complete
	################################

	while True:
		
		check_job = ""
		
		if args["job_system"] == "slurm":
			check_job = "squeue -j " + active_job
			
		elif args["job_system"] == "torque":
			print "ERROR: torque support not yet implemented in wait_for_job"
		else:
			print "ERROR: Unknown job_system: ", args["job_system"]
			exit()
			
		if active_job in  run_bash_cmnd(check_job):
		
			if args["verbose"]:
				print "Sleeping for 60 more seconds while waiting for job ", active_job, "...", args["job_name"]
		
			time.sleep(60) # sleep for 60 seconds
		else:		
			print "Breaking ... "
			break				
			
	return
	

def wait_for_jobs(*argv, **kwargs):

	# Notes:
	#
	# Accepts a jobid and queries the queueing system to determine
	# whether the job is still active ... does not return until job is complete
	# will likely need to be modified for ~parallel learning~

	################################
	# 0. Set up an argument parser
	################################

	default_keys   = [""]*3
	default_values = [""]*3
	
	default_keys  [0] = "job_system" ; default_values[0] = "slurm"
	default_keys  [1] = "verbose"    ; default_values[1] = False 
	default_keys  [2] = "job_name"   ; default_values[2] = "unspecified" 
	
	args = dict(zip(default_keys, default_values))
	args.update(kwargs)
	
	active_jobs = argv[0] # Pointer!
	
	
	################################
	# 1. Determine job status, hold until complete
	################################

	njobs  = len(active_jobs)
	
	active = [True]*njobs
	
	while True:
	
		for i in xrange(njobs):
			check_job = ""
			
			if type(active_jobs[i]) == type(1):
				active_jobs[i] = str(active_jobs[i])
		
			if args["job_system"] == "slurm":
				check_job = "squeue -j " + active_jobs[i]
			
			elif args["job_system"] == "torque":
				print "ERROR: torque support not yet implemented in wait_for_job"
			else:
				print "ERROR: Unknown job_system: ", args["job_system"]
				exit()
			
			if active_jobs[i] in  run_bash_cmnd(check_job):
				active[i] = True
			else:
				active[i] = False
			
		
		if True in active:
		
			if args["verbose"]:
				print "Sleeping for 60 more seconds while waiting for jobs ", active_jobs, "...", args["job_name"]
		
			time.sleep(60) # sleep for 60 seconds
		else:		
			print "Breaking ... "
			break					
	return
	

def str2bool(v):
    return v.lower() in ("true")

def break_apart_xyz(*argv):

	# Takes as input a number of frames and a .xyz or .xyzf file, and breaks it apart 
	# into frames. 

	# Optional: break into chunks of (3rd arg)


	print "Breaking apart file: ", argv[1]
	print "WARNING: Converting forces from Hartree/bohr to simulation units (kca/mol/Ang)"
	

	# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
	FRAMES = int(argv[0])

	#########

	# What is the input .xyz file?
	XYZFILE = open(argv[1],"r")

	CHUNK_LEN = 1
	if len(argv) >= 3:
	    CHUNK_LEN = int(argv[2])
	    
	FIRST_ONLY = False
	if len(argv) >= 4:
	    FIRST_ONLY = argv[3]

	#########

	ZEROES = len(str(FRAMES))+1


	for f in xrange(FRAMES):

	    if f%CHUNK_LEN == 0:
	    
	        if f > 1:
	            OFSTREAM.close()
	            FRSTREAM.close()

	        # Generate the output filename

	        TAG = ""
	        for i in xrange(ZEROES):
	            if f == 0:
	                if f+1 < pow(10.0,i):
	                    for j in xrange(ZEROES-i):
	                        TAG += "0"
	                    TAG += `f`
	                    break
	                
	            elif f < pow(10.0,i):
	                for j in xrange(ZEROES-i):
	                    TAG += "0"
	                TAG += `f`
	                break

	        OUTFILE  = argv[1]
	        FORCES   = argv[1]
	        TESTER   = OUTFILE [0:-4]
	        TESTER   = TESTER  [-1]

	        if TESTER == ".":
	            FORCES   = OUTFILE[0:-5] + "_FORCES_#" + TAG + ".xyzf" 
	            OUTFILE  = OUTFILE[0:-5] + "_#"        + TAG + ".xyzf" 
	            
	        else:
	            FORCES   = OUTFILE[0:-4] + "_FORCES_#" + TAG + ".xyz" 
	            OUTFILE  = OUTFILE[0:-4] + "_#"        + TAG + ".xyz" 
	            

	        OFSTREAM = open(OUTFILE,"w")
	        FRSTREAM = open(FORCES,"w")
	        
	    if FIRST_ONLY and f%CHUNK_LEN > 0: # We still need to filter through ignored frames
	    
	        # Read the first line to get the number of atoms in the frame,
	        # print back out to the xyzf file
	    
	        ATOMS = XYZFILE.readline()
	    
	        ATOMS = ATOMS.split()
		
	        ATOMS = int(ATOMS[0])
	    
	        # Read/print the comment line

	        XYZFILE.readline()

	        # Now, read/print each atom line
	    
	        for j in xrange(ATOMS):
	            XYZFILE.readline()
	            
	    else:
	    
	        # Read the first line to get the number of atoms in the frame,
	        # print back out to the xyzf file
	    
	        ATOMS = XYZFILE.readline()
	        
	        OFSTREAM.write(ATOMS)
	    
	        ATOMS = ATOMS.split()
	        
	        ATOMS = int(ATOMS[0])
	    
	        # Read/print the comment line

	        OFSTREAM.write( XYZFILE.readline())

	        # Now, read/print each atom line
	    
	        for j in xrange(ATOMS):
	    
	            LINE = XYZFILE.readline()

	            OFSTREAM.write(LINE)
	        
	            LINE = LINE.split()

	            if len(LINE)>4:
	                FRSTREAM.write(`float(LINE[4])*(627.50960803*1.889725989)` + '\n')
	                FRSTREAM.write(`float(LINE[5])*(627.50960803*1.889725989)` + '\n')
	                FRSTREAM.write(`float(LINE[6])*(627.50960803*1.889725989)` + '\n')    
	return


	
import sys


def dftbgen_to_xyz(*argv):

	#NOTE: Assumes an orthorhombic box

	# How many frames are there? ... just do grep -F "Step" <file> | wc -l to find out
	FRAMES = int(argv[0])

	# What is the input file?
	IFSTREAM = open(argv[1],"r")

	SKIP = 1
	if len(argv) == 3:
		SKIP = int(argv[2])

	# What is the outputfile
	OUTFILE  = argv[1]
	OUTFILE  = OUTFILE[0:-4] + ".xyz" # replace ".gen" with ".xyz"
	OFSTREAM = open(OUTFILE,"w")

	BOXFILE  = argv[1]
	BOXFILE  = BOXFILE[0:-4] + ".box" # replace ".gen" with ".xyz"
	BOXSTREAM = open(BOXFILE,"w")


	for i in xrange(FRAMES):
		
		# Read the first line to get the number of atoms in the frame
		
		ATOMS = IFSTREAM.readline()
		ATOMS = ATOMS.split()
		ATOMS = int(ATOMS[0])
		
		# Read the next line to get the atom types
		
		SYMBOLS = IFSTREAM.readline()
		SYMBOLS = SYMBOLS.split()
		
		# Print the header bits of the xyz file
		
		if (i+1)%SKIP == 0:
		
			OFSTREAM.write(`ATOMS` + '\n')
			OFSTREAM.write("Frame " + `i+1` + '\n')
		
		# Now read/print all the atom lines in the present frame
		
		for j in xrange(ATOMS):
		
			LINE = IFSTREAM.readline()
			LINE = LINE.split()
			
			# Replace the atom type index with a chemical symbol
		
			for k in xrange(len(SYMBOLS)):
				if k+1 == int(LINE[1]):
					LINE[1] = SYMBOLS[k]
					break
					
			# Print out the line
			
			if (i+1)%SKIP == 0:
			
				OFSTREAM.write(' '.join(LINE[1:len(LINE)]) + '\n')
			
		# Finally, read the box lengths... assume cubic
		
		LINE = IFSTREAM.readline()	# Cell angles?
		
		LINE = IFSTREAM.readline().split()	
		X = LINE[0]
		
		LINE = IFSTREAM.readline().split()	
		Y = LINE[1]
		
		LINE = IFSTREAM.readline().split()	
		Z = LINE[2]	
		
		if (i+1)%SKIP == 0:
		
			BOXSTREAM.write(X + " " + Y + " " + Z + '\n')
			
	return

