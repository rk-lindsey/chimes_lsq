"""

	A configuration file containing all user-specified parameters.
	
	Notes: Relative paths will almost certainly break the code.
	
"""


################################
# Set user specified variables 
################################

#WARNING: This code probably won't work if you use relative paths!

##### General variables

EMAIL_ADD     = "lindsey11@llnl.gov" # driver will send updates on the status of the current run ... If blank (""), no emails are sent

SEED          = 1

ATOM_TYPES    = ["C", "O"]
NO_CASES      = 3
THIS_SMEAR    = 2400       # All self-consistent vasp calcs will use this smearing (K) 

# FYI:
#
# CASE-0: 2400 K; 1.79 gcc
# CASE-1: 6500 K; 2.50 gcc 
# CASE-2: 9350 K; 2.56 gcc



DRIVER_DIR    = "/p/lustre1/rlindsey/RC4B_RAG/11-12-18/PYTHON_DRIVER/"			# This driver's src location
WORKING_DIR   = DRIVER_DIR + "/tests/CASES-3_SOLVER-DLASSO_MATRIX-SINGLE_FIT-F+E/run_test/"	# Directory from which all ALCs will be run
CHIMES_SRCDIR = "/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/SVN_SRC-11-12-18/"		# Location of ChIMES src files ... expects also post_proc_lsq2.py there


##### General HPC options

HPC_PPN       = 36
HPC_ACCOUNT   = "pbronze"
HPC_SYSTEM    = "slurm"
HPC_PYTHON    = "/usr/tce/bin/python"
HPC_EMAIL     = False # Adds "MSUB -m abe" to slurm scripts


##### ChIMES LSQ

ALC0_FILES    = WORKING_DIR + "ALL_BASE_FILES/ALC-0_BASEFILES/"

CHIMES_LSQ    = CHIMES_SRCDIR + "chimes_lsq"
CHIMES_SOLVER = CHIMES_SRCDIR + "lsq2.py"
CHIMES_POSTPRC= CHIMES_SRCDIR + "post_proc_lsq2.py"

WEIGHTS_FORCE = 1.0 
WEIGHTS_ENER  = 5.0

REGRESS_ALG   = "dlasso" # "lassolars"
REGRESS_VAR   = "1.0E-4"

CHIMES_BUILD_NODES = 2
CHIMES_BUILD_QUEUE = "pdebug"
CHIMES_BUILD_TIME  = 1		# Hours

CHIMES_SOLVE_NODES = 1
CHIMES_SOLVE_QUEUE = "pdebug"
CHIMES_SOLVE_TIME  = 1		# Hours


##### ChIMES MD

CHIMES_MD     = CHIMES_SRCDIR + "chimes_md"
CHIMES_MOLANAL= "/g/g17/rlindsey/CURR_TRACKED-GEN/contrib/molanlal/"  
CHIMES_MDFILES= WORKING_DIR + "ALL_BASE_FILES/CHIMESMD_BASEFILES/"

CHIMES_PEN_PREFAC = 1.0E6
CHIMES_PEN_DIST   = 0.02

CHIMES_MD_NODES = 8
CHIMES_MD_QUEUE = "pbatch"
CHIMES_MD_TIME  = 1		# Hours


##### Cluster specific paths/variables

DO_CLUSTER    = True
TIGHT_CRIT    = WORKING_DIR + "ALL_BASE_FILES/tight_bond_crit.dat"
LOOSE_CRIT    = WORKING_DIR + "ALL_BASE_FILES/loose_bond_crit.dat"
CLU_CODE      = DRIVER_DIR  + "/utilities/new_ts_clu.cpp"


##### ALC-specific variables ... Note: Hardwired for partial memory mode

MEM_BINS = 40
MEM_CYCL = MEM_BINS/10
MEM_NSEL = 10
MEM_ECUT = 100.0

CALC_REPO_ENER_QUEUE = "pbatch"
CALC_REPO_ENER_TIME  = 4


##### VASP Specific variables

VASP_FILES   = WORKING_DIR + "ALL_BASE_FILES/VASP_BASEFILES"
VASP_POSTPRC = CHIMES_SRCDIR + "vasp2xyzf.py"
VASP_NODES   = 6
VASP_TIME    = 1	# Hours
VASP_QUEUE   = "pbatch"
VASP_EXE     = "/usr/gapps/emc-vasp/vasp.5.4.1/build/std/vasp"
VASP_MODULES = "mkl" # Currently unused, but should be included



