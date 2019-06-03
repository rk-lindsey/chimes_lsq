#!/bin/bash	

# Run with: 
#
#	NOTE: 	Expects run_md.base in working directory	
#
#	1.
#		msub -l nodes=<#nodes>:ppn=<#procs> <this .cmd file> <parameter file>
# 		Expects xyzlist.dat and ts_xyzlist.dat in the working directory
#	2.
#		msub -l nodes=<#nodes>:ppn=<#procs> <this .cmd file> <parameter file> REPO
# 		Expects only xyzlist.dat in the working directory
	  		  	  
#MSUB -V		  
#MSUB -o stdoutmsg	  

PARAMS=$1 # PARAMS="GEN_FF/params.txt.reduced"
EXEC=$2   # EXEC="/p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/SVN_SRC-11-12-18/chimes_md"
BASE=$3   # "run_md.base"
DRIVER=$4

echo "Using parameter file: $PARAMS "
echo "Using executable: $EXEC "
echo "Using base file: $BASE "

# Error checks

if [ ! -e new-get_dumb_ener_subjob.sh ] ; then 
	echo "ERROR: This script requires new-get_dumb_ener_subjob.sh to run."
	exit 0
fi


NODES=$SLURM_JOB_NUM_NODES
PROCS=$SLURM_JOB_CPUS_PER_NODE
PROCS=${PROCS%(*}
NPROC=$[ $NODES * $PROCS ]

echo "Computed nproc, ppn, and total procs: $NODES $PROCS $NPROC"

rm -rf DUMB_ENER_SUBJOB-*	  

awk -v pars="../${PARAMS}" '/PRMFILE/{print;getline;$0=pars}{print}' $BASE > tmp
mv tmp $BASE

TARGETS="ts tight"

if [ "$#" -eq 5 ]; then
	TARGETS="tight"
fi 


for crit in $TARGETS
do

		echo "Working on $crit"
	
		# Read the list of configurations.. format is: <natoms> <n_c> <n_o> </path/to/xyz/file>
		
		TARG=""; if [[ "$crit" == "ts" ]] ; then TARG="${crit}_"; fi
		
		if [ ! -e ${TARG}xyzlist.dat ] ; then

			echo "Warning: Cannot find file ${TARG}xyzlist.dat ... skipping this loop."
			continue
		fi

		readarray -t REPO_FILES < ${TARG}xyzlist.dat
		
		
		# Prepare to break the job into $NPROC separate jobs
		
		NJOBS=${#REPO_FILES[@]} # Number of configurations to compute energies for
		
		if [ $NJOBS -lt $NPROC ] ; then
		
			NPROC=$NJOBS
			
			echo "NJOBS < NPROCS ... will only use $NPROC procs."			
		fi
			
		
		JPERN=`echo " $NJOBS / $NPROC    " | bc`
		NJOB1=`echo " $JPERN * $NPROC + 1" | bc`
		
		echo "	...Processing $NJOBS configurations"
		echo "	...(approx. $JPERN configurations per proc)" 
		
		SUBJOB=0
		LFTOVR=0
		
		rm -rf  DUMB_ENER_SUBJOB-*

		for (( i=0; i<$NJOBS; i++))
		do
			# Which subjob are we on? Set up the subjob .sh driver scripts and the .cmd script
		
			NEWSUBJOB=`echo " ($i ) % $JPERN " | bc` # Returns modulo
			
			CURRJOB=$[ $i + 1 ]
			
			if   [ $CURRJOB -eq $NJOB1 ] ; then
			
				LFTOVR=$SUBJOB
			
				SUBJOB=1
				
			elif [ $CURRJOB -gt $NJOB1 ] ; then
			
				SUBJOB=$[ $SUBJOB + 1 ]
			
			elif [[ $NEWSUBJOB -eq 0 ]] ; then 

				SUBJOB=$[ $SUBJOB + 1 ]
				
				mkdir DUMB_ENER_SUBJOB-${SUBJOB}
				cp ${DRIVER}/utilities/new-get_dumb_ener_subjob.sh DUMB_ENER_SUBJOB-${SUBJOB}
				cp $PARAMS DUMB_ENER_SUBJOB-${SUBJOB}				  
			fi
			
			echo "		...Assigning configuration $i to proc $SUBJOB / $NPROC:"

			echo ${REPO_FILES[$i]} >> DUMB_ENER_SUBJOB-${SUBJOB}/xyzlist.dat
			
		done
		
		if [ $LFTOVR -gt 0 ] ; then
			SUBJOB=$LFTOVR
		fi
		
		echo "	...Counted $SUBJOB subjobs"

		# Run this "TARG's" NPROCS dumb energy calculations 
		
		CURRNODE=-1

		for (( i=1; i<=$SUBJOB; i++))
		do
			ITERATOR=`echo "( $i - 1 ) % $PROCS" | bc` 

			if [ $ITERATOR -eq 0 ] ; then let CURRNODE=CURRNODE+1; fi
			
			cd DUMB_ENER_SUBJOB-${i}
			
			TASK=`pwd`
			TASK="${TASK}/new-get_dumb_ener_subjob.sh 1 ${EXEC} ${BASE} ${i} "

			${TASK} & #srun -N 1 -n 1 ${TASK}

			cd ..
		done

		
		# Process output of TARG's runs
			
		wait	  
		
		rm -f ${TARG}xyzlist.energies ${TARG}xyzlist.energies_normed
		
		for (( i=1; i<=$SUBJOB; i++))	   
		do				   
			cat DUMB_ENER_SUBJOB-${i}/xyzlist.energies		>> ${TARG}xyzlist.energies		  
			cat DUMB_ENER_SUBJOB-${i}/xyzlist.energies_normed 	>> ${TARG}xyzlist.energies_normed	  
			rm -rf DUMB_ENER_SUBJOB-${i}	  
		done				  

done

rm -f all.energies all.energies_normed all.xyzlist.dat

for crit in $TARGETS
do
	TARG=""; if [[ "$crit" == "ts" ]] ; then TARG="${crit}_"; fi

	cat ${TARG}xyzlist.dat			>> all.xyzlist.dat
	cat ${TARG}xyzlist.energies		>> all.energies
	cat ${TARG}xyzlist.energies_normed	>> all.energies_normed 
done

	

	
	
	
