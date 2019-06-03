#!/bin/bash							
#MSUB -N ALC-X-MD						
#MSUB -l nodes=8:ppn=36						
#MSUB -l walltime=2:00:00					
#MSUB -q pbatch							
#MSUB -m abe							
#MSUB -A pbronze							
#MSUB -V								
#MSUB -o stdoutmsg						
srun -n 288 /p/lscratchrza/rlindsey/RC4B_RAG/11-12-18/SVN_SRC-11-12-18//chimes_md run_md.in > run_md.out
