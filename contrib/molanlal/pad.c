/* Parallel Analysis Driver (pad) 
   
   This program reads a list of commands, and executes the commands
   across a number of MPI processes.  

   Usage: srun -n <nprocs> pad <infile>

   infile is an input file with a list of shell commands to execute.
   The first line of inline gives the number of commands per processor,
   which is fixed.  Usually this is the number of sequential steps in
   the molecular analysis.
   
   The next lines give commands to be executed in parallel.  It is
   assumed that each analysis program is MPI serial, but it may be
   thread parallel.  If the analysis program is thread parallel, an
   appropriate number of CPU's per mpi process must be allocated
   in the msub script.

   Here is an example input file:
---------------------------------
2
molanal.new Benzene.md.gen > molanal.1.out
findmolecules.pl molanal.1.out > findmolecules.1.out
molanal.new Benzene.md.gen > molanal.2.out
findmolecules.pl molanal.2.out > findmolecules.2.out
molanal.new Benzene.md.gen > molanal.3.out
findmolecules.pl molanal.3.out > findmolecules.3.out
---------------------------------
This could be executed in parallel by up to 3 mpi processes.  Each
process will execute 2 commands in sequence.

*/

#include <stdio.h>
#include <mpi.h>
#include <ctype.h>

#define BUFSIZE 512

int main(int argc, char *argv[])
{
  int numprocs, rank ;
  int nlines ;
  FILE *fcom ;
  int j, k ;
  char command[BUFSIZE] ;
  int ncommand ;
  int execute = 0 ;
  int blank ;

  if ( argc == 2 ) {
    fcom = fopen(argv[1],"r") ;
    if ( fcom == NULL || ferror(fcom) ) {
      printf("Could not open %s\n", argv[0]) ;
      exit(1) ;
    }
  } else {
    printf("argc = %d\n", argc) ;
    printf("Usage:  pad <driver file>\n") ;
    exit(1) ;
  }

  MPI_Init(&argc, &argv) ;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs) ;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank) ;

  /* printf("Process %d out of %d\n", rank, numprocs) ; */

  fgets(command, BUFSIZE, fcom) ;

  if ( sscanf(command,"%d", &ncommand) <= 0 ) {
    printf("Could not read the number of commands per process\n") ;
    exit(1) ;
  }
  if ( ncommand <= 0 ) {
    printf("The number of commands per process must be positive\n" ) ;
    exit(1) ;
  }

  /** Read in lines.  **/
  for ( j = 0 ; fgets(command, BUFSIZE,fcom) != NULL ; j++ ) {

    /* Decide whether to execute the line */
    if ( (j / ncommand) % numprocs == rank ) {
      execute = 1 ;
    } else {
      execute = 0 ;
    }
    blank = 1 ;
    for ( k = 0 ; k < BUFSIZE && command[k] != '\0' ; k++ ) {
      if ( ! isspace(command[k]) ) {
	blank = 0 ;
	break ;
      }
    }
    if ( blank ) {
      printf("Warning: skipping a blank input line\n") ;
    }
    if ( execute ) {
      /* Execute the line read in */
      printf("Executing command: %s", command) ;
      if ( system(command) != 0 ) {
	printf("Error in executing: %s", command) ;
	break ;
      }
    }
  }
  
  MPI_Finalize() ;
}

