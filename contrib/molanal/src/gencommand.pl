#!/usr/bin/perl
# gencommand.pl generates a list of molecular analysis commands
# for use by pad.c in parallel analysis.
#
# This is meant to be a simple program that is modified
# on a case by case basis for analysis.

use strict ;
use warnings ;

###############################################################
#
# USER-DEFINED DATA
#
#
my $start=1;                     # Starting file number.
my $end=21 ;                     # Ending file number.
my $step=2 ;                      # Step between file number to process
my $bindir="/g/g17/fried/bin" ;  # Name of executable directory.
my $jobname="tatb.mss2t.9kpms" ;  # Prefix for dump and output files.
my $num_threads=4 ;              # Number of OpenMP threads.

# Commands to execute.  Note the printf-style format, and the end of line characters "\n".
#my @commands = 
#(
#"$bindir/lammpsxyztogen.pl $jobname.%03d.dump.x\n",
#"export OMP_NUM_THREADS=$num_threads ; $bindir/molanal.new $jobname.%03d.x.gen >$jobname.%03d.mol\n", 
#"$bindir/findmolecules.pl $jobname.%03d.mol > $jobname.%03d.find\n"
#) ;

my @commands = 
    (
#"$bindir/lammpsxyztogen.pl $jobname.%03d.dump.x\n",
"$bindir/fixgen.pl $jobname.%03d.x.gen $jobname.%03d.x.gen >& fixgen.%03d.log\n" 
     ) ;
# END USER-DEFINED DATA
################################################################
my $commands_per_cpu = scalar(@commands) ;       # Number of commands per cpu.

print "$commands_per_cpu\n" ;

for ( my $j = $start ; $j <= $end ; $j += $step ) { 
    foreach my $com ( @commands ) {
#       Use this for other commands.
#	printf($com, $j, $j) ;

#       Use this for fixgen.pl
	printf($com, $j, $j+1, $j) ;
    }
}
