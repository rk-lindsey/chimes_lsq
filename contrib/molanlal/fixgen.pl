#! /usr/bin/perl

# Fix a number of gen files so they start sequentially.
# A list of all files to be fixed needs to be given in sequential order
# (e.g. dump.001, dump.002, etc.). 
# The first dump file (dump.001) is truncated
# if need be so that it ends before the second (dump.002) begins.
# The original version of fixed dump files are saved with a ".old" prefix.
# It is strongly recommended that you save backup copies of the gen files 
# before running this script to prevent possible data loss.
#
# Usage: fixgen.pl --fast <gen files>
#        fixgen.pl --slow <gen files>
#
use strict ;
use warnings ;
use File::Copy ;
use File::Temp qw / tempfile tempdir / ;

sub check_tail ;


my $arg = $ARGV[0] ;
my $fast = 1 ;

if ( $arg eq "--fast" ) {
    $fast = 1 ;
    shift @ARGV ;
} elsif ( $arg eq "--slow" ) {
    $fast = 0 ;
    shift @ARGV ;
}

my @files = @ARGV ;
my $tfirst = -1 ;
my $tlast ;
my $dt ;
my $dtlast ;
my $natom ;

for ( my $j = 0 ; $j < scalar(@files) - 1 ; $j++ ) {

    if ( ! (-f $files[$j]) ) {
	die "$files[$j] is not a file: can't fix\n" ;
    }
    if ( ! (-f $files[$j+1]) ) {
	die "$files[$j+1] is not a file: can't fix\n" ;
    }
    if ( ! (-w $files[$j] ) ) {
	die "$files[$j] is not writable: can't fix\n" ;
    }
    if ( ! (-r $files[$j+1] ) ) {
	die "$files[$j+1] is not readable\n" ;
    }

    open(IN2, "< $files[$j+1]") || die "Could not open $files[$j]\n" ;

#   Grab the time from the next file.
    $_ = <IN2> ;
#   print ;
    $_ = lc($_) ;
    $_ =~ s/time[\s]*step/time/ ;

    # Grab number of atoms from the regular expression.
    if ( /([0-9]+)/ ) {
	$natom = $1 ;
	print "Number of atoms in $files[$j+1] = $natom\n" ;
    } else {
	die "Could not find the number of atoms in $files[$j+1]\n" ;
    }

    if ( /time[\s]*[=]*[\s]*([0-9e\.\+\-]+)/ ) {
#	print ;
	$tfirst = $1 ;
	print "First time for $files[$j+1] = $tfirst \n" ;
    } else {
	die "Could not find a starting time for $files[$j+1]\n" ;
    }
    close(IN2) ;

    if ( check_tail($tfirst, $natom, $files[$j]) ) {
	print "Last time found in $files[$j] was OK\n" ;
	next ;
    }

    # Scan through the current file, and compare the current time
    # to the time of the next file.  Stop when the current time
    # is greater than the time of the next file.
    fix_file($tfirst, $files[$j], $files[$j+1]) ;

}


sub check_tail
# Check the tail of the dump file to see if it matches with
# the beginning of the next one.
# Arguments : 
#   tfirst : first time in next file.
#   filename ; name of current file that is being checked.
{
    my $tfirst   = $_[0] ;
    my $natom    = $_[1] ;
    my $filename = $_[2] ;
    
    my $tnow = 0.0 ;
    my $tlast = -1.0 ;
    my $natom2 = 0 ;
    my $tailfile ;

    my $ntail = $natom + 6 ;
    
    $tailfile = "$filename.tail" ;
    if ( system("tail -$ntail $filename > $tailfile") ) {
	die "Could not find the tail of $filename\n" ;
    }

    open(IN3, "< $tailfile") || die "Could not open $tailfile\n" ;

    $_ = <IN3> || die "$tailfile is empty\n" ;

    my $time_line = $_ ;

#    print "TIME LINE: $time_line\n" ;

    $time_line = lc($time_line) ;
    $time_line =~ s/time[\s]*step/time/ ;
    if ( /([0-9]+)/ ) {
	$natom2 = $1 ;
	if ( $natom != $natom2 ) {
	    print "Number of atoms in $filename = $natom2\n" ;
	    die "$filename has the wrong number of atoms\n" ;
	}
    } else {
	die "Could not find the number of atoms in $filename\n" ;
    }
    if ( $time_line =~ /time[\s]*[=]*[\s]*([0-9e\.\+\-]+)/ ) {
	$tlast = $tnow ;
	$tnow = $1 ;
	print "Last time in $filename = $tnow\n" ;
    }
    if ( $tlast > $tnow ) {
	die "Error: times were out of order in $filename at $tlast\n" ;
    }
    if ( $tnow > $tfirst ) {
	# Need to fix the file.
	print "Need to fix $filename\n" ;
	print "tnow = $tnow, tfirst = $tfirst\n" ;
	unlink $tailfile ;
	return 0 ;
    }
    unlink $tailfile ;
    return(1) ;
}

sub fix_file
# Read through a gen file, and truncate it at the appropriate time.
{
    my $tfirst = $_[0] ;
    my $file1  = $_[1] ;
    my $file2  = $_[2] ;

    my $tnow = 0.0 ;
    my $tlast = -1.0 ;
    my $dt = 0.0 ;
    my $dtlast = -1.0 ;
    my  $printed = 0 ;
    my $tmpname = "tmp.out" ;
    my $tmpfile ;
    my $template ;

    ($tmpfile, $tmpname) = tempfile($template, DIR => ".") ;

    open($tmpfile, ">", $tmpname) ;

    open(IN1, "< $file1") || die "Could not open $file1\n" ;

    while ( <IN1> ) {
	my $time_line = $_ ;
	$time_line = lc($time_line) ;
	$time_line =~ s/time[\s]*step/time/ ;

#	print "TIME LINE: $time_line\n" ;
	if ( $time_line =~ /time[\s]*[=]*[\s]*([0-9e\.\+\-]+)/ ) {
	    $tlast = $tnow ;
	    $tnow = $1 ;
#	    print "Current time = $tnow\n" ;
	}
	if ( $tlast > $tnow ) {
	    die "Error: times were out of order in $file1 at $tlast\n" ;
	}
	if ( $tlast > 0.0 ) {
	    $dtlast = $dt ;
	    $dt = $tnow - $tlast ;
	    if ( $dtlast > 0.0 ) {
		if ( abs($dt - $dtlast) > 0.1 ) {
		    print "Warning: Time step changed\n" ;
		}
	    }
	}
	if ( $tnow > $tfirst ) {
	    # Need to fix the file.
	    if ( $printed ) {
		print "Truncating $file1\n" ;
#		close(TMP) ;
		close($tmpname) ;
		close(IN1) ;

		print "Saving old file to $file1.old\n" ;
		move($file1, "$file1.old")  || 
		    die "Moving $file1 to $file1.old failed\n" ;
		my $tmp = $tmpname ;
		print("Moving $tmp to $file1\n") ;
		move($tmpname, $file1) || 
		    die "Moving $tmpname to $file1 failed\n" ;
		return ;
	    } else {
		print "Error: File $file2 and $file1 do not overlap in time.\n" ;
		print "Check to see if the files are out of order\n" ;
		exit(1) ;
	    }
	} else {
	    $printed = 1 ;
	    print { $tmpfile } $_ ;
	}
    }
    
    if ( $tfirst - $tnow > $dt ) {
	printf("Warning: time gap of %11.4e found between files %s and %s\n",
	       $tfirst - $tnow, $file1, $file2) ;
    }
    print "Last time for $file1 =  $tnow\n" ;
    print "$file1 did not need to be fixed.\n" ;
    close($tmpfile) ;
    close(IN1) ;

}
