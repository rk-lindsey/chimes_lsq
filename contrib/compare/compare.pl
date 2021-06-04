#! /usr/bin/perl
##
## This provides a finite tolerance comparison for Cheetah outputs.
## The Cheetah output file must end in cho or out.
## String case, whitespace, and number format are ignored.
## Certain words or phrases can be ignored.  In that case,
## the line is skipped over in the file containing the word/phrase
## to be ignored, and matching is tried again.
use strict ;
use warnings ;
use File::Basename ;

# Tolerance for number comparison.
my $eps = 1.0e-04 ;
my $verbose = 0 ;

# Skip over uninteresting lines 
my @ignore_words = qw(time date revision cheetah copyright parallel
                     chunk process thread openmp memory mem_track) ;
#my @skip_phrase_lines = ("Size of data/stack/text segments");
my @ignore_phrases = ("allocation limits", 
                      "clock cycles", 
                      "global chebyshev", 
                      "low gas volume", 
                      "reactant library title", 
                      "replicate suffix",
                      "startup input", 
                      "static allocations", 
                      "temperature is over the limit",
                      "the cache file", 
                      "too many iterations",
                      "version tag",
		      "Date",
		      "Warning: r < rmin",
		      "Penalty potential",
		      "..options read.",
		      "..xmat read.",
		      "..yvec read.",
		      "..xmat normalized.",
		      "..yvec normalized.",
		      "..reading single",
		      "..weights read",
		      "cholesky error test",
		      "..reading split",
		      "finished",
		      "using the lasso algorithm",
		      "srun -N"
		      );

my $read1 = 1 ;
my $read2 = 1 ;

if ( scalar(@ARGV) == 0 ) {
  die "Error:  must specify the file to compare\n" ;
} elsif ( scalar(@ARGV) == 1 ) {
  die "Error: must specify the reference output file\n" ;
}

my $file1 = $ARGV[0] ;
my $file2 = $ARGV[1] ;

my $diff_file = $file1 ;
my @suffix_list = qw(.chi .cho .log .out .txt) ;
my ($name, $path, $suffix) = fileparse($diff_file, @suffix_list) ;
$diff_file =~ s/$suffix$/\.diff/ ;

if ( ! open(F1, "<", $file1 ) ) {
  die "Could not open $file1\n" ;
}
if ( ! open(F2, "<", $file2 ) ) {
  die "Could not open $file2\n" ;
}

my $more = 1 ;
my $count = 0 ;
my $different = 0 ;
my ($l1, $l2, @t1, @t2, $tok1, $tok2) ;

LINE: while ( $more ) {
  # Get 1 line from each file.
  $count++ ;
  if ( $read1 > 0 ) {
    if ( $_ = <F1> ) {
      $l1 = $_ ;
      chop $l1 ;
    }
    else {
      last LINE ;
    }
  }
  if ( $read2 >  0 ) {
    if ( $_ = <F2> ) {
      $l2 = $_ ;
      chop $l2 ;
    }
    else {
      last LINE ;
    }
  }
  $read1 = 1 ;
  $read2 = 1 ;
  if ( $verbose ) {
    printf("Comparing lines\n") ;
    printf(":  %s\n", $l1) ;
    printf(":  %s\n", $l2) ;
  }
  # Translate the lines to lower case.
  $l1 =~ tr/A-Z/a-z/ ;
  $l2 =~ tr/A-Z/a-z/ ;

  foreach my $word ( @ignore_words ) {
    if ( ($l1 =~ /$word/) || ($l2 =~ /$word/) ) {
      # Read only from the file that matches the word that
      # should be ignored.
      if ( $l1 =~ /$word/ ) {
        if ( $verbose ) {
          print "Ignoring line $l1\n" ;
        }
        $read1 = 1 ;
      } else {
        $read1 = 0 ;
      }
      if ( $l2 =~ /$word/ ) {
        $read2 = 1 ;
        if ( $verbose ) {
          print "Ignoring line $l2\n" ;
        }
      } else {
        $read2 = 0 ;
      }
      next LINE ;
    }
  }
  
  foreach my $phrase ( @ignore_phrases ) {
    if ( ($l1 =~ /$phrase/) || ($l2 =~ /$phrase/) ) {
      # Read only from the file that matches the phrase that
      # should be ignored.
      if ( $l1 =~ /$phrase/ ) {
        if ( $verbose ) {
          print "Ignoring line $l1\n" ;
        }
        $read1 = 1 ;
      } else {
        $read1 = 0 ;
      }
      if ( $l2 =~ /$phrase/ ) {
        if ( $verbose ) {
          print "Ignoring line $l2\n" ;
        }
        $read2 = 1 ;
      } else {
        $read2 = 0 ;
      }
      next LINE ;
    }
  }

  @t1 = split(" ", $l1 ) ;
  @t2 = split(" ", $l2 ) ;
  
  if ( $#t1 != $#t2 ) {
    $different = 1 ;
    print "Difference found on line $count\n" ;
    print "   :$l1\n" ;
    print "   :$l2\n" ;
    last ;
  }
  
  for ( my $i = 0 ; $i < scalar(@t1) ; $i++ ) {
    $tok1 = $t1[$i] ;
    $tok2 = $t2[$i] ;
  
    print "TOKEN $tok1 $tok2\n" if $verbose ;
  
    $tok1 =~ s/,//g ;
    $tok2 =~ s/,//g ;
  
    if ( $tok1 eq $tok2 ) {
      next ;
    }
    elsif ( $tok1 =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && 
            $tok2 =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ ) {
      # Do a finite tolerance comparison on numbers 
      printf("Allowed number comparison: %s %s\n", $tok1, $tok2) 
      if $verbose ;
      if ( abs($tok1 - $tok2)/(abs($tok1)+1.0e-08) < $eps ) {
        next ;
      } else {
        $different = 1 ;
        printf("Different number comparison: %s %s\n", $tok1, $tok2) ;
        last ;
      }
    }
    else {
      $different = 1 ;
      last ;
    }
  }
  if ( $different != 0 ) {
    print "Difference found on line $count\n" ;
    print "   :$l1\n" ;
    print "   :$l2\n" ;
    last ;
  }
}
if ( $different ) {
  print("Differences were found\n") ;
  system("diff -c $file1 $file2 > $diff_file") ;
}
elsif ( $verbose ) {
  print("No differences were found\n") ;
}
exit($different) ;


