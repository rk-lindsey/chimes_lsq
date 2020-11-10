#! /usr/bin/perl
# This is a program to compare the A matrix of the current code to the older Koziol code, for H2O
# test problems.
# The codes use different orderings of columns of the A matrix.  This comparison program
# compares reordered columns.
#
# Usage:  ordered comparison: compare_A_mat.pl --coulomb --fit_pover <Lucas A Mat> <Generalized A mat>
#         unordered comparison: compare_A_mat.pl --unordered <Lucas A Mat> <Generalized A mat>
#
# Example 1:  Compare 2-body Cheby A matrix between codes fitting to Coulomb and overcoordination:
#  run test_suite-lsq/h2o-2bcheby in new code.
#  run "make cheby3" in example/h2o-lsq directory of old code.
#  fmatch.layer/contrib/compare/compare_A_mat.pl --coulomb --fit_pover fmatch/example/h2o-lsq/A.txt fmatch.layer/test_suite-lsq/h2o-2bcheby/A.txt.
#
# Example 2:  Compare 3-Body Cheby A matrix between codes fitting to Coulomb but not overcoordination:
#  run test_suite-lsq/h2o-3bcheby2 in new code directory.
#  run "make cheby_fit3" in old code directory.
#  fmatch.layer/contrib/compare/compare_A_mat.pl --unordered fmatch/example/h2o-3b/A.txt fmatch.layer/test_suite-lsq/h2o-3b/A.txt
#
use strict ;
use warnings ;
use Getopt::Long qw(GetOptions) ;

# Set coulomb to 1 if charges are fit.
my $coulomb = 0 ;

# Set fit_pover to 1 if overcoordination parameter is fit.
my $fit_pover = 0 ;

# Set unordered to 1 if an unordered comparison is desired (example: 3-body chebyshev).
my $unordered = 0 ;

GetOptions(  'coulomb' => \$coulomb,
			  'fit_pover' => \$fit_pover,
				 'unordered' => \$unordered) ;

open(IN1, "<", $ARGV[0]) or 
  die "Could not open $ARGV[0]\n" ;

open(IN2, "<", $ARGV[1]) or
  die "Could not open $ARGV[1]\n" ;

my $eps = 1.0e-05 ;

my $count = 1 ;
while ( <IN1> ) {
  
  my (@f1, @f2) ;

  @f1 = split(" ") ;

  $_ = <IN2> ;

  @f2 = split(" ") ;

  if ( scalar(@f1) != scalar(@f2) ) {
	 my $sf1 = scalar(@f1) ;
	 my $sf2 = scalar(@f2) ;
    die "Different number of columns found: $sf1 vs. $sf2\n" ;
  }
  my $cols = scalar(@f1) - 3 * $coulomb - $fit_pover ;
  print "Number of short-range columns found = $cols\n" ;


  if ( $unordered ) {
	 # Give up on reordering the numbers.  Simply sort and compare in order.
	 @f1 = sort { $a <=> $b } @f1 ;
	 @f2 = sort { $a <=> $b } @f2 ;
	 for ( my $i = 0 ; $i < scalar(@f1) ; $i++ ) {
		if ( abs($f1[$i] - $f2[$i]) > $eps ) {
		  die "Difference found in number (|$f1[$i]| vs. |$f2[$i]\) of row $count\n" ;
		}
    }	
	 print "Passed unordered check\n" ;
  } else {

	 # Reorder the columns so that we can compare between codes.
	 my $s1 = $cols / 3 ;
	 print "Cheby order = $s1\n" if ( $count == 1 ) ;

	 if ( $s1 <= 1 ) {
		die "Bad Cheby order found: $s1\n" ;
	 }

	 for ( my $i = 0 ; $i < $s1 ; $i++ ) {
		if ( abs($f1[$i] - $f2[$i]) > $eps ) {
		  die "Difference found in column $i\n" ;
		}
	 }

	 for ( my $i = 0 ; $i < $s1 ; $i++ ) {
		if ( abs($f1[$i+$s1] - $f2[$i+2*$s1]) > $eps ) {
		  my $col = $i + $s1 ;
		  my $col2 = $i + 2*$s1 ;
		  die "Difference found in file 1 column $col file 2 $col2\n" ;
		}
	 }

	 for ( my $i = 0 ; $i < $s1 ; $i++ ) {
		if ( abs($f1[$i+2 * $s1] - $f2[$i+$s1]) > $eps ) {
		  my $col = $i + 2 * $s1 ;
		  my $col2 = $i + $s1 ;
		  die "Difference found in file 1 column $col file 2 $col2\n" ;
		}
	 }

	 if ( $coulomb ) {
		my $coul_start = 3 * $s1 ;
		if ( abs($f1[$coul_start] - $f2[$coul_start]) > $eps ) {
		  die "Difference found in coulomb column $coul_start \n" ;
		}

		if ( abs($f1[$coul_start+1] - $f2[$coul_start+2]) > $eps ) {
		  die "Difference found in coulomb column $coul_start \n" ;
		}

		if ( abs($f1[$coul_start+2] - $f2[$coul_start+1]) > $eps ) {
		  die "Difference found in coulomb column $coul_start \n" ;
		}
		print "Coulomb parameters for line $count passed\n" ;
	 }

	 if ( $fit_pover ) {
		my $pover_col = 3 * $s1 + 3 * $coulomb ;
		if ( abs($f1[$pover_col] - $f2[$pover_col]) > $eps ) {
		  die "Difference found in overcoord column $pover_col \n" ;
		}
		print "Over-coord parameters for line $count passed\n" ;
	 }
  }

  print "Line $count passed\n" ;
  $count++ ;
}
