#! /usr/bin/perl

# A program to run the molecular analyzer on all dftb output files
# and collect the results.

# Usage: molanal.pl <gen name> <chr name> <maximum file>
# 
# The molanal.new program will be run on each gen file of the
# form <gen name>.000 ... <gen name>.<maximum file>
#
# xyz files are created with names molanal.001.out, etc.
#
# Larry Fried, 5/30/2005
#
$gen_base = $ARGV[0] ;
$chr_base = $ARGV[1] ;
$max = $ARGV[2] ;

for ( $i = 1 ; $i <= $max ; $i++ ) {
    $ext = sprintf("%05d",$i) ;
    $gen_name = $gen_base . $ext ;
    $chr_name = $chr_base . $ext ;
    system("molanal.new $gen_name $chr_name > molanal.$ext.out") ;
    rename ("molanal.xyz", "molanal.${ext}.xyz") ;
}

