#! /usr/bin/perl
# Program to skip the output header.
# nskip is the number of lines to skip.

$nskip = 3 ;
for ( $j = 1 ; $j <= $nskip ; $j++ ) {
    $_ = <> ;
}
while ( <> ) {
    print ;
}
