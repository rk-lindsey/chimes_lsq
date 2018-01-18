#! /usr/bin/perl
$cfiles = `find . -name '*.C'` ;
$hfiles = `find . -name '*.h'` ;
$cfiles2 = $cfiles ;
$hfiles2 = $hfiles ;
$cfiles2 =~ s/\n/ /g ;
$hfiles2 =~ s/\n/ /g ;
$command = "etags $cfiles2 $hfiles2" ;
print $command, "\n" ;
$status =  system($command) ;
if ( $status > 0 ) {
    exit $status ;
}
open(FILES,">files.txt") || die "Could not open files.txt\n" ;
print FILES $cfiles, $hfiles ;
close(FILES) ;
$command = "gtags -f files.txt" ;
$status =  system($command) ;
if ( $status > 0 ) {
    exit $status ;
}

