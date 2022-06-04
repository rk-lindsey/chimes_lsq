#! /usr/bin/perl
#
# A program to join several molecular analyzer outputs together after
# parallel evaluation.  Files must be given in sequential time order !
#
#
use strict ;
use warnings ;

sub print_props ;
sub read_reactions ;
sub print_rxns ;
sub read_conc_history ;

{
    my $name ;
    my @rxn_header ;
    my $molname ;
    my @mol_prop_list ;
    my @rxn_prop_list ;
    my @tst_prop_list ;
    my @conc_prop_list ;
    my @poly_prop_list ;
    my @molwt_hist_list ;

    my @frame_count_list ;

    my %conc_props ;
    my $conc_header ;
    my $conc_start  = "Concentration history for" ;
    my $conc_return = "Number of Polymer Molecules" ;
    my $conc_stop = "^\$" ;

    my $filecount ;

    my $molecule_start = "Averages for molecules" ;
    my $molecule_stop = "^\$" ;

    my $tst_start = "Averages for transition" ;
    my $tst_stop = "^\$" ;

    my $rxn_start = "Rate constants in units of" ;
    my $rxn_stop = "^\$" ;

    my @mol_header ;
    my @tst_header ;

    my @molwt_header ;

    my $poly_start = 'Polymer ([\w\s]+) History' ;
    my $poly_stop  = "^\$" ;

    my $molwt_start = "Molecular weight histogram" ;
    my $molwt_stop = "^\$" ;

    $filecount = 0 ;

    foreach my $file ( @ARGV ) {
	@mol_header = () ;
	@tst_header = () ;
	@rxn_header = () ;
	@molwt_header = () ;

	# Set references to empty hashes for each file.
	my $conc_props = {} ;
	my $mol_props  = {} ;
	my $tst_props  = {} ;
	my $rxn_props  = {} ;
	my $poly_props = {} ;

        # molwt_props is a reference to an empty array.
	my $molwt_props = [] ;
	my $frame_count ;

	open(IN, "< $file" ) || die "Could not open $file\n" ;

	read_properties($mol_props, $molecule_start, $molecule_stop,
	               \@mol_header,\$frame_count) ;
	$frame_count_list[$filecount] = $frame_count ;

#	print "MOL_HEADER = $mol_header\n" ;

	# Store an array of references to the molecular properties.

	$mol_prop_list[$filecount] = $mol_props ;

	read_properties($tst_props, $tst_start, $tst_stop,
	                \@tst_header,\$frame_count) ;

	# Store an array of references to the transition state properties.
	$tst_prop_list[$filecount] = $tst_props ;

	read_molwt_histo($molwt_props, $molwt_start, $molwt_stop, \@molwt_header) ;
	$molwt_hist_list[$filecount] = $molwt_props ;

	read_reactions($rxn_props, $rxn_start, $rxn_stop, \@rxn_header) ;

	$rxn_prop_list[$filecount] = $rxn_props ;

	read_conc_history($conc_props, $conc_start, $conc_stop,
			  $conc_return, \$conc_header) ;
	$conc_prop_list[$filecount] = $conc_props ;

	
	read_poly_history($poly_props, $poly_start, $poly_stop) ;
	$poly_prop_list[$filecount] = $poly_props ;

	$filecount++ ;
    }
    print_props(\@mol_prop_list, $filecount,\@mol_header,\@frame_count_list) ;
    print_props(\@tst_prop_list, $filecount,\@tst_header,\@frame_count_list) ;
    print_molwt_histo(\@molwt_hist_list, $filecount, \@molwt_header,\@frame_count_list) ;
    print_rxns(\@rxn_prop_list, $filecount,\@rxn_header,\@frame_count_list) ;
    print_conc_history(\@conc_prop_list, $filecount, $conc_header) ;
    print_poly_history(\@poly_prop_list, $filecount) ;
}

sub print_props
## Print the properties in the given property hash.
#
# Arguments:  mol_prop_list
#             filecount
#             header
#             frame_count_list
{
    my @avg_prop ;
    my $molecule ;
    my $fields ;
    my %avg_count ;

    my $mol_prop_list  = $_[0] ;
    my $filecount      = $_[1] ;
    my $header         = $_[2] ;
    my $frame_ref      = $_[3] ;

#    print "\n" ;
    print @$header ;

    my %local_props ;
    my $field ;
    die "Filecount needs to be positive" if $filecount <= 0 ;

#    print "PRINT_PROPS: filecount = $filecount\n" ;

    my $local_props ;
    my @molecules ;

    # Loop over files.  Average input data.
    my $total_frames = 0 ;

    for ( my $k = 0 ; $k < $filecount ; $k++ ) {
	my $prop2 = $mol_prop_list->[$k] ;
	@molecules = keys ( %{$prop2} ) ;

#	print "PROP_LIST: $mol_prop_list->[$k]\n" ;

	foreach $molecule ( @molecules ) {
#	    printf ("PRINT_PROP: %s", $molecule) ;
	    my $field2 = $prop2->{$molecule} ;

	    if ( defined ( $field2 ) ) {
		# Loop over each data entry.
		$avg_count{$molecule} += $frame_ref->[$k] ;
		for ( my $l = 0 ; $l < scalar(@{$field2}) ; $l++ ) {
		    $local_props{$molecule}[$l] += $field2->[$l] * $frame_ref->[$k];
		}
	    }
	}
	$total_frames += $frame_ref->[$k] ;
    }

    # Print out averaged data.
    @molecules = keys ( %local_props ) ;
    @molecules =  sort { prop_byconc($a, $b, \%local_props) } @molecules ;

    # Average concentrations across all files, whether the species is present or not.
    # Concentrations are zero if the species is not present in a file.
    # Average intensive properties (size, mol wt., etc.) by only those files where the
    # species is present.

    foreach $molecule ( @molecules ) {
	$field = $local_props{$molecule} ;
	for ( my $l = 0 ; $l < 2 ; $l++ ) {
	    $field->[$l] /= $total_frames
	}
	for ( my $l = 2 ; $l < scalar(@{$field}) ; $l++ ) {
	    $field->[$l] /= $avg_count{$molecule} ;
	}
    }

    # Print out molecular properties.
    foreach $molecule ( @molecules ) {
	printf ("%70s", $molecule) ;
	$field = $local_props{$molecule} ;
	next if ! defined($field) ;
	if ( $avg_count{$molecule} > 0 ) {
	    for ( my $l = 0 ; $l < scalar(@{$field}) ; $l++ ) {
		if ( ! defined($field->[$l]) ) {
		    # Sanity check.
		    die "Field not defined for $molecule\n"; 
		} 
		# Report average.
		printf("%9.4f ", $field->[$l]) ;
	    }
	}
	print "\n" ;
    }
}


sub print_rxns
## Print the reactions
#
# Arguments:  mol_prop_list
#             filecount
#             rxn_header
{
    my $reaction ;
    my $fields ;

    my $rxn_prop_list  = $_[0] ;
    my $filecount      = $_[1] ;
    my $rxn_header     = $_[2] ;
    my $frame_ref      = $_[3] ;

## TO DO: Average over quantities in reaction header.
###    print @{$rxn_header} ;
    my $line = scalar(@{$rxn_header}-3) ;
    for ( my $j = $line ; $j < scalar(@$rxn_header) ; $j++ ) {
	print $rxn_header->[$j] ;
    }

    die "Filecount needs to be positive" if $filecount <= 0 ;

#    print "PRINT_PROPS: filecount = $filecount\n" ;

    # Loop over reactions.
    my %total_count ;
    my %rate ;
    
    ## Average across all reactions found in files.  Keep track of the number
    ## of times a reaction has appeared.
    my $total_frames = 0.0 ;
    for ( my $k = 0 ; $k < $filecount ; $k++ ) {
	my $prop2 = $rxn_prop_list->[$k] ;

	foreach $reaction ( keys ( %{$prop2} ) ) {	
	    # Average across each file.
	    my $field2 = $prop2->{$reaction} ;
	    
	    if ( defined ( $field2 ) ) {
		$total_count{$reaction} += $field2->[0] ;
		$rate{$reaction} += $field2->[1] * $frame_ref->[$k] ;
	    }
	}
	$total_frames += $frame_ref->[$k] ;
    }

    # Use customized sort routine syntax.
    my @sorted_reactions = sort { rxn_bycount($a, $b, \%total_count) } keys ( %total_count ) ;
    foreach $reaction ( @sorted_reactions ) {
	if ( $total_frames > 0 ) {
	    $rate{$reaction} /= $total_frames ;
	} else {
	    die "Error in print_rxn: Number of frames = 0\n" ;
	}
	# Report average.
	printf("%9d %11.4e %s\n", $total_count{$reaction}, $rate{$reaction}, $reaction) ;
    }

}


sub print_conc_history
## Print the concentration history
#
# Arguments:  conc_prop_list
#             filecount
#             conc_header
{
    my $molecule ;
    my $fields ;

    my $conc_prop_list  = $_[0] ;
    my $filecount      = $_[1] ;
    my $conc_header     = $_[2] ;
    my $l ;

    my %maxconc ;       # Maximum concentration across all files.

    die "Filecount needs to be positive" if $filecount <= 0 ;

#    print "PRINT_PROPS: filecount = $filecount\n" ;

    my %local_props ;

    # Loop over molecules.
    # Use customized sort routine syntax.


    my $count = 0 ;
    my $last_time = 0 ;
    my $lmax ;
    my $tmax ;

    for ( my $k = 0 ; $k < $filecount ; $k++ ) {
#	print "FILE = $k\n" ;
	my $prop2 = $conc_prop_list->[$k] ;
	my @molecules = keys ( %{$prop2} ) ;	
	$lmax = 0 ;
	$tmax = 0 ;
	foreach $molecule ( @molecules ) {
	    my $field2 = $prop2->{$molecule} ;

	    # Loop over each data entry.
	    if ( defined ( $field2 ) ) {
#		print "FIELD2 $field2\n" ;
		my $f2sz = scalar(@$field2) ;
#		print "FIELD2 size $f2sz\n" ;
		for ( $l = 0 ; $l < scalar(@$field2) ; $l++ ) {
		    # Increment the time from the last file.
		    $field2->[$l][0] += $last_time ;
		    for ( my $j = 0 ; defined($field2->[$l][$j]) ;
			  $j++ ) {
			$local_props{$molecule}[$count+$l][$j] = 
			    $field2->[$l][$j] ;
#			printf("CONC LOCAL_PROPS: Molecule = %s Count = %d Index1 = %d Index2 = %d Value = %11.4f\n",
#			       $molecule, $count, $count+$l, $j, $field2->[$l][$j]) ;
		    }
		    if ( defined($maxconc{$molecule}) ) {
			if ( $field2->[$l][1] > $maxconc{$molecule} ) {
			    $maxconc{$molecule} = $field2->[$l][1] ;
			}
		    } else {
			$maxconc{$molecule} = $field2->[$l][1] ;
		    }
		}
		if ( $l > $lmax ) {
		    $lmax = $l ;
		    $tmax = $field2->[$l-1][0] - $last_time ;
		}
	    }
#	    printf("LMAX = %d\n", $lmax) ;
	}
	## Calculate the time and array index offset for each file.
	$count += $lmax ;
	$last_time += $tmax ;
#	print "COUNT = $count, LAST_TIME = $last_time\n" ;
    }

    my @sorted_molecules = sort { hist_byconc($a, $b, \%maxconc) } keys ( %local_props ) ;

    foreach $molecule ( @sorted_molecules ) {

	printf ("\nConcentration history for %s\n", $molecule) ;
	printf ("\nMaximum concentration = %12.4e moles/cc\n", $maxconc{$molecule}) ;
	$fields = $local_props{$molecule} ;

	# Loop over each data entry.
	my $last_time = 0 ;
	print $conc_header ;

	# Concatenate results from each file.
	my $field2 = $local_props{$molecule} ;

	if ( defined ( $field2 ) ) {
#		print "FIELD2 $field2\n" ;
	    my $f2sz = scalar(@$field2) ;
#		print "FIELD2 size $f2sz\n" ;
	    for ( $l = 0 ; $l < scalar(@$field2) ; $l++ ) {
		# Increment the time from the last file.
		# Could have undefined values if a molecule is missing from some files.
		if ( defined($field2->[$l]) ) {
		    for ( my $j = 0 ; defined($field2->[$l][$j]) ;
			  $j++ ) {
			if ( $j <= 2 ) {
			    printf("%11.4e ", $field2->[$l][$j]) ;
			} else {
			    printf("%7.4f", $field2->[$l][$j]) ;
			}
		    }
		    printf("\n") ;
		}
	    }
	}
    }

}


sub print_molwt_histo
## Print the molwtentration histo
#
# Arguments:  molwt_prop_list
#             filecount
#             molwt_header
#             frame_count_list
{
    my $molecule ;
    my $fields ;

    my $molwt_prop_list  = $_[0] ;
    my $filecount        = $_[1] ;
    my $molwt_header     = $_[2] ;
    my $frame_ref        = $_[3] ;
    my $total_frames ;

    my $l ;


    die "Filecount needs to be positive" if $filecount <= 0 ;

#    print "PRINT_PROPS: filecount = $filecount\n" ;

    my @local_props ;

    # Loop over molecules.
    # Use customized sort routine syntax.


    $total_frames = 0 ;
    for ( my $k = 0 ; $k < $filecount ; $k++ ) {
#	print "FILE = $k\n" ;
	my $field2 = $molwt_prop_list->[$k] ;
	# Loop over each data entry.
	if ( defined ( $field2 ) ) {
	    my $f2sz = scalar(@$field2) ;
#		print "FIELD2 size $f2sz\n" ;
	    for ( $l = 0 ; $l < scalar(@$field2) ; $l++ ) {
		my $idx = int($field2->[$l][0]) ;
		# Bin by molecular weight.
		for ( my $j = 1 ; defined($field2->[$l][$j]) ;
		      $j++ ) {
		    $local_props[$idx][$j] += $field2->[$l][$j] * $frame_ref->[$k] ;
#			printf("MOLWT LOCAL_PROPS: Count = %d Index1 = %d Index2 = %d Value = %11.4f\n",
#			       $molecule, $count, $count+$l, $j, $field2->[$l][$j]) ;
		}
	    }
	}
	$total_frames += $frame_ref->[$k] ;
    }
    ## Calculate the time and array index offset for each file.

    # Loop over each data entry.

    print "\n" ;
    print @$molwt_header ;

#   Average results from each file.

#		print "LOCAL_PROPS $local_props\n" ;
    my $f2sz = scalar(@local_props) ;
#		print "LOCAL_PROPS size $f2sz\n" ;
    for ( $l = 0 ; $l < scalar(@local_props) ; $l++ ) {
	# Increment the time from the last file.
	# Could have undefined values if a molecule is missing from some files.
	if ( defined($local_props[$l]) ) {
	    for ( my $j = 1 ; defined($local_props[$l][$j]) ;
		  $j++ ) {
		$local_props[$l][$j] /= $total_frames ;
		printf("%11.4e %11.4e ", $l, $local_props[$l][$j]) ;
	    }
	    printf("\n") ;
	}
    }
    print "\n" ;
}



sub print_poly_history
## Print the polymer concentration history
#
# Arguments:  poly_prop_list
#             filecount
{
    my $property ;
    my $fields ;

    my $poly_prop_list  = $_[0] ;
    my $filecount      = $_[1] ;
    my $l ;

    die "Filecount needs to be positive" if $filecount <= 0 ;

#    print "PRINT_PROPS: filecount = $filecount\n" ;

    my $local_props = $poly_prop_list->[0] ;

    # Loop over properties.
    # Use customized sort routine syntax.
    my @sorted_properties = sort keys ( %$local_props ) ;
    my $last_time ;

    foreach $property ( @sorted_properties ) {

	if ( $property eq "Number" ) {
	    print "\nNumber of polymer molecules\n" ;
	} else {
	    printf ("\nPolymer history for %s\n", $property) ;
	}
	$fields = $local_props->{$property} ;

	# Loop over each data entry.
	$last_time = 0 ;
	my $poly_header     = $fields->[0] ;

	print $poly_header ;

	for ( my $k = 0 ; $k < $filecount ; $k++ ) {
	    # Average across each file.
	    my $prop2 = $poly_prop_list->[$k] ;
	    my $field2 = $prop2->{$property} ;

	    if ( defined ( $field2 ) ) {
#		print "FIELD2 $field2\n" ;
		my $f2sz = scalar(@$field2) ;
#		print "FIELD2 size $f2sz\n" ;
		for ( $l = 1 ; $l < scalar(@$field2) ; $l++ ) {
                    # Increment the time from the last file.
		    $field2->[$l][0] += $last_time ;
		    for ( my $j = 0 ; defined($field2->[$l][$j]) ;
			  $j++ ) {
			printf("%12.5e", $field2->[$l][$j]) ;
		    }
		    printf("\n") ;
		}
		$last_time += $field2->[$l-1][0] - $last_time ;
#		print "LAST_TIME = $last_time\n" ;
	    }
		    
	}
    }

}


sub read_properties
#
{
    my $molname ;
    my $mol_props  =  $_[0] ;
    my $start_line =  $_[1] ;
    my $stop_line  =  $_[2] ;
    my $header     =  $_[3] ;
    my $frame_count = $_[4] ;

    while ( <IN> ) {
	push @$header, $_ ;
	if ( /frames were read in/ ) {
	    my @f = split(" ") ;
	    ${$frame_count} = $f[0] ;
	}
	last if ( /$start_line/ ) ;

    }
    $_ = <IN> ;
    push @$header, $_ ;

#    print "HEADER = @$header\n" ;


    while ( <IN> ) {
	last if /$stop_line/ ;
#	print "READ_PROPERTIES: $_" ;

#       This regular expression finds the molecule name and
#       puts it into $1.
	if ( /(^[a-zA-Z0-9\s\(\-\)]+)[0,1]\./) {
	    $molname = $1 ;
#	    print "Molecule name = |$molname|\n" ;
	    my $field = $_ ;

	    my $len = length($molname) ;
#	    print "LEN = $len\n" ;

	    $field = substr($_,  $len) ;
	    
#	    print "Data file = $field\n" ;
	    my @data = split(" ", $field) ;
	    $molname =~ s/  //g ;
	    $molname =~ s/^ // ;
	    $molname =~ s/ $// ;
	    $mol_props->{$molname} = [ @data ] ;
	} else {
#	    print "NO MOLECULE\n" ;
#	    print ;
	    ;
        }
    }
}



sub read_reactions
#
{
    my $rxnname ;
    my $rxn_props  = $_[0] ;
    my $start_line = $_[1] ;
    my $stop_line  = $_[2] ;
    my $rxn_head   = $_[3] ;
    my @data ;
	    

    while ( <IN> ) {
	push @{$rxn_head}, $_ ;
	last if ( /$start_line/ ) ;
    }
    $_ = <IN> ;
    push @{$rxn_head}, $_ ;

#    print "RXN HEADER = |@{$rxn_head}|\n" ;


    while ( <IN> ) {
	last if /$stop_line/ ;
#	print "READ_REACTIONS: $_" ;

#	if ( /(^[0-9\s+-e\.]+)([A-Z0-9a-z]+$)/) {
#	if ( /(^[0-9\s\.e+q]+)/) {
	my @f = split(" ",$_) ;

	$data[0] = shift(@f) ;
	$data[1] = shift(@f) ;

	$rxnname = join " ", @f ;
#	print "Reaction name = |$rxnname|\n" ;
#	print "DATA FIELD    = |@data|\n" ;
	$rxnname =~ s/  //g ;
	$rxnname =~ s/^ // ;
	$rxnname =~ s/ $// ;
	$rxn_props->{$rxnname} = [ @data ] ;
    }
}

sub read_conc_history
# Read the concentration histories for each molecule, putting
# each one in a separate hash.
{
    my $concname ;
    my $conc_props  = $_[0] ;
    my $start_line  = $_[1] ;
    my $stop_line   = $_[2] ;
    my $return_line = $_[3] ;
    my $conc_head   = $_[4] ;
    my @data ;
    my $molecule ;
    my @props ;

#    print "RETURN LINE = |$return_line|\n" ;
    while ( ! eof IN ) {
	if ( /$return_line/ ) {
#	    print "FOUND RETURN LINE\n" ;
	    return ;
	}
	while ( <IN> ) {
	    if ( /$return_line/ ) {
#		print "FOUND CONC RETURN LINE\n" ;
		return ;
	    } elsif ( /$start_line/ ) {
		$molecule = $_ ;
		chomp $molecule ;
		$molecule =~ s/$start_line// ;
		$molecule =~ s/^ // ;
		$molecule =~ s/ $// ;
		for ( my $j = 0 ; $j < 4 ; $j++ ) {
		    $_ = <IN> ;
		}
		last ;
	    } else {
#		print ;
		;
	    }
	}
	$$conc_head = <IN> ;

#	print "MOLECULE = $molecule\n" ;
#	print "CONC HEADER = $$conc_head\n" ;

	my $count = 0 ;
	while ( <IN> ) {
	    last if /$stop_line/ ;
	    if ( /$return_line/ ) {
		print "FOUND CONC RETURN LINE\n" ;
		return ;
	    }

#	    print "READ_CONC_HISTORY: $_" ;

	    my @f = split(" ",$_) ;

#	    print "DATA FIELD    = |@f|\n" ;
	    for ( my $k = 0 ; $k < scalar(@f) ; $k++ ) {
		$props[$count][$k] = $f[$k] ;
	    }
	    $count++ ;
	}
	for ( my $j = 0 ; $j < $count ; $j++ ) {
	    for ( my $k = 0 ; defined($props[$j][$k]) ; $k++ ) {
		$conc_props->{$molecule}[$j][$k] = $props[$j][$k] ;
	    }
	}
#
#	my $prop_ref = $conc_props->{$molecule} ;
#	for ( my $j = 0 ; $j < $count ; $j++ ) {
#	    print "PROPS: $prop_ref->[$j]\n" ;
#	    for ( my $k = 0 ; defined($prop_ref->[$j][$k]) ; $k++ ) {
#		print "PROPS VAL: $prop_ref->[$j][$k] " ;
#	    }
#	    print ("\n") ;
#	}
    }
}


sub read_molwt_histo
# Read the molecular weight histogram.
{
    my $molwtname ;
    my $molwt_props  = $_[0] ;
    my $start_line   = $_[1] ;
    my $return_line  = $_[2] ;
    my $molwt_head   = $_[3] ;
    my @data ;
    my @props ;

#    print "RETURN LINE = |$return_line|\n" ;
    while ( <IN> ) {
	if ( /$start_line/ ) {
	    push @$molwt_head, $_ ;
	    last ;
	} else {
#		print ;
	    ;
	}
    }
    $_ = <IN> ;
    push @$molwt_head, $_ ;

#    print "MOLWT HEADER = @$molwt_head\n" ;

    my $count = 0 ;
    while ( <IN> ) {
	if ( /$return_line/ ) {
#	    print "FOUND MOLWT RETURN LINE\n" ;
	    last ;
	}

#	print "READ_MOLWT_HISTORY: $_" ;

	my @f = split(" ",$_) ;

#	print "DATA FIELD    = |@f|\n" ;
	for ( my $k = 0 ; $k < scalar(@f) ; $k++ ) {
	    $molwt_props->[$count][$k] = $f[$k] ;
	}
	$count++ ;
    }
}



sub read_poly_history
# Read the polyentration histories for each molecule, putting
# each one in a separate hash.
{
    my $polyname ;
    my $poly_props  = $_[0] ;
    my $start_line  = $_[1] ;
    my $stop_line   = $_[2] ;
    my $last_line   = $_[3] ;
    my $poly_head   ;
    my @data ;
    my $property ;
    my @props ;

    my $first_time = 0 ;

    $property = "Number" ;
    $first_time =  1 ;
    while ( ! eof IN ) {
	while ( <IN> ) {
	    if ( /$start_line/ || $first_time ) {
#               Grab the property from the regular expression.
		if ( ! $first_time ) {
		    $property = $1 ;
		    if ( ($property =~ /Molecular Weight/) ) {
			$_ = <IN> ;
		    }
		} else {
		    $first_time = 0 ;
		}
		$_ = <IN> ;
		chomp $property ;
		$property =~ s/^ // ;
		$property =~ s/ $// ;
		last ;
	    } else {
		print ;
	    }
	}
	$poly_head = $_ ;

#	print "PROPERTY = $property\n" ;
#	print "POLY HEADER = $poly_head\n" ;

	my $count = 1 ;
	while ( <IN> ) {
	    last if /$stop_line/ ;

#	    print "READ_POLY_HISTORY: $_" ;

	    my @f = split(" ",$_) ;

#	    print "DATA FIELD    = |@f|\n" ;
	    for ( my $k = 0 ; $k < scalar(@f) ; $k++ ) {
		$props[$count][$k] = $f[$k] ;
	    }
	    $count++ ;
	}
	$poly_props->{$property}[0] = $poly_head ;
	for ( my $j = 1 ; $j < $count ; $j++ ) {
	    for ( my $k = 0 ; defined($props[$j][$k]) ; $k++ ) {
		$poly_props->{$property}[$j][$k] = $props[$j][$k] ;
	    }
	}
    }
}

sub prop_byconc 
# Used for sorting property hashes by concentration.
{
    $a = $_[0] ;
    $b = $_[1] ;
    my $local_props = $_[2] ;

    my $fielda = $local_props->{$a} ;
    my $fieldb = $local_props->{$b} ;

    $$fieldb[0] <=> $$fielda[0] ;
}


sub hist_byconc 
# Used for sorting concentration histories.
{
    $a = $_[0] ;
    $b = $_[1] ;
    my $maxconc = $_[2] ;

    my $fielda = $maxconc->{$a} ;
    my $fieldb = $maxconc->{$b} ;

    $fieldb <=> $fielda ;
}


sub rxn_bycount 
# Used for sorting property hashes by concentration.
{
    $a = $_[0] ;
    $b = $_[1] ;
    my $bycount = $_[2] ;

    my $fielda = $bycount->{$a} ;
    my $fieldb = $bycount->{$b} ;

    $fieldb <=> $fielda ;
}

