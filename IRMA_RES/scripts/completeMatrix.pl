#!/usr/bin/env perl
# completeMatrix.pl
#
# Samuel Shepard - 2013
#
# Description: complete the lower diagonal to make a matrix square

use Getopt::Long;
GetOptions( 'annot-col|A' => \$hasAnnot );

if ( scalar(@ARGV) != 1 ) {
    die("Usage:\n\tperl $0 <lower_diagonal.txt> [-A|--annot-col]\n");
}

open( IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0]\n");
$/     = "\n";
@lines = <IN>;
chomp(@lines);
close(IN);

@h = @a = ();
@M = ();
$N = scalar(@lines);
for ( $i = 0; $i < $N; $i++ ) {
    $M[$i] = [split( "\t", $lines[$i] )];
    $h[$i] = shift( @{ $M[$i] } );
    if ( defined($hasAnnot) ) {
        $a[$i] = shift( @{ $M[$i] } );
    }
}

open( OUT, '>', $ARGV[0] ) or die("Cannot open $ARGV[1] for writing.\n");
for ( $i = 0; $i < $N; $i++ ) {
    print OUT $h[$i];
    if ( defined($hasAnnot) ) {
        print OUT "\t", $a[$i];
    }
    for ( $j = 0; $j < $N; $j++ ) {
        if ( $j <= $i ) {
            print OUT "\t", $M[$i][$j];
        } else {
            print OUT "\t", $M[$j][$i];
        }
    }
    print OUT "\n";
}
close(OUT);
