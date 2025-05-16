#!/usr/bin/env perl
# makePluralityConsensus
#
# Sam Shepard - 2016
#
# Description: make a simple plurality consensus for use in module creation (reference seeds)

use Getopt::Long;
GetOptions( 'name|N=s' => \$name, 'allow-deletions|D' => \$allowDeletions );

if ( scalar(@ARGV) != 1 ) {
    $message = "Usage:\n\tperl $0 <input.fasta> [-N|--name <STRING>] [-D|--allow-deletions]\n";
    die( $message . "\n" );
}

# PROCESS fasta data
open( IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
$/     = ">";
%count = ();
$L     = 0;
while ( $record = <IN> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = uc( join( '', @lines ) );
    $length   = length($sequence);

    if ( $length == 0 ) {
        next;
    }

    # unify gap characters
    $sequence =~ tr/:~./-/;
    @col = split( '', $sequence );
    foreach $p ( 0 .. ( $length - 1 ) ) {
        $count[$p]{ $col[$p] }++;
    }
}
close(IN);
$L = $length;

if ($name) {
    print '>', $name, "\n";
} else {
    print ">consensus\n";
}

foreach $p ( 0 .. ( $length - 1 ) ) {
    @nts = sort { $count[$p]{$b} <=> $count[$p]{$a} } keys( %{ $count[$p] } );
    foreach $nt ( 0 .. $#nts ) {
        if ( $nts[$nt] ne '-' || $allowDeletions ) {
            print $nts[$nt];
            last;
        }
    }
}
print "\n";
