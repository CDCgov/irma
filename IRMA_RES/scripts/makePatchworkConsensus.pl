#!/usr/bin/env perl
# makePatchworkConsensus
#
# Sam Shepard - 2019
#
# Description: patches potentially missing data in a pHMM alignment for
# secondary assembly. See: https://wonder.cdc.gov/amd/flu/irma/secondary_residual.html

use strict;
use warnings;
use Getopt::Long;

my ( $name, $message );
GetOptions( 'name|N=s' => \$name );

if ( scalar(@ARGV) != 2 ) {
    $message = "Usage:\n\tperl $0 <a2m> <keyword_for_pairing>\n";
    die( $message . "\n" );
}

# PROCESS a2m data
my @count   = ();
my @default = ();
my @lines   = ();
my @col     = ();
my @nts     = ();
my ( $header, $sequence, $length, $nt, $p );
my $L       = 0;
my $keyword = $ARGV[1];

$/ = ">";
open( IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
while ( my $record = <IN> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = shift(@lines);
    $sequence = join( '', @lines );
    $length   = length($sequence);
    if ( $length == 0 ) { next; }

    @col = split( '', $sequence );
    if ( $header =~ /\Q$keyword\E/ ) {
        foreach $p ( 0 .. ( $length - 1 ) ) {
            $default[$p]{ $col[$p] }++;
        }

    } else {
        foreach $p ( 0 .. ( $length - 1 ) ) {
            $count[$p]{ $col[$p] }++;
        }
    }
}
close(IN);
$L = $length;

if ( defined($name) && $name ne '' ) {
    print STDOUT '>', $name, "\n";
} else {
    print STDOUT ">consensus\n";
}

$sequence = '';
foreach my $p ( 0 .. ( $length - 1 ) ) {
    @nts = sort { $count[$p]{$b} <=> $count[$p]{$a} } keys( %{ $count[$p] } );
    $nt  = defined( $nts[0] ) ? $nts[0] : '.';
    if ( $nt =~ /[.-]/ ) {
        @nts = sort { $default[$p]{$b} <=> $default[$p]{$a} } keys( %{ $default[$p] } );
        $nt  = defined( $nts[0] ) ? $nts[0] : '.';
    }
    $sequence .= uc($nt);
}
$sequence =~ tr/.-//d;
print STDOUT $sequence, "\n";
