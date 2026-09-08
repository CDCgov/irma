#!/usr/bin/env perl
# finalRefs.pl
#
# Sam Shepard - 2014
#
# Description: report the latest read-gathering reference round available for
# each gene.

use English qw(-no_match_vars);
use File::Basename;
use Getopt::Long;
GetOptions( 'ignore-annotation|G' => \$ignoreAnnotation );

my %maxRoundByGene = ();
local $RS = '>';
foreach my $file (@ARGV) {
    open( my $IN, '<', $file ) or die("Cannot open $file for reading.\n");
    my $round = basename( $file, '.refs' );
    if ( $round =~ /R(\d+)/ ) {
        $round = $1;
    }

    while ( my $record = <$IN> ) {
        chomp($record);
        my @lines    = split( /\r\n|\n|\r/, $record );
        my $gene     = shift(@lines);
        my $sequence = lc( join( '', @lines ) );

        if ( length $sequence <= 0 ) {
            next;
        }

        if ( $ignoreAnnotation && $gene =~ /^([^{]+)\{[^}]*}/ ) {
            $gene = $1;
        }

        if ( !defined $maxRoundByGene{$gene} || $maxRoundByGene{$gene} < $round ) {
            $maxRoundByGene{$gene} = $round;
        }
    }
    close($IN);
}

my @genes = sort( keys(%maxRoundByGene) );

if ( !@genes ) {
    die("No genes were found for final assembly!\n");
}

print 'R', $maxRoundByGene{ $genes[0] }, '-', $genes[0];
for my $i ( 1 .. $#genes ) {
    print ' R', $maxRoundByGene{ $genes[$i] }, '-', $genes[$i];
}
