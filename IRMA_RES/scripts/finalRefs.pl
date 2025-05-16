#!/usr/bin/env perl
# finalRefs.pl
#
# Sam Shepard - 2014
#
# Description: get rid of alternative references after read gathering (if applicable)

use File::Basename;
use Getopt::Long;
GetOptions( 'ignore-annotation|G'   => \$ignoreAnnotation,
            'exclude-alternative|X' => \$excludeAlt );

$/ = '>';
foreach $file (@ARGV) {
    open( IN, '<', $file ) or die("Cannot open $file for reading.\n");
    $round = basename( $file, '.refs' );
    if ( $round =~ /R(\d+)/ ) {
        $round = $1;
    }

    while ( $record = <IN> ) {
        chomp($record);
        @lines    = split( /\r\n|\n|\r/, $record );
        $gene     = shift(@lines);
        $sequence = lc( join( '', @lines ) );

        if ( length($sequence) <= 0 ) {
            next;
        }

        if ( $excludeAlt && $gene =~ /{alt}/ ) {
            next;
        }

        if ( $ignoreAnnotation && $gene =~ /^([^{]+)\{[^}]*}/ ) {
            $gene = $1;
        }

        if ( !defined($maxRoundByGene) || $maxRoundByGene{$gene} < $round ) {
            $seqByGene{$gene}      = $sequence;
            $maxRoundByGene{$gene} = $round;
        }
    }
    close(IN);
}

@genes = sort( keys(%seqByGene) );
print 'R', $maxRoundByGene{ $genes[0] }, '-', $genes[0];
for ( $i = 1; $i < scalar(@genes); $i++ ) {
    print ' R', $maxRoundByGene{ $genes[$i] }, '-', $genes[$i];
}
