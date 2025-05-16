#!/usr/bin/env perl
# Filename:         getPairingStats
#
# Author:           Samuel S. Shepard - 2022
#
# Description:      Converts `storable` Illumina pairing statistics to a table.

use strict;
use warnings;
use Carp    qw(croak);
use English qw(-no_match_vars);

if ( scalar(@ARGV) < 1 ) {
    die("Usage:\n\tperl $PROGRAM_NAME <stats1> <...>\n\n");
}

my @keys  = qw(dmv fmv tmv obs insObs insErr);
my %table = ();
foreach my $file (@ARGV) {
    open( my $STATS, '<', $file ) or croak("Could not open $file: $OS_ERROR\n");
    while ( my $line = <$STATS> ) {
        chomp($line);
        my ( $rn, $key, $value ) = split( "\t", $line );
        if ( defined $rn && defined $key ) {
            $table{$rn}{$key} += $value // 0;
        }
    }
    close($STATS) or croak("Could not close file: $OS_ERROR\n");
}

foreach my $rn ( keys %table ) {
    my $obs    = $table{$rn}{'obs'}    // 0;
    my $tmv    = $table{$rn}{'tmv'}    // 0;
    my $fmv    = $table{$rn}{'fmv'}    // 0;
    my $dmv    = $table{$rn}{'dmv'}    // 0;
    my $insObs = $table{$rn}{'insObs'} // 0;
    my $insErr = $table{$rn}{'insErr'} // 0;

    my $TMJ  = $obs - $fmv - $tmv;
    my $hObs = $TMJ + $tmv;

    if ( $obs > 0 ) {
        print STDOUT $rn, "\tObservations\t$obs\n";
        print STDOUT $rn, "\tExpectedErrorRate\t", ( $fmv / $obs ),        "\n";
        print STDOUT $rn, "\tMinimumExpectedVariation\t", ( $tmv / $obs ), "\n";
        print STDOUT $rn, "\tMinimumDeletionErrorRate\t", ( $dmv / $obs ), "\n";
    } else {
        print STDOUT $rn, "\tObservations\t$obs\n";
        print STDOUT $rn, "\tExpectedErrorRate\t0\n";
        print STDOUT $rn, "\tMinimumExpectedVariation\t0\n";
        print STDOUT $rn, "\tMinimumDeletionErrorRate\t0\n";
    }

    if ( $insObs > 0 ) {
        if ( $insObs == $insErr ) {
            $insObs++;
        }
        $insErr /= $insObs;
        print STDOUT $rn, "\tMinimumInsertionErrorRate\t$insErr\n";
    } else {
        print STDOUT $rn, "\tMinimumInsertionErrorRate\t0\n";
    }
}
