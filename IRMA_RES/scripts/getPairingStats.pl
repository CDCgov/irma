#!/usr/bin/env perl

# Filename:         getPairingStats
# Description:      Converts `storable` Illumina pairing statistics to a table.
#
# Date dedicated:   2022-10-21
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Shepard SS, Meno S, Bahl J, Wilson MM, Barnes J, Neuhaus E.
#                   Viral deep sequencing needs an adaptive approach: IRMA, the
#                   iterative refinement meta-assembler. BMC Genomics.
#                   2016;17(1). doi:10.1186/s12864-016-3030-6
#
# =============================================================================
#
#                            PUBLIC DOMAIN NOTICE
#
#  This source code file or script constitutes a work of the United States
#  Government and is not subject to domestic copyright protection under 17 USC §
#  105. This file is in the public domain within the United States, and
#  copyright and related rights in the work worldwide are waived through the CC0
#  1.0 Universal public domain dedication:
#  https://creativecommons.org/publicdomain/zero/1.0/
#
#  The material embodied in this software is provided to you "as-is" and without
#  warranty of any kind, express, implied or otherwise, including without
#  limitation, any warranty of fitness for a particular purpose. In no event
#  shall the Centers for Disease Control and Prevention (CDC) or the United
#  States (U.S.) government be liable to you or anyone else for any direct,
#  special, incidental, indirect or consequential damages of any kind, or any
#  damages whatsoever, including without limitation, loss of profit, loss of
#  use, savings or revenue, or the claims of third parties, whether or not CDC
#  or the U.S. government has been advised of the possibility of such loss,
#  however caused and on any theory of liability, arising out of or in
#  connection with the possession, use or performance of this software.
#
#  Please provide appropriate attribution in any work or product based on this
#  material.

use strict;
use warnings;
use Carp qw(croak);
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
