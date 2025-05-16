#!/usr/bin/env perl
#
# Filename:         annotateSubstitutions
# Description:      Annotations variant calling effects on translation (assumes
#                   already in-frame) and rewrites the table.
#
# Date dedicated:   2022-09-30
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
use English qw(-no_match_vars);
use Carp    qw(croak);

if ( scalar @ARGV != 2 ) {
    die("Usage:\n\tperl $PROGRAM_NAME <reference> <alleles.txt> [OPTIONS]\n");
}

my %gc = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    '---' => '-'     # Deletion
);

# physio-chemical factors
# Atchley et al. 2008
# "Solving the protein sequence metric problem."
# Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Epub 2005 Apr 25.
my %pcf = (
    'A' => [-0.59, -1.3,  -0.73,  1.57, -0.15],
    'C' => [-1.34,  0.47, -0.86, -1.02, -0.26],
    'D' => [ 1.05,  0.3,  -3.66, -0.26, -3.24],
    'E' => [ 1.36, -1.45,  1.48,  0.11, -0.84],
    'F' => [-1.01, -0.59,  1.89, -0.4,   0.41],
    'G' => [-0.38,  1.65,  1.33,  1.05,  2.06],
    'H' => [ 0.34, -0.42, -1.67, -1.47, -0.08],
    'I' => [-1.24, -0.55,  2.13,  0.39,  0.82],
    'K' => [ 1.83, -0.56,  0.53, -0.28,  1.65],
    'L' => [-1.02, -0.99, -1.51,  1.27, -0.91],
    'M' => [-0.66, -1.52,  2.22, -1.01,  1.21],
    'N' => [ 0.95,  0.83,  1.3,  -0.17,  0.93],
    'P' => [ 0.19,  2.08, -1.63,  0.42, -1.39],
    'Q' => [ 0.93, -0.18, -3.01, -0.5,  -1.85],
    'R' => [ 1.54, -0.06,  1.5,   0.44,  2.9],
    'S' => [-0.23,  1.4,  -4.76,  0.67, -2.65],
    'T' => [-0.03,  0.33,  2.21,  0.91,  1.31],
    'V' => [-1.34, -0.28, -0.54,  1.24, -1.26],
    'W' => [-0.6,   0.01,  0.67, -2.13, -0.18],
    'Y' => [ 0.26,  0.83,  3.1,  -0.84,  1.51],
    '-' => [ 0,     0,     0,     0,    0]
);

my @aa      = sort( keys(%pcf) );
my $SEQ_N   = scalar(@aa);
my %distMat = ();
foreach my $i ( 0 .. $SEQ_N - 1 ) {
    foreach my $j ( 0 .. $SEQ_N - 1 ) {
        my $dist = 0;
        foreach my $k ( 0 .. 4 ) {
            $dist += ( $pcf{ $aa[$i] }[$k] - $pcf{ $aa[$j] }[$k] )**2;
        }
        $dist = sqrt($dist);
        $distMat{ $aa[$i] }{ $aa[$j] } = $dist;
    }
}
my $maxDist = $distMat{'S'}{'Y'};
shift(@aa);    # remove gaps
open( my $IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
local $RS = ">";
my $first_record = <$IN>;
$first_record = <$IN>;
chomp($first_record);
close $IN or croak("Cannot close file: $OS_ERROR\n");

my @lines    = split( /\r\n|\n|\r/smx, $first_record );
my $header   = shift(@lines);
my $sequence = uc( join( q{}, @lines ) );
my $length   = length($sequence);
my @aaSeq    = ();

( $aaSeq[0], $aaSeq[1], $aaSeq[2] ) = ( q{}, q{}, q{} );
my @seq = split( q{}, $sequence );

## Perl does not have step by 3 for range operator
## no critic (ControlStructures::ProhibitCStyleForLoops)
for ( my $i = 0; $i < $length; $i += 3 ) {
    if ( ( $i + 0 + 2 ) < $length ) {
        my $codon1 = substr( $sequence, $i, 3 );
        $aaSeq[0] .= $gc{$codon1};
    }

    if ( ( $i + 1 + 2 ) < $length ) {
        my $codon2 = substr( $sequence, $i + 1, 3 );
        $aaSeq[1] .= $gc{$codon2};
    }

    if ( ( $i + 2 + 2 ) < $length ) {
        my $codon3 = substr( $sequence, $i + 2, 3 );
        $aaSeq[2] .= $gc{$codon3};
    }
}

my $sc1 = ( $aaSeq[0] =~ tr/*// );
my $sc2 = ( $aaSeq[1] =~ tr/*// );
my $sc3 = ( $aaSeq[2] =~ tr/*// );
my ( $aaRef, $offset );
if ( $sc1 <= $sc2 && $sc1 <= $sc3 ) {
    $aaRef  = $aaSeq[0];
    $offset = 0;
    print STDERR "Frame 1 chosen.\n";
} elsif ( $sc2 <= $sc1 && $sc2 <= $sc3 ) {
    $aaRef  = $aaSeq[1];
    $offset = 1;
    print STDERR "Frame 2 chosen.\n";
} else {
    $aaRef  = $aaSeq[2];
    $offset = 2;
    print STDERR "Frame 3 chosen.\n";
}
my $N   = $length - $offset - ( $length % 3 );
my $UB  = $N - 1;
my @ref = split( q{}, $aaRef );

print STDERR $UB, "\t", $offset, "\n";
open( my $TBL_IN, '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");
local $RS = "\n";
$header = <$TBL_IN>;
chomp($header);
print STDOUT $header, "\tCN\tCP\thasSubstitution\tMutation\tPCD\taaAllele\tSignificant\tNonsense\tDeletion\n";
while ( my $tbl_record = <$TBL_IN> ) {
    chomp($tbl_record);
    my @fields = split( "\t", $tbl_record );
    my ( $position,  $allele,   $type )      = ( $fields[1], $fields[2], $fields[10] );
    my ( $frequency, $pairedUB, $qualityUB ) = ( $fields[5], $fields[8], $fields[9] );

    $position--;
    my $CN = int( ( $position - $offset ) / 3 );
    my $CP = ( $position - $offset ) % 3;

    my $aaAllele = $ref[$CN];

    my $sig = 'No';
    if ( $frequency > $pairedUB && $frequency > $qualityUB ) {
        $sig = 'Yes';
    }

    my $pcd      = 0;
    my $sub      = 'No';
    my $nonsense = 'No';
    my $mutation = q{};
    my $deletion = 'No';
    if ( $type eq 'Majority' ) {
        $mutation = $ref[$CN];
    } elsif ( $position < $offset || $position > $UB ) {
        $sub = 'OOB';
    } else {
        my $codon = substr( $sequence, $position - $CP, 3 );
        substr( $codon, $CP, 1, $allele );
        my $aa = defined $gc{$codon} ? $gc{$codon} : 'X';

        # Partial deletion
        if ( $codon !~ /^[A-Za-z]{3}$/smx ) {
            $aa = '~';
        }

        if ( $aa ne $ref[$CN] ) {
            $pcd      = defined $distMat{$aa}{ $ref[$CN] } ? $distMat{$aa}{ $ref[$CN] } : 0;
            $mutation = $ref[$CN] . ( $CN + 1 ) . $aa;

            if ( $aa eq '-' || $aa eq '~' ) {
                $deletion = 'Yes';
            } else {
                $sub = 'Yes';
                if ( $aa eq '*' ) {
                    $nonsense = 'Yes';
                }
            }
        }
        $aaAllele = $aa;
    }

    print STDOUT $tbl_record, "\t", ( $CN + 1 ), "\t", ( $CP + 1 ), "\t", $sub, "\t", $mutation, "\t", $pcd, "\t",
      $aaAllele, "\t", $sig, "\t", $nonsense, "\t", $deletion, "\n";
}
close $TBL_IN or croak("Cannot close file: $OS_ERROR\n");
