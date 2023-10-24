#!/usr/bin/env perl

# Filename:         fastQ_converter
# Description:      Read FastQ files, applies QC filtering (quality and length),
#                   adapter trimming, and format conversion as requested.
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

use warnings;
use strict;

use English qw(-no_match_vars);
use Carp qw(croak);
use Getopt::Long;

# "$baseQuality" masking not implemented and removed.
my $readSide;

## Heuristics
my $qualityThreshold = 0;
my $minLength        = 0;

## IO Parameters
my $fileID          = q{};
my $saveQualityFile = q{};
my $saveStats       = q{};
my $logFile         = q{};
my $logID           = q{};

## All Flags
my $useMedian        = 0;
my $fastQformat      = 0;
my $complementAndAdd = 0;
my $ordinal          = 0;
my $skipRemaining    = 0;    # Skips the rest of the loop when using saveStats
my $keepHeader       = 0;
my $maskAdapter      = 0;
my $clipAdapter      = 0;
my $uracilToThymine  = 0;
my $clippedMinLength = 0;
my $fuzzyAdapter     = 0;

Getopt::Long::Configure('no_ignore_case');
GetOptions(
            'read-quality|T=i'         => \$qualityThreshold,
            'use-median|M'             => \$useMedian,
            'fastQ-output|Q'           => \$fastQformat,
            'min-length|L=i'           => \$minLength,
            'complement-add|C'         => \$complementAndAdd,
            'ordinal-headers|O'        => \$ordinal,
            'file-id|F=s'              => \$fileID,
            'save-quality|S=s'         => \$saveQualityFile,
            'save-stats|A=s'           => \$saveStats,
            'skip-remaining|K'         => \$skipRemaining,
            'log-file|G=s'             => \$logFile,
            'log-id|g=s'               => \$logID,
            'keep-header|H'            => \$keepHeader,
            'mask-adapter|m=s'         => \$maskAdapter,
            'clip-adapter|c=s'         => \$clipAdapter,
            'fuzzy-adapter|Z'          => \$fuzzyAdapter,
            'uracil-to-thymine|U'      => \$uracilToThymine,
            'enforce-clipped-length|E' => \$clippedMinLength,
            'read-side|R=i'            => \$readSide
);

if ( -t STDIN && scalar @ARGV != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <fastQ> [options]\n"
         . "\t\t-T|--read-quality <threshold>\t\tSpecify the read quality threshold (geometric mean, median).\n"
         . "\t\t-M|--use-median\t\t\t\tInterprets the threshold (-T) as the median, not average.\n"
         . "\t\t-Q|--fastQ-output\t\t\tOutputs fastQ instead of fastA format.\n"
         . "\t\t-L|--min-length <threshold>\t\tMinimum length of sequence data, default = 0.\n"
         . "\t\t-C|--complement-add\t\t\tTake the reverse complement and add to data.\n"
         . "\t\t-O|--ordinal-headers\t\t\tReplace header with strict ordinals.\n"
         . "\t\t-F|--file-id <STR>\t\t\tFile id for ordinals.\n"
         . "\t\t-S|--save-quality <STR>\t\t\tSave quality file for back-mapping.\n"
         . "\t\t-A|--save-stats <STR>\t\t\tSave quality vs. length statistics file for analysis.\n"
         . "\t\t-K|--skip-remaining\t\t\tDo not output data FASTA/FASTQ data (assumes -A).\n"
         . "\t\t-H|--keep-header\t\t\tKeep header as usual.\n"
         . "\t\t-c|--clip-adapter <STR>\t\t\tClip adapter.\n"
         . "\t\t-m|--mask-adapter <STR>\t\t\tMask adapter.\n"
         . "\t\t-Z|--fuzzy-adapter\t\t\tAllow one mismatch.\n"
         . "\t\t-U|--uracil-to-thymine\t\t\tCovert uracil to thymine.\n"
         . "\t\t-E|--enforce-clipped-length\t\tThe minimum length threshold (-L) is enforced when adapter clipped (-c).\n"
         . "\t\t-R|--read-side <INT>\t\t\tIf FASTQ header is SRA format and missing a read identifier, alter the header."
         . "\n" );
}

my $LOG_OUT;
if ( $logFile ne q{} ) {
    if ( $logID ne q{} ) {
        $logID = ':' . $logID;
    }
    open( $LOG_OUT, '>>', $logFile ) or die("Cannot open $logFile for appending.\n");
}

if ( $fileID ne q{} ) {
    if ($complementAndAdd) {
        $fileID = '|' . $fileID . '|';
    } else {
        $fileID = $fileID . '|';
    }
}

my $STAT;
if ( $saveStats ne q{} ) {
    open( $STAT, '>', $saveStats ) or die("Cannot open $saveStats for writing.\n");
    $saveStats = 1;
}

my $QUA;
if ( $saveQualityFile ne q{} ) {
    open( $QUA, '>', $saveQualityFile ) or die("Cannot open $saveQualityFile.\n");
    $saveQualityFile = 1;
}

my $forwardAdapter = q{};
my $reverseAdapter = q{};
my $adapterMask    = q{};

if ($clipAdapter) {
    $forwardAdapter = $clipAdapter;
    $reverseAdapter = reverse($forwardAdapter);
    $reverseAdapter =~ tr/ATGC/TACG/;
}

if ($maskAdapter) {
    $forwardAdapter = $maskAdapter;
    $reverseAdapter = reverse($forwardAdapter);
    $reverseAdapter =~ tr/ATGC/TACG/;
    $adapterMask = 'N' x length($forwardAdapter);
}

my $fuzzyAdaptersFWD = q{};
my $fuzzyAdaptersREV = q{};
if ( $fuzzyAdapter && ( $maskAdapter || $clipAdapter ) ) {
    my $L = length $forwardAdapter;
    if ( $L > 0 ) {
        my $N = '[ATCGN]';
        $L--;
        foreach my $i ( 0 .. $L ) {
            my $tmp = $forwardAdapter;
            substr( $tmp, $i, 1, $N );
            $fuzzyAdaptersFWD .= $tmp . '|';
            $tmp = $reverseAdapter;
            substr( $tmp, $i, 1, $N );
            $fuzzyAdaptersREV .= $tmp . '|';
        }

        chop($fuzzyAdaptersFWD);
        chop($fuzzyAdaptersREV);
        $fuzzyAdaptersFWD = qr/$fuzzyAdaptersFWD/ismx;
        $fuzzyAdaptersREV = qr/$fuzzyAdaptersREV/ismx;
    }
} else {
    $fuzzyAdapter = 0;
}

if ( $clipAdapter || $maskAdapter ) {
    $forwardAdapter = qr/$forwardAdapter/ismx;
    $reverseAdapter = qr/$reverseAdapter/ismx;
}

if ( $minLength < 0 ) {
    die("ERROR: minimum length must be a non-negative integer.\n");
} else {
    $minLength = int($minLength);
}

if ( defined $readSide && ( $readSide < 0 || $readSide > 3 ) ) {
    die("ERROR: --read-side must be in range [0, 3]\n");
}

local $RS = "\n";

my $reads_passing_qc = 0;
my $id               = 0;
my $pp               = q{};
my $nn               = q{};
my $dnp              = q{};

if ($complementAndAdd) {
    $pp  = 'P';    # position strand
    $nn  = 'N';    # negative strand
    $dnp = '_';    # delimiter
}

my $give_warning_for_long_fastq = 1;

sub fix_SRA_format {
    my $s         = shift // q{};
    my $read_side = shift // '0';
    my $delim     = q{ };
    my $pre       = substr( $s, 0, 4 );

    if ( index( $s, q{ } ) == -1 ) {
        $delim = '_';
    }

    my ( $mol_id, @remainder ) = split $delim, $s;
    $mol_id //= q{};

    if ( ( $pre eq '@SRR' || $pre eq '@DRR' || $pre eq '@ERR' ) && index( $mol_id, '.' ) != -1 ) {

        # SRA format
        my @fields = split '\.', $mol_id;
        if ( scalar @fields == 2 ) {
            return join( $delim, ( $mol_id . '.' . $read_side, @remainder ) );
        }
    }

    return $s;
}

while ( my $hdr = <> ) {
    chomp($hdr);

    if ( defined $readSide ) {
        $hdr = fix_SRA_format( $hdr, $readSide );
    }

    my $seq = <>;
    chomp($seq);

    my $junk = <>;
    chomp($junk);

    my $quality = <>;
    chomp($quality);

    if ( length $seq < $minLength ) {
        next;
    }

    if ($uracilToThymine) {
        $seq =~ tr/uU/tT/;
    }

    if ($clipAdapter) {

        # Required for backwards compatible behavior
        ## no critic (ControlStructures::ProhibitCascadingIfElse)
        if ( $seq =~ $reverseAdapter ) {

            # Trim 3' end
            $seq     = substr( $seq,     0, $LAST_MATCH_START[0] );
            $quality = substr( $quality, 0, $LAST_MATCH_START[0] );
        } elsif ( $seq =~ $forwardAdapter ) {

            # Trim 5' end
            $seq     = substr( $seq,     $LAST_MATCH_END[0] );
            $quality = substr( $quality, $LAST_MATCH_END[0] );
        } elsif ( $fuzzyAdapter && $seq =~ $fuzzyAdaptersREV ) {

            # Trim 3' end
            $seq     = substr( $seq,     0, $LAST_MATCH_START[0] );
            $quality = substr( $quality, 0, $LAST_MATCH_START[0] );
        } elsif ( $fuzzyAdapter && $seq =~ $fuzzyAdaptersFWD ) {

            # Trim 5' end
            $seq     = substr( $seq,     $LAST_MATCH_END[0] );
            $quality = substr( $quality, $LAST_MATCH_END[0] );
        }

        if ( $clippedMinLength && length $seq < $minLength ) {
            next;
        }
    } elsif ($maskAdapter) {

        # Logic short circuits to prevent forward search if reverse found.
        # Result is saved to aid readability.
        my $wasMasked = ( $seq =~ s/$reverseAdapter/$adapterMask/ismx or $seq =~ s/$forwardAdapter/$adapterMask/ismx );
        if ( $fuzzyAdapter && !$wasMasked ) {
            $seq =~ s/$fuzzyAdaptersREV/$adapterMask/ismx or $seq =~ s/$fuzzyAdaptersFWD/$adapterMask/ismx;
        }
    }

    my @a = unpack( "c* i*", $quality );
    my $q = 0;
    my $n = scalar(@a);
    $id++;

    if ($useMedian) {
        my @sorted = sort(@a);
        if ( $n % 2 == 0 ) {
            $q = ( $sorted[$n / 2] + $sorted[( $n / 2 ) - 1] ) / 2 - 33;
        } else {
            $q = $sorted[( $n - 1 ) / 2] - 33;
        }
    } else {
        ## use geometric mean error otherwise (average QS)
        foreach my $x (@a) {
            $q += $x;
        }
        $q = ( $q - $n * 33 ) / $n;
    }

    if ($saveStats) {
        print $STAT $id, "\t", $q, "\t", $n, "\n";
        if ($skipRemaining) {
            next;
        }
    }

    if ( $q >= $qualityThreshold ) {
        $reads_passing_qc++;
        $hdr = substr( $hdr, 1 );

        if ( !$keepHeader ) {
            if ( length $hdr > 254 ) {

                # Default is to truncate bytes. If a unicode character is on the
                # truncation boundary, could cause issues for UTF-8 readers.
                $hdr = substr( $hdr, 0, 254 );
                if ($give_warning_for_long_fastq) {
                    $give_warning_for_long_fastq = 0;
                    print STDERR "$PROGRAM_NAME WARNING: FASTQ headers truncated,",
                      " downstream BAM format expects no more than 254 bytes!\n";
                }
            }
            $hdr =~ tr/ /_/;
        }

        if ($ordinal) {
            if ($fastQformat) {
                print STDOUT '@', $pp, $fileID, $id, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
            } else {
                print STDOUT '>', $pp, $fileID, $id, "\n", $seq, "\n";
                if ($saveQualityFile) {
                    print $QUA $pp, $fileID, $id, "\t", $hdr, "\t", $quality, "\n";
                }
            }
        } else {
            if ($fastQformat) {
                print STDOUT '@', $hdr, $dnp, $pp, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
            } else {
                print STDOUT '>', $hdr, $dnp, $pp, '|', sprintf( "%.1f", $q ), '|', $n, "\n", $seq, "\n";
                if ($saveQualityFile) {
                    print $QUA $hdr, $dnp, $pp, "\t", $quality, "\n";
                }
            }
        }

        # Take the reverse complement and add it
        # courtesy http://reverse-complement.com/ambiguity.html
        if ($complementAndAdd) {
            $seq = reverse($seq);
            $seq =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
            $quality = reverse($quality);
            if ($ordinal) {
                if ($fastQformat) {
                    print STDOUT '@', $nn, $fileID, $id, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
                } else {
                    print STDOUT '>', $nn, $fileID, $id, "\n", $seq, "\n";
                    if ($saveQualityFile) {
                        print $QUA $nn, $fileID, $id, "\t", $hdr, "\t", $quality, "\n";
                    }
                }
            } else {
                if ($fastQformat) {
                    print STDOUT '@', $hdr, $dnp, $nn, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
                } else {
                    print STDOUT '>', $hdr, $dnp, $nn, '|', sprintf( "%.1f", $q ), '|', $n, "\n", $seq, "\n";
                    if ($saveQualityFile) {
                        print $QUA $hdr, $dnp, $nn, "\t", $quality, "\n";
                    }
                }
            }
        }
    }
}

if ( $logFile ne q{} ) {
    my $reads_passing_length = $id;
    print $LOG_OUT $logFile, "$logID\t", $reads_passing_length, "\t", $reads_passing_qc, "\t", $qualityThreshold, "\t",
      $minLength,
      "\t", $useMedian, "\n";
    close $LOG_OUT or croak("Cannot close file: $OS_ERROR\n");
}

if ($saveStats) {
    close $STAT or croak("Cannot close file: $OS_ERROR\n");
}

if ($saveQualityFile) {
    close $QUA or croak("Cannot close file: $OS_ERROR\n");
}
