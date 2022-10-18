#!/usr/bin/env perl
#
# Filename:         alignPositionByPSA
# Description:      Re-aligns tables using a pairwise sequence alignment to a
#                   reference.
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
use Carp qw(croak);
use Getopt::Long;

my ( $skipHeader, $insertsToRef, $showOriginal, $skipComments, $reprintComments, $checkComments ) = ( 0, 0, 0, 0, 0, 0 );
my ( $field, $delim, $name, $prefix );
GetOptions(
            'no-header|H'              => \$skipHeader,
            'table-field|F=i'          => \$field,
            'table-delim|D=s'          => \$delim,
            'ref-name|N=s'             => \$name,
            'prefix|P=s'               => \$prefix,
            'inserts-to-ref|I'         => \$insertsToRef,
            'show-original-position|O' => \$showOriginal,
            'skip-comments|S'          => \$skipComments,
            'reprint-comments|R'       => \$reprintComments
);

if ( scalar(@ARGV) != 2 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <pairwiseAlignment.fasta> <table> [options]\n"
         . "\t\t-H|--no-header\t\t\tDo not print input table header line. (first line)\n"
         . "\t\t-F|--table-field <INT>\t\tField number in table (1-base) containing positional information. DEFAULT = 1\n"
         . "\t\t-D|--table-delim <CHA>\t\tField delimiter for table. DEFAULT = <TAB>\n"
         . "\t\t-N|--ref-name <STR>\t\tName of the reference sequence. DEFAULT = first record\n"
         . "\t\t-P|--prefix <STR>\t\tPrefix to table.\n"
         . "\t\t-I|--inserts-to-ref\t\tInserts relative to reference.\n"
         . "\t\t-O|--show-original-position\tShow original position column.\n"
         . "\t\t-S|--skip-comments\t\tSkip lines beginning with #.\n"
         . "\t\t-R|--reprint-comments\t\tSimply reprint lines beginning with #.\n"
         . "\n" );
}

if ( $skipComments || $reprintComments ) {
    $checkComments = 1;
}

if ( !defined $delim ) {
    $delim = "\t";
}

if ( !defined $field ) {
    $field = 0;
} elsif ( $field > 0 ) {
    $field = int($field) - 1;
} else {
    $field = 0;
}

my $firstRecord = 0;
if ( !defined $name ) {
    $firstRecord = 1;
}

if ( defined $prefix ) {
    $prefix = $prefix . "\t";
} else {
    $prefix = q{};
}

local $RS = ">";
open( my $FASTA_IN, '<', $ARGV[0] ) or die("ERROR: Cannot open $ARGV[0] for reading.\n");
my $K          = 0;
my %matchByRef = ();
my ( $refSeq, $refHeader, $refLength );
my ( $altSeq, $altHeader, $altLength );
while ( my $fasta_record = <$FASTA_IN> ) {
    chomp($fasta_record);
    my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
    my $header   = shift(@lines);
    my $sequence = uc( join( q{}, @lines ) );
    my $length   = length($sequence);
    if ( $length == 0 ) {
        next;
    } else {
        $K++;
    }

    if ( ( $firstRecord && $K == 1 ) || ( !$firstRecord && $header =~ /$name/smx ) ) {
        $refSeq    = $sequence;
        $refHeader = $header;
        $refLength = $length;
    } else {
        $altSeq    = $sequence;
        $altHeader = $header;
        $altLength = $length;
    }
}
close $FASTA_IN or croak("Cannot close file: $OS_ERROR\n");

if ( $refLength != $altLength ) {
    die("ERROR $PROGRAM_NAME: Unequal lengths ($refLength <> $altLength).\n");
}

my @ref      = split( q{}, $refSeq );
my @alt      = split( q{}, $altSeq );
my %refByAlt = ();
my ( $aCoord, $rCoord, $insertIndex ) = ( 0, 0, 0 );
foreach my $i ( 0 .. $refLength - 1 ) {
    my $r = $ref[$i];
    my $a = $alt[$i];

    if ( $r ne '-' ) {
        $insertIndex = 0;
        $rCoord++;
    }

    if ( $a ne '-' ) {
        $aCoord++;
    }

    # MATCH
    if ( $a ne '-' && $r ne '-' ) {
        $refByAlt{$aCoord} = $rCoord;

        # INSERTION relative to reference
    } elsif ( $a ne '-' ) {
        $insertIndex++;
        $refByAlt{$aCoord} = sprintf( '%d.%04d', $rCoord, $insertIndex );
    }
}

local $RS = "\n";
open( my $TBL_IN, '<', $ARGV[1] ) or die("ERROR: Cannot open $ARGV[1] for reading.\n");
my $header = <$TBL_IN>;
if ( $checkComments && substr( $header, 0, 1 ) eq '#' ) {
    if ($reprintComments) {
        print STDOUT $header;
    }
} else {
    if ( !$skipHeader ) {
        if ($showOriginal) {
            chomp($header);
            my @H = split( $delim, $header );
            $H[$field] = $H[$field] . '_original' . $delim . $H[$field] . '_revised';

            print STDOUT $prefix, join( $delim, @H ), "\n";
        } else {
            print STDOUT $header;
        }
    }
}
my @data = <$TBL_IN>;
chomp(@data);
close $TBL_IN or croak("Cannot close file: $OS_ERROR\n");

my $removeInserts = 1;
if ($insertsToRef) {
    $removeInserts = 0;
}

( $aCoord, $rCoord ) = ( 0, 0 );
foreach my $line (@data) {
    if ( $checkComments && substr( $line, 0, 1 ) eq '#' ) {
        if ($reprintComments) {
            print $line, "\n";
        }
        next;
    }
    my @fields = split( $delim, $line );
    $aCoord = $fields[$field];

    if ( !defined $refByAlt{$aCoord} ) {
        print STDERR $aCoord, " not found\n";
        next;
    } else {
        $rCoord = $refByAlt{$aCoord};
    }

    if ( $removeInserts && $rCoord =~ /\./smx ) { next; }

    $fields[$field] = $showOriginal ? $aCoord . $delim . $rCoord : $rCoord;
    print STDOUT $prefix, join( $delim, @fields ), "\n";
}
