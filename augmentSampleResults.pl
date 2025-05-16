#!/usr/bin/env perl
#
# Filename:         augmentSampleResuults
# Description:      Globs, concatentates, and adds columns to IRMA output tables.
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
use File::Basename;
use Getopt::Long;

my ( $nameField, $field );
my ( $idOnly, $ignoreIDfield, $noHeader, $noKeyHeader, $timeStamp, $exactFilePath ) = ( 0, 0, 0, 0, 0, 0 );

my $word   = 'all';
my $infix  = 'tables';
my $suffix = q{};

GetOptions(
            'word|W=s'          => \$word,
            'suffix|S=s'        => \$suffix,
            'ID-only|I'         => \$idOnly,
            'infix|X=s'         => \$infix,
            'dir-field|F=i'     => \$field,
            'ignore-dir|G'      => \$ignoreIDfield,
            "no-header|H"       => \$noHeader,
            'no-key-header|K'   => \$noKeyHeader,
            'timestamp|T'       => \$timeStamp,
            'exact-file-path|E' => \$exactFilePath,
            'name|N=s'          => \$nameField
);

if ( scalar @ARGV != 1 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <key.txt> [options]\n"
         . "\t\t-H|--no-header\t\tDo not output header.\n"
         . "\t\t-K|--no-key-header\tDoes not contain a key header (implies -H).\n"
         . "\t\t-W|--word <STR>]\tContains <word> in table filename.\n"
         . "\t\t-S|--suffix <STR>\tTable filename ends in suffix.\n"
         . "\t\t-I|--ID-only\t\tKey just contains the ID.\n"
         . "\t\t-X|--infix <STR>\tInfix is the intermediate path. Default is 'tables/'\n"
         . "\t\t-F|--dir-field <INT>\tField containing dirname. Default = 2\n"
         . "\t\t-G|--ignore-dir\t\tDo not print out dirname into collated data.\n"
         . "\t\t-T|--timestamp\t\tAdd timestamp to the output as the last column.\n"
         . "\t\t-E|--exact-file-path\tExact file path is contained, nullifies infix & suffix.\n"
         . "\t\t-N|--name <STR>\t\tAdd a fixed name field between the key and table.\n"
         . "\n" );
}

my ( $now, $nowH ) = ( q{}, q{} );
if ($timeStamp) {
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime(time);
    $now  = sprintf( "\t%04d-%02d-%02d %02d:%02d:%02d", $year + 1900, $mon + 1, $mday, $hour, $min, $sec );
    $nowH = "\ttimestamp";
}

my $sampleDirField = 2;
if ( defined $field && $field > 0 ) {
    $sampleDirField = $field;
}

my ( $nameHdr, $nameVal ) = ( q{}, q{} );
if ( defined $nameField ) {
    $nameHdr = "\tName";
    $nameVal = "\t$nameField";
}

open( my $IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0].\n");
my $hdrs;
my @tmp = ();
if ($noKeyHeader) {
    $hdrs     = q{};
    $noHeader = 1;
} else {
    $hdrs = <$IN>;
    chomp($hdrs);
    @tmp = split( "\t", $hdrs );
    if ($idOnly) {
        $hdrs = $tmp[0];
    }
}

my $first = 1;
if ($ignoreIDfield) {
    my $line = q{};
    $first = 1;
    foreach my $i ( 0 .. $#tmp ) {
        if ( $i != ( $sampleDirField - 1 ) ) {
            if ($first) {
                $first = 0;
                $line  = $tmp[$i];
            } else {
                $line .= "\t" . $tmp[$i];
            }
        }
    }
    $hdrs = $line;
}

my $firstData = 1;
if ($noHeader) {
    $firstData = 0;
}

while ( my $line = <$IN> ) {
    chomp($line);
    my @fields = split( "\t", $line );
    if ($idOnly) { $line = $fields[0]; }

    my @files = ();
    my $id    = q{};
    if ($exactFilePath) {
        @files = ( $fields[$sampleDirField - 1] );
    } else {
        $id    = $fields[$sampleDirField - 1];
        @files = glob("$id/$infix/*$word*$suffix");
    }

    if ($ignoreIDfield) {
        $line  = q{};
        $first = 1;
        foreach my $i ( 0 .. $#fields ) {
            if ( $i != ( $sampleDirField - 1 ) ) {
                if ($first) {
                    $first = 0;
                    $line  = $fields[$i];
                } else {
                    $line .= "\t" . $fields[$i];
                }
            }
        }
    }

    foreach my $file (@files) {
        open( my $DAT, '<', $file ) or die("Cannot open $file for reading.\n");
        my $header = <$DAT>;
        chomp($header);
        my $numCols  = scalar( split( "\t", $header ) );
        my $basename = basename($file);

        if ($firstData) {
            print $hdrs, $nameHdr, "\t", $header, $nowH, "\n";
            $firstData = 0;
        }

        while ( my $line2 = <$DAT> ) {
            chomp($line2);
            print STDOUT $line, $nameVal, "\t", $line2, $now, "\n";
        }

        close $DAT or croak("Cannot close file: $OS_ERROR\n");
    }
}
close $IN or croak("Cannot close file: $OS_ERROR\n");
