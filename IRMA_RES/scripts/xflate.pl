#!/usr/bin/env perl
#
# Filename:         xflate
# Description:      Turns reads into deflated read patterns and inflates them
# 					back to reads.
#
# Date dedicated:   2023-07-24
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

my ( $inflate, $clusterAll, $fastQ, $label, $rci );

use Carp qw(croak);
use English qw(-no_match_vars);
use Getopt::Long;
GetOptions(
            'inflate|I'                    => \$inflate,
            'cluster-all|C'                => \$clusterAll,
            'file-label|L=s'               => \$label,
            'reverse-complement-inflate|R' => \$rci,
            'fastQ|Q'                      => \$fastQ
);

if ( scalar(@ARGV) < 2 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <xflate.txt> [options] <fasta1.fa> ...\n"
         . "\t\t-I|--inflate\t\tInflate sequence files.\n"
         . "\t\t-C|--cluster-all\tUse a cluster for all sequences.\n"
         . "\t\t-L|--file-label <STR>\tLabels the clusters with an identifier.\n"
         . "\t\t-R|--rev-comp-inf\tReverse complement inflation.\n"
         . "\t\t-Q|--fastQ\t\tExpects fastQ input (deflate) and fastQ output (inflate).\n"
         . "\n" );
}

if ( !defined($label) ) {
    $label = q{};
} else {
    $label = '|' . $label;
}

my $tableFile = shift(@ARGV);

## INFLATE
if ($inflate) {
    open( my $TABLE, '<', $tableFile ) or die("Cannot open $tableFile for reading (trying to inflate without a file?).\n");

    %sequenceByCluster = ();
    foreach my $file (@ARGV) {
        local $RS = '>';
        open( my $IN, '<', $file ) or die("Cannot open $file.\n");
        while ( my $fasta_record = <$IN> ) {
            chomp($fasta_record);
            my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
            my $header   = shift(@lines);
            my $sequence = lc( join( q{}, @lines ) );

            my $length = length($sequence);
            if ( $length == 0 ) {
                next;
            }

            if ( $header =~ /^(C\d+%\d+%[^{]*)/smx ) {
                my $cluster = $1;
                if ( $rci && $header =~ /{c}/smx ) {
                    $sequence = reverse($sequence);
                    $sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
                }

                $sequenceByCluster{$cluster} = $sequence;
            } else {
                if ($fastQ) {
                    my $quality = '?' x length $sequence;
                    print STDOUT '@', $header, "\n", $sequence, "\n+\n", $quality, "\n";
                } else {
                    print STDOUT '>', $header, "\n", $sequence, "\n";
                }
            }
        }
    }

    local $RS = "\n";
    if ($fastQ) {
        while ( my $line = <$TABLE> ) {
            chomp($line);
            my @fields  = split( "\t", $line );
            my $cluster = $fields[0];

            if ( defined $sequenceByCluster{$cluster} ) {
                ## Must increment 2 at a time.
                ## no critic (ControlStructures::ProhibitCStyleForLoops)
                for ( my $i = 1; $i < scalar(@fields); $i += 2 ) {
                    print STDOUT $fields[$i], "\n", $sequenceByCluster{$cluster}, "\n+\n", $fields[$i + 1], "\n";
                }
            }
        }
    } else {
        while ( my $line = <$TABLE> ) {
            chomp($line);
            my @fields  = split( "\t", $line );
            my $cluster = $fields[0];

            if ( defined $sequenceByCluster{$cluster} ) {
                ## Must increment 2 at a time.
                ## no critic (ControlStructures::ProhibitCStyleForLoops)
                for ( my $i = 1; $i < scalar(@fields); $i += 2 ) {
                    print STDOUT '>', $fields[$i], "\n", $sequenceByCluster{$cluster}, "\n";
                }
            }
        }
    }
    close $TABLE or croak("Cannot close file: $OS_ERROR");

## DEFLATE
} else {
    if ( -e $tableFile ) {
        print STDERR "WARNING, using ${tableFile}2 since file exists.\n";
        $tableFile .= '2';
    }
    open( my $TABLE, '>', $tableFile ) or die("Cannot open $tableFile for reading (trying to inflate without a file?).\n");

    # Adding a 'my' declaration adds a severe performance penalty for some reason.
    %clustersBySequence = ();
    %qualityByHeader    = ();

    foreach my $file (@ARGV) {
        open( my $IN, '<', $file ) or die("Cannot open $file.\n");
        if ($fastQ) {
            local $RS = "\n";
            while ( my $header = <$IN> ) {
                chomp($header);
                my $sequence = <$IN>;
                chomp($sequence);
                my $junk = <$IN>;
                chomp($junk);
                my $quality = <$IN>;
                chomp($quality);

                if ( !defined $clustersBySequence{$sequence} ) {
                    $clustersBySequence{$sequence}[0] = $header;
                } else {
                    push( @{ $clustersBySequence{$sequence} }, $header );
                }

                $qualityByHeader{$header} = $quality;
            }
        } else {
            local $RS = '>';
            while ( my $fasta_record = <$IN> ) {
                chomp($fasta_record);
                my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
                my $header   = shift(@lines);
                my $sequence = lc( join( q{}, @lines ) );

                my $length = length $sequence;
                if ( $length == 0 ) {
                    next;
                }

                if ( !defined $clustersBySequence{$sequence} ) {
                    $clustersBySequence{$sequence}[0] = $header;
                } else {
                    push( @{ $clustersBySequence{$sequence} }, $header );
                }
            }
        }
        close $IN or croak("Cannot close file: $OS_ERROR\n");
    }

    my $i = 0;
    foreach my $seq ( keys %clustersBySequence ) {
        $clusterSize = scalar( @{ $clustersBySequence{$seq} } );
        if ( $clusterAll || $clusterSize > 1 || $clustersBySequence{$seq}[0] =~ /^C\d+?%/smx ) {
            print $TABLE "C${i}%", $clusterSize, '%', $label;
            if ( $clustersBySequence{$seq}[0] =~ /\{(.+?)\}/smx ) {
                my $annot = $1;
                print STDOUT '>C', $i, '%', $clusterSize, "%$label\{$annot\}\n", $seq, "\n";
            } else {
                print STDOUT '>C', $i, '%', $clusterSize, "%$label\n", $seq, "\n";
            }

            if ($fastQ) {
                foreach my $hdr ( @{ $clustersBySequence{$seq} } ) {
                    print $TABLE "\t", $hdr, "\t", $qualityByHeader{$hdr};
                }
            } else {
                foreach my $hdr ( @{ $clustersBySequence{$seq} } ) {
                    print $TABLE "\t", $hdr;
                }
            }
            print $TABLE "\n";
            $i++;
        } else {
            print STDOUT '>', $clustersBySequence{$seq}[0], "\n", $seq, "\n";
        }
    }
    close $TABLE or croak("Cannot close table: $OS_ERROR");
}
