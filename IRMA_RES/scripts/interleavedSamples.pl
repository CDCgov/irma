#!/usr/bin/env perl
#
# Description:      Interleaved sampling of a FASTA file.
#
# Date dedicated:   2022-07-06
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
#
# Citation:         Shepard SS, Davis CT, Bahl J, Rivailler P, York IA, Donis
#                   RO. LABEL: fast and accurate lineage assignment with
#                   assessment of H5N1 and H9N2 influenza A hemagglutinins. PLoS
#                   One. 2014;9(1):e86921. Published 2014 Jan 23.
#                   doi:10.1371/journal.pone.0086921
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
#  Please cite the manuscript  or author in any work or product based on this
#  material.

use warnings;
use strict;
use English qw(-no_match_vars);
use File::Basename;
use Getopt::Long qw(:config no_auto_abbrev);
use Carp qw(croak);

my ( $numberGroups, $fraction, $extension );
my ( $byReadPairs, $readZipped, $underscoreHeader, $fastQ, $removeEmptyFiles, $silenceSummary ) = ( 0, 0, 0, 0, 0, 0 );

GetOptions(
            'groups|G=i'           => \$numberGroups,
            'fraction|F=i'         => \$fraction,
            'fastQ|Q'              => \$fastQ,
            'by-read-pairs|P'      => \$byReadPairs,
            'read-zipped|Z'        => \$readZipped,
            'underscore-header|U'  => \$underscoreHeader,
            'extension|X:s'        => \$extension,
            'remove-empty-files|R' => \$removeEmptyFiles,
            'silence|S'            => \$silenceSummary
);

if ( scalar @ARGV < 2 ) {
    die(   "\n$PROGRAM_NAME <input.fasta> <out_prefix> [-G <#groups>|-F <denom-fraction>] [OPTIONS]\n"
         . "\t-F|--fraction POSITIVE_NUMBER\t\tFraction of dataset, using denominator D: 1/D.\n"
         . "\t-G|--groups POSITIVE_NUMBER\t\tNumber of datasets required.\n"
         . "\t-X|--extension\t\t\t\tExtension for output samplings.\n"
         . "\t-Q|--fastQ\t\t\t\tFastQ format for input and output.\n"
         . "\t-P|--by-read-pairs\t\t\tFastQ format for IN/OUT, interleave by read molecular ID (implies -Q).\n"
         . "\t-R|--remove-empty-files\t\t\tDeletes files at the end if no data was written.\n"
         . "\t-S|--silence\t\t\tSilence summary output of samples written.\n"
         . "\n" );
}

my $PROGRAM = basename( $PROGRAM_NAME, '.pl' );
if ($byReadPairs) {
    $fastQ = 1;
}

if ($fastQ) {
    $extension = 'fastq';
}

if ( defined $fraction && defined $numberGroups ) {
    die("$PROGRAM ERROR: specify Fraction OR the number of Groups.\n");
} elsif ( defined($numberGroups) ) {
    my $GROUP_LIMIT = 9999;
    if ( $numberGroups < 1 ) {
        die("ERROR: The number of groups must be more than zero.\n");
    } elsif ( $numberGroups > $GROUP_LIMIT ) {
        print STDERR "$PROGRAM WARNING: groups currently capped to $GROUP_LIMIT.\n";
        $numberGroups = $GROUP_LIMIT;
    }
    $fraction = 0;
} else {
    $numberGroups = $fraction;
    if ( $numberGroups < 2 ) {
        die("$PROGRAM ERROR: The denominator must be more than one.\n");
    }
    $fraction = 1;
}

my ( $fraction_handle, $fraction_filename );
my @handles = ();
my @count   = ();
my @files   = ();
if ($fraction) {
    my $suffix = 'th';
    if ( $numberGroups =~ /1$/smx ) {
        $suffix = 'st';
    } elsif ( $numberGroups =~ /2$/smx ) {
        $suffix = 'nd';
    } elsif ( $numberGroups =~ /3$/smx ) {
        $suffix = 'rd';
    }

    my $filename = sprintf( "%s_%d%s", $ARGV[1], $numberGroups, $suffix );
    if ( defined $extension ) {
        $filename .= ".$extension";
    } else {
        $filename .= '.fasta';
    }
    open( $fraction_handle, '>', $fraction_filename ) or die("$PROGRAM ERROR: cannot open $filename\n");
} else {
    for my $i ( 0 .. $numberGroups - 1 ) {
        my $filename = sprintf( "%s_%04d", $ARGV[1], ( $i + 1 ) );
        if ( defined $extension ) {
            $filename .= ".$extension";
        } else {
            $filename .= '.fasta';
        }
        open( $handles[$i], '>', $filename ) or die("$PROGRAM ERROR: Cannot open $filename\n");
        $files[$i] = $filename;
        $count[$i] = 0;
    }
}

# process parameters
chomp(@ARGV);
my $IN;
if ($readZipped) {
    open( $IN, "zcat $ARGV[0] |" ) or die("$PROGRAM ERROR: Could not open $ARGV[0].\n");
} else {
    open( $IN, '<', $ARGV[0] ) or die("$PROGRAM ERROR: Could not open $ARGV[0].\n");
}

sub get_mol_id {
    my $s       = shift // q{};
    my $default = shift // '0';
    my $pre     = substr( $s, 0, 4 );
    my $id;

    if ( index( $s, q{ } ) != -1 ) {
        my @remainder = ();
        ( $id, @remainder ) = split q{ }, $s;
        if ( ( $pre eq '@SRR' || $pre eq '@DRR' || $pre eq '@ERR' ) && index( $id, '.' ) != -1 ) {

            # SRA format
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $id = join( '.', @fields[0 .. 1] );
            }
        }

    } elsif ( index( $s, '/' ) != -1 ) {

        # Legacy Illumina
        ($id) = split '/', $s;
    } else {

        ## ASSUME SRA with phony Illumina
        if ( ( $pre eq '@SRR' || $pre eq '@DRR' || $pre eq '@ERR' ) && index( $s, '.' ) != -1 ) {
            my @remainder = ();
            ( $id, @remainder ) = split '_', $s;
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $id = join( '.', @fields[0 .. 1] );
            }
        } else {

            # IRMA Illumina output
            my @fields = split ':', $s;
            if ( scalar @fields > 6 ) {
                my $umi;
                ($umi) = split '_', $fields[6];
                $id = join( ':', ( @fields[0 .. 5], $umi // q{} ) );
            }
        }
    }

    return $id;
}

my $id = 0;
if ($fastQ) {
    local $RS = "\n";
    if ($byReadPairs) {
        my %indexByMolID = ();
        while ( my $hdr = <$IN> ) {
            my $seq     = <$IN>;
            my $junk    = <$IN>;
            my $quality = <$IN>;
            chomp($quality);

            my $index;
            my $molID = get_mol_id($hdr);
            if ( length $molID > 0 ) {
                if ( defined $indexByMolID{$molID} ) {
                    $index = $indexByMolID{$molID};
                } else {
                    $index = $id % $numberGroups;
                    $indexByMolID{$molID} = $index;
                    $id++;
                    $count[$index]++;
                }
            } else {
                die("Irregular header for fastQ read pairs: $hdr\n");
            }

            if ($underscoreHeader) { $hdr =~ tr/ /_/; }

            if ( !$fraction ) {
                my $handle = $handles[$index];
                print $handle $hdr, $seq, $junk, $quality, "\n";
            } elsif ( $index == 0 ) {
                print $fraction_handle $hdr, $seq, $junk, $quality, "\n";
            }
        }
    } else {
        while ( my $hdr = <$IN> ) {
            my $seq     = <$IN>;
            my $junk    = <$IN>;
            my $quality = <$IN>;
            chomp($quality);

            if ($underscoreHeader) { $hdr =~ tr/ /_/; }

            my $index = $id % $numberGroups;
            $id++;
            $count[$index]++;
            if ( !$fraction ) {
                my $handle = $handles[$index];
                print $handle $hdr, $seq, $junk, $quality, "\n";
            } elsif ( $index == 0 ) {
                print $fraction_handle $hdr, $seq, $junk, $quality, "\n";
            }
        }
    }
} else {
    local $RS = ">";
    while ( my $fasta_record = <$IN> ) {
        chomp($fasta_record);
        my @lines    = split( /\r\n|\n|\r/smx, $fasta_record );
        my $header   = shift(@lines);
        my $sequence = lc( join( q{}, @lines ) );

        if ( length($sequence) == 0 ) {
            next;
        }

        if ($underscoreHeader) { $header =~ tr/ /_/; }
        my $index = $id % $numberGroups;
        $id++;
        $count[$index]++;
        if ( !$fraction ) {
            my $handle = $handles[$index];
            print $handle '>', $header, "\n", $sequence, "\n";
        } elsif ( $index == 0 ) {
            print $fraction_handle '>', $header, "\n", $sequence, "\n";
        }
    }
}
close($IN) or croak("Cannot close $ARGV[0]");

if ($fraction) {
    close($fraction_handle) or croak("Cannot close: $fraction_filename\n");
    if ( !$silenceSummary ) {
        print STDOUT "\n Total\t  Got\tSample Name\n";
        print STDOUT '----------------------------------------------------', "\n";
        printf( "%6d\t%5d\t%s\n", $id, $count[0], $fraction_filename );
        print STDOUT '----------------------------------------------------', "\n";
    }
} else {
    foreach my $handle (@handles) {
        close($handle) or die("Cannot close $handle\n");
    }

    if ( !$silenceSummary ) {
        print STDOUT "\n Total\t  Got\tSample Name\n";
        print STDOUT '----------------------------------------------------', "\n";
        for my $i ( 0 .. $numberGroups - 1 ) {
            printf( "%6d\t%5d\t%s\n", $id, $count[$i], $files[$i] );

        }
        print STDOUT '----------------------------------------------------', "\n";
    }

    for my $i ( 0 .. $numberGroups - 1 ) {
        if ( defined $removeEmptyFiles && $count[$i] == 0 ) {
            unlink( $files[$i] );
        }

    }
}