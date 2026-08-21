#!/usr/bin/env perl
# Filename:         call
# Description:      IRMA variant calling and final consensus generation.
#
# Author:           Samuel S. Shepard, Centers for Disease Control and Prevention

## no critic (ControlStructures::ProhibitCascadingIfEls,Subroutines::RequireArgUnpacking)
use 5.016001;
use warnings;
use strict;

use Storable;
use English qw(-no_match_vars);
use Getopt::Long;
use Carp qw(croak);

#use Data::Dumper;

my ( $printAllAlleles, $sigLevel, $pairedStats, $autoFreq );

my $noGap      = 0;        # no gap allele
my $minCount   = 2;        # minimum allele count
my $minFreq    = 0.005;    # minimum allele frequency
my $minFreqIns = 0.005;    # minimum insertion frequency
my $minFreqDel = 0.005;    # minimum deletion frequency
my $minConf    = 0.5;      # minimum confidence not machine error
my $minQuality = 20;       # minimum average allele quality
my $minTotal   = 2;        # minimum total coverage depth

GetOptions(
            'no-gap-allele|G'            => \$noGap,
            'min-freq|F=f'               => \$minFreq,
            'min-insertion-freq|I=f'     => \$minFreqIns,
            'min-deletion-freq|D=f'      => \$minFreqDel,
            'min-count|C=i'              => \$minCount,
            'min-quality|Q=i'            => \$minQuality,
            'min-total-col-coverage|T=i' => \$minTotal,
            'print-all-sites|P'          => \$printAllAlleles,
            'conf-not-mac-err|M=f'       => \$minConf,
            'sig-level|S=f'              => \$sigLevel,
            'paired-error|E=s'           => \$pairedStats,
            'auto-min-freq|A'            => \$autoFreq
);

if ( $minCount < 0 )   { $minCount   = 0; }
if ( $minFreq < 0 )    { $minFreq    = 0; }
if ( $minFreqIns < 0 ) { $minFreqIns = 0; }
if ( $minFreqDel < 0 ) { $minFreqDel = 0; }
if ( $minConf < 0 )    { $minConf    = 0; }
if ( $minQuality < 0 ) { $minQuality = 0; }
if ( $minTotal < 0 )   { $minTotal   = 2; }

if ( scalar(@ARGV) < 3 ) {
    die(   "Usage:\n\tperl $PROGRAM_NAME <ref> <prefix> <aln.sto> <...>\n"
         . "\t\t-G|--no-gap-allele\t\t\tDo not count gaps alleles as variants.\n"
         . "\t\t-F|--min-freq <FLT>\t\t\tMinimum frequency for a variant to be processed. Default = 0.01.\n"
         . "\t\t-C|--min-count <INT>\t\t\tMinimum count of variant. Default = 2.\n"
         . "\t\t-Q|--min-quality <INT>\t\t\tMinimum average variant quality, preprocesses data. Default = 20.\n"
         . "\t\t-T|--min-total-col-coverage <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n"
         . "\t\t-P|--print-all-vars\t\t\tPrint all variants.\n"
         . "\t\t-M|--conf-not-mac-err <FLT>\t\tConfidence not machine error allowable minimum. Default = 0.5\n"
         . "\t\t-S|--sig-level <FLT>\t\t\tSignificance test (90, 95, 99, 99.9) variant is not machine error.\n"
         . "\t\t-E|--paired-error <FILE>\t\tFile with paired error estimates.\n"
         . "\t\t-A|--auto-min-freq\t\t\tAutomatically find minimum frequency heuristic.\n"
         . "\n" );
}

# FUNCTIONS #
sub calcProb($$) {
    my ( $w, $e ) = @_;
    if ( $e > $w ) {
        return 0;
    } else {
        return ( ( $w - $e ) / $w );
    }
}

sub lgg($) {
    return log( $_[0] ) / log(10);
}

sub max($$) {
    if ( $_[0] > $_[1] ) {
        return $_[0];
    } else {
        return $_[1];
    }
}

sub min($$) {
    if ( $_[0] < $_[1] ) {
        return $_[0];
    } else {
        return $_[1];
    }
}

sub avg($$) {
    return ( ( $_[0] + $_[1] ) / 2 );
}

sub toIndicesZero($) {
    my ($aln)  = @_;
    my $coords = q{};
    my $first  = q{};
    my $index  = 0;
    my $length = 0;

    while ( $aln =~ m/(\.+|[^.]+)/gsmx ) {
        ( $length, $first ) = ( length($1), substr( $1, 0, 1 ) );
        if ( $first ne '.' ) {
            $coords .= $index . ',';
            $index += $length;
            $coords .= ( $index - 1 ) . ';';
        } else {
            $index += $length;
        }
    }
    chop($coords);
    return $coords;
}
#############

my $takeSig = 0;
my ( $kappa, $kappa2, $eta, $gamma1, $gamma2 );
if ( defined($sigLevel) ) {
    $takeSig = 1;
    if ( $sigLevel >= 1 ) {
        $sigLevel /= 100;
    }

    if ( $sigLevel >= .999 ) {
        $kappa = 3.090232;
    } elsif ( $sigLevel >= .99 ) {
        $kappa = 2.326348;
    } elsif ( $sigLevel >= .95 ) {
        $kappa = 1.644854;
    } elsif ( $sigLevel >= .90 ) {
        $kappa = 1.281552;
    } else {
        $kappa = 3.090232;
    }

    ### second order correction ###
    $kappa2 = $kappa**2;
    $eta    = $kappa2 / 3 + 1 / 6;
    $gamma1 = $kappa2 * ( 13 / 18 ) + 17 / 18;
    $gamma2 = $kappa2 * ( 1 / 18 ) + 7 / 36;
    ###############################

    sub UB($$) {
        my $p = $_[0];
        my $N = $_[1];

        if ( $N <= 0 ) {
            print STDERR "Unexpected error: $N coverage depth.\n";
            return 0;
        }

        if ( $p == 1 ) {
            return 1;
        }

        # Let b = -1, so V = u - u^2
        my $V = $p - $p**2;

        # And N + 2*eta
        my $u2 = ( $p * $N + $eta ) / ( $N + 2 * $eta );

        my $inRoot = $V + ( $gamma2 - $gamma1 * $V ) / $N;

        # N < gamma1 - 4*gamma2
        # N < kappa^2/2 + 31/18
        if ( $inRoot < 0 ) {
            return 1;
        } else {
            my $UB = $u2 + $kappa * sqrt($inRoot) / sqrt($N);
            return max( min( $UB, 1 ), 0 );
        }
    }
}

# Consider implementing multiple references
local $RS = ">";
my ( $REF, $REF_NAME, $REF_SEQ, $REF_LEN );
open( $REF, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( my $fasta_record = <$REF> ) {
    chomp($fasta_record);
    my @lines = split( /\r\n|\n|\r/smx, $fasta_record );
    $REF_NAME = shift(@lines);
    $REF_SEQ  = join( q{}, @lines );
    if ( length($REF_SEQ) < 1 ) {
        next;
    }
    $REF_LEN = length($REF_SEQ);
    last;
}
close $REF or croak("Cannot close file: $OS_ERROR\n");
if ( !defined $REF_LEN ) { die("No reference found.\n"); }

my ( $DE, $PE, $IE, $is_paired ) = ( 0, 0, 0, 0 );
if ( defined $pairedStats ) {
    local $RS = "\n";
    my %pStats = ();
    my $PSF;
    open( $PSF, '<', $pairedStats ) or die("Cannot open $pairedStats for reading.\n");
    while ( my $line = <$PSF> ) {
        chomp($line);
        my ( $rn, $type, $value ) = split( "\t", $line );
        $pStats{$rn}{$type} = $value;
    }
    close $PSF or croak("Cannot close file: $OS_ERROR\n");

    $DE        = $pStats{$REF_NAME}{'MinimumDeletionErrorRate'};
    $PE        = $pStats{$REF_NAME}{'ExpectedErrorRate'};
    $IE        = $pStats{$REF_NAME}{'MinimumInsertionErrorRate'};
    $is_paired = 1;
}

my %icTable    = ();
my %iqTable    = ();
my @cTable     = ();
my @qTable     = ();
my %alignments = ();
my @data       = ();
my %varLine    = ();
my %variants   = ();
my %dcTable    = ();

foreach my $i ( 2 .. $#ARGV ) {
    @data = @{ retrieve( $ARGV[$i] ) };

    # combine alignments
    foreach my $aln ( keys( %{ $data[5] } ) ) {
        $alignments{$aln} += $data[5]{$aln};
    }

    # combine allele and quality counts
    foreach my $p ( 0 .. ( $REF_LEN - 1 ) ) {

        # Missing sites or alleles stay missing: no zero-valued entries are created.
        my $counts    = $data[0][$p] // next;
        my $qualities = $data[2][$p] // {};

        foreach my $allele ( keys %{$counts} ) {
            if ( defined $counts->{$allele} ) {
                $cTable[$p]{$allele} += $counts->{$allele};
            }

            if ( defined $qualities->{$allele} ) {
                $qTable[$p]{$allele} += $qualities->{$allele};
            }
        }
    }

    # handle insertion data
    foreach my $p ( keys( %{ $data[1] } ) ) {
        foreach my $insert ( keys( %{ $data[1]{$p} } ) ) {
            $icTable{$p}{$insert} += $data[1]{$p}{$insert};
            $iqTable{$p}{$insert} += $data[3]{$p}{$insert};
        }
    }

    # handle deletion data
    foreach my $p ( keys( %{ $data[4] } ) ) {
        foreach my $inc ( keys( %{ $data[4]{$p} } ) ) {
            $dcTable{$p}{$inc} += $data[4]{$p}{$inc};
        }
    }
}

#print STDERR Dumper(@cTable);

my $prefix = $ARGV[1];

my $ALLA;
if ($printAllAlleles) {
    open( $ALLA, '>', $prefix . '-allAlleles.txt' ) or die("ERROR: cannot open $prefix-allAlleles.txt for writing.\n");
    print $ALLA 'Reference_Name',  "\t",       'Position', "\t";
    print $ALLA 'Allele',          "\t",       'Count',    "\t", 'Total', "\t", 'Frequency', "\t";
    print $ALLA 'Average_Quality', "\t",       'ConfidenceNotMacErr';
    print $ALLA "\t",              'PairedUB', "\t", 'QualityUB', "\t", 'Allele_Type', "\n";
}

open( my $VARS, '>', $prefix . '-variants.txt' ) or die("ERROR: cannot open $prefix-variants.txt for writing.\n");
print $VARS 'Reference_Name', "\t",                        'Position', "\t", 'Total';
print $VARS "\t",             'Consensus_Allele',          "\t",       'Minority_Allele';
print $VARS "\t",             'Consensus_Count',           "\t",       'Minority_Count';
print $VARS "\t",             'Consensus_Frequency',       "\t",       'Minority_Frequency';
print $VARS "\t",             'Consensus_Average_Quality', "\t",       'Minority_Average_Quality';
print $VARS "\t",             'ConfidenceNotMacErr',       "\t",       'PairedUB', "\t", 'QualityUB', "\n";

open( my $COVG, '>', $prefix . '-coverage.txt' ) or die("Cannot open $prefix-coverage.txt for writing.\n");
open( my $CONS, '>', $prefix . '.fasta' )        or die("Cannot open $prefix.fasta for writing.\n");
print $COVG
  "Reference_Name\tPosition\tCoverage Depth\tConsensus\tDeletions\tAmbiguous\tConsensus_Count\tConsensus_Average_Quality\n";
print $CONS '>', $REF_NAME, "\n";

my @alpha      = split( q{}, 'AaCcGgTtUuRrYySsWwKkMmBbDdHhVvNn-.' );
my %base_order = ();
@base_order{@alpha} = 0 .. $#alpha;

my $hFreq        = 0;
my %totals       = ();
my $consensusSeq = q{};
my $cons_p       = 0;
foreach my $p ( 0 .. ( $REF_LEN - 1 ) ) {
    my @site_alleles =
      sort { ( $base_order{$a} // scalar(@alpha) ) <=> ( $base_order{$b} // scalar(@alpha) ) or $a cmp $b }
      keys( %{ $cTable[$p] } );
    my $nAlleles = scalar(@site_alleles);

    my $consensus       = '.';
    my $conCount        = 0;
    my $canonical_total = 0;

    if ( $nAlleles > 0 ) {

        # consensus for nAlleles = 1
        $consensus       = $site_alleles[0];
        $conCount        = $cTable[$p]{$consensus};
        $canonical_total = $conCount;
    }

    # find the consensus allele for nAlleles ≥ 2
    foreach my $b ( 1 .. $#site_alleles ) {
        my $base = $site_alleles[$b];
        if ( $base ne '-' ) {
            if (    $consensus eq '-'
                 || $cTable[$p]{$base} > $conCount
                 || ( $cTable[$p]{$base} == $conCount && $qTable[$p]{$base} > $qTable[$p]{$consensus} ) ) {
                $conCount  = $cTable[$p]{$base};
                $consensus = $base;
            }
        }
        $canonical_total += $cTable[$p]{$base};
    }

    # Account for ambiguous (do not count as part of coverage)
    if ( defined $cTable[$p]{'N'} ) {
        $canonical_total -= $cTable[$p]{'N'};
    }

    #if ( $consensus eq '.' && $consensusSeq eq '' ) {
    #    next;
    #} else {
    $cons_p++;

    #}

    $consensusSeq .= $consensus;
    print $CONS $consensus;

    my $conFreq = 0;
    if ( $consensus eq 'N' ) {
        $conFreq = 'NA';
    } elsif ( $canonical_total != 0 ) {
        $conFreq = $conCount / $canonical_total;
    }

    my $conQuality = 0;
    if ( $conCount > 0 ) {
        $conQuality = ( $qTable[$p]{$consensus} - $conCount * 33 ) / $conCount;
    }

    if ( defined $cTable[$p]{'-'} ) {
        print $COVG $REF_NAME, "\t", ($cons_p), "\t", ( $canonical_total - $cTable[$p]{'-'} ), "\t", $consensus, "\t",
          $cTable[$p]{'-'};
    } else {
        print $COVG $REF_NAME, "\t", ($cons_p), "\t", $canonical_total, "\t", $consensus, "\t", 0;
    }

    if ( !defined $cTable[$p]{'N'} ) {
        print $COVG "\t", 0;
    } else {
        print $COVG "\t", $cTable[$p]{'N'};
    }
    print $COVG "\t", $conCount, "\t", $conQuality, "\n";

    foreach my $base (@site_alleles) {

        # plurality allele can be ATGC + "N" + "-"
        if ( $base eq $consensus ) {
            if ($printAllAlleles) {

                my ( $confidence, $quality, $pairedUB, $qualityUB );
                if ( $base eq 'N' ) {
                    $confidence = 'NA';
                    $quality    = $conQuality;
                    $pairedUB   = 'NA';
                    $qualityUB  = 'NA';
                } elsif ( $base eq '-' ) {
                    $confidence = 'NA';
                    $quality    = 'NA';
                    $pairedUB   = UB( $DE, $canonical_total );
                    $qualityUB  = 0;
                } else {
                    my $ee = 1 / ( 10**( $conQuality / 10 ) );

                    $confidence = calcProb( $conFreq, $ee );
                    $quality    = $conQuality;
                    $pairedUB   = UB( $PE, $canonical_total );
                    $qualityUB  = UB( $ee, $canonical_total );
                }
                print $ALLA $REF_NAME, "\t", $cons_p, "\t", $base, "\t", $conCount, "\t", $canonical_total, "\t", $conFreq,
                  "\t", $quality, "\t", $confidence, "\t", $pairedUB, "\t", $qualityUB, "\t", 'Consensus', "\n";
            }
        } elsif ( $base ne 'N' ) {

            # no zero count columns
            my $count = $cTable[$p]{$base};
            if ( $count == 0 ) {
                next;
            }

            my $freq = $count / $canonical_total;
            my $quality;
            if ( $base ne '-' ) {
                $quality = ( $qTable[$p]{$base} - $count * 33 ) / $count;
            } else {
                $quality = $minQuality;
            }

            my ( $confidence, $pairedUB, $qualityUB );
            if ( $base eq '-' ) {
                $quality    = 'NA';
                $confidence = 'NA';
                $pairedUB   = UB( $DE, $canonical_total );
                $qualityUB  = 0;
            } elsif ( $consensus eq 'N' ) {
                $freq       = 'NA';
                $confidence = 'NA';
                $pairedUB   = 'NA';
                $qualityUB  = 'NA';
            } else {

                # quality-based estimated error
                my $ee = 1 / ( 10**( $quality / 10 ) );

                $confidence = calcProb( $freq, $ee );
                $pairedUB   = UB( $PE, $canonical_total );
                $qualityUB  = UB( $ee, $canonical_total );

                # Deletions have no expected error from quality scores Even if
                # were deletion minor variants were allowed, they would be
                # unikely to contribute to the  auto-frequency heuristic
                # threshold.
                if ( $freq <= $ee && $freq > $hFreq ) {
                    $hFreq = $freq;
                }
            }

            # Valid called variant: ATGC + "-"
            # IRMA v1.1.0 does not allow minor variants with ambiguous consensus
            if (    $consensus ne 'N'
                 && !( $noGap && $base eq '-' )
                 && $freq >= $minFreq
                 && $count >= $minCount
                 && $quality >= $minQuality
                 && $canonical_total >= $minTotal ) {

                if ($printAllAlleles) {
                    print $ALLA $REF_NAME, "\t", $cons_p, "\t", $base, "\t", $count, "\t", $canonical_total, "\t", $freq,
                      "\t", $quality, "\t", $confidence, "\t", $pairedUB, "\t", $qualityUB, "\t", 'Minority', "\n";
                }

                if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) {
                    next;
                }

                $variants{$p}{$base} = $freq;
                $varLine{$p}{$base}  = $REF_NAME . "\t" . $cons_p . "\t" . $canonical_total . "\t";
                $varLine{$p}{$base} .= $consensus . "\t" . $base . "\t" . $conCount . "\t" . $count . "\t";
                $varLine{$p}{$base} .= $conFreq . "\t" . $freq . "\t" . $conQuality . "\t" . $quality . "\t";
                $varLine{$p}{$base} .= $confidence . "\t" . $pairedUB . "\t" . $qualityUB . "\n";

            } elsif ($printAllAlleles) {

                # any minor variant: ATGC + "-"
                print $ALLA $REF_NAME, "\t", $cons_p, "\t", $base, "\t", $count, "\t", $canonical_total, "\t", $freq, "\t",
                  $quality, "\t", $confidence, "\t", $pairedUB, "\t", $qualityUB, "\t", 'Minority', "\n";
            }
        }
    }
}
print $CONS "\n";
close $CONS or croak("Cannot close file: $OS_ERROR\n");
close $COVG or croak("Cannot close file: $OS_ERROR\n");
close $ALLA or croak("Cannot close file: $OS_ERROR\n");

# Revise variants according to the heuristic auto-frequency
# and print the variant table.
foreach my $p ( sort { $a <=> $b } keys(%varLine) ) {
    foreach my $base ( sort { $varLine{$p}{$a} cmp $varLine{$p}{$b} } keys( %{ $varLine{$p} } ) ) {
        if ($autoFreq) {
            if ( $variants{$p}{$base} > $hFreq ) {
                print $VARS $varLine{$p}{$base};
            } else {
                delete( $variants{$p}{$base} );
            }
        } else {
            print $VARS $varLine{$p}{$base};
        }
    }
}
close $VARS or croak("Cannot close file: $OS_ERROR\n");

my %coordSupport = ();
my %coordList    = ();
foreach my $aln ( keys(%alignments) ) {
    $coordList{ toIndicesZero($aln) } += $alignments{$aln};
}

foreach my $listOfCoords ( keys(%coordList) ) {
    my @coords = split( ';', $listOfCoords );
    foreach my $coord (@coords) {
        my ( $start, $stop ) = split( ',', $coord );
        $coordSupport{$start}{$stop} += $coordList{$listOfCoords};
    }
}

my %coordStops  = ();
my @coordStarts = sort { $a <=> $b } keys(%coordSupport);
foreach my $start (@coordStarts) {
    $coordStops{$start} = [sort { $b <=> $a } keys( %{ $coordSupport{$start} } )];
}

open( my $INSV, '>', $prefix . '-insertions.txt' ) or die("ERROR: cannot open $prefix-insertions.txt for writing.\n");
print $INSV "Reference_Name\tUpstream_Position\tInsert\tContext\tCalled\tCount\tTotal\tFrequency\tAverage_Quality",
  "\tConfidenceNotMacErr\tPairedUB\tQualityUB", "\n";
foreach my $p ( sort { $a <=> $b } keys(%icTable) ) {
    my $pp    = $p + 1;
    my $total = 0;
    foreach my $start (@coordStarts) {
        if ( $start <= $p ) {
            foreach my $stop ( @{ $coordStops{$start} } ) {
                if ( $pp <= $stop ) {
                    $total += $coordSupport{$start}{$stop};
                } else {
                    last;
                }
            }
        } else {
            last;
        }
    }

    foreach my $insert ( sort { $a cmp $b } keys( %{ $icTable{$p} } ) ) {
        my $count = $icTable{$p}{$insert};
        if ( $count < $minCount ) { next; }

        my $called  = "TRUE";
        my $quality = 0;
        my $freq    = 0;
        if ( $count > 0 )             { $quality = $iqTable{$p}{$insert} / $count; }
        if ( $quality < $minQuality ) { $called  = "FALSE"; }
        if ( $total > 0 )             { $freq    = $count / $total; }

        if ( $freq < $minFreqIns || $total < $minTotal ) { $called = "FALSE"; }

        my $EE         = 1 / ( 10**( $quality / 10 ) );
        my $confidence = calcProb( $freq, $EE );
        my $pairedUB   = UB( $IE, $total );
        my $qualityUB  = UB( $EE, $total );

        if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) { $called = "FALSE"; }
        my ( $left_flanking, $right_flanking ) = ( q{}, q{} );
        if ( $p < 5 ) {
            $left_flanking = substr( $consensusSeq, 0, $p + 1 );
        } else {
            $left_flanking = substr( $consensusSeq, $p - 4, 5 );
        }

        if ( $p > ( $REF_LEN - 6 ) ) {
            $right_flanking = substr( $consensusSeq, $pp, $REF_LEN - $pp );
        } else {
            $right_flanking = substr( $consensusSeq, $pp, 5 );
        }

        print $INSV $REF_NAME, "\t", ( $p + 1 ), "\t", uc($insert), "\t", lc($left_flanking), uc($insert),
          lc($right_flanking), "\t", $called, "\t", $count, "\t", $total, "\t", $freq, "\t", $quality, "\t",
          $confidence, "\t", $pairedUB, "\t", $qualityUB, "\t", "\n";
    }
}
close $INSV or croak("Cannot close file: $OS_ERROR\n");

open( my $DELV, '>', $prefix . '-deletions.txt' ) or die("ERROR: cannot open $prefix-deletions.txt for writing.\n");
print $DELV "Reference_Name\tUpstream_Position\tLength\tContext\tCalled\tCount\tTotal\tFrequency\tPairedUB\n";
foreach my $p ( sort { $a <=> $b } keys(%dcTable) ) {
    foreach my $inc ( sort { $a <=> $b } keys( %{ $dcTable{$p} } ) ) {
        my $count = $dcTable{$p}{$inc};
        if ( $count < $minCount ) { next; }
        my $total  = 0;
        my $pp     = $p + $inc + 1;
        my $called = "TRUE";

        # get depth
        foreach my $start (@coordStarts) {
            if ( $start <= $p ) {
                foreach my $stop ( @{ $coordStops{$start} } ) {
                    if ( $pp <= $stop ) {
                        $total += $coordSupport{$start}{$stop};
                    } else {
                        last;
                    }
                }
            } else {
                last;
            }
        }

        my $freq = 0;
        if ( $total > 0 ) {
            $freq = $count / $total;
        }

        if ( $freq < $minFreqDel || $total < $minTotal ) { $called = "FALSE"; }

        my $pairedUB = UB( $DE, $total );
        if ( $freq <= $pairedUB ) { $called = "FALSE"; }

        my ( $left_flanking, $right_flanking ) = ( q{}, q{} );
        if ( $p < 5 ) {
            $left_flanking = substr( $consensusSeq, 0, $p + 1 );
        } else {
            $left_flanking = substr( $consensusSeq, $p - 4, 5 );
        }

        if ( $p > ( $REF_LEN - 6 - $inc ) ) {
            $right_flanking = substr( $consensusSeq, $pp, $REF_LEN - $pp );
        } else {
            $right_flanking = substr( $consensusSeq, $pp, 5 );
        }

        my $mid = '-' x $inc;
        print $DELV $REF_NAME, "\t", ( $p + 1 ), "\t", $inc, "\t", $left_flanking, $mid, $right_flanking;
        print $DELV "\t", $called, "\t", $count, "\t", $total, "\t", $freq, "\t", $pairedUB, "\n";
    }
}
close $DELV or croak("Cannot close file: $OS_ERROR\n");

my $variantCount = 0;
foreach my $variantPosition ( keys(%variants) ) {
    $variantCount += scalar( keys( %{ $variants{$variantPosition} } ) );
}

if ( $variantCount > 1 ) {
    my $varFile = $prefix . '-vars.sto';
    my $patFile = $prefix . '-pats.sto';
    store( \%variants, $varFile );

    my %readPats = ();
    my @vars     = sort { $a <=> $b } keys(%variants);
    foreach my $sequence ( keys(%alignments) ) {
        my $aln = q{};
        foreach my $pos (@vars) {
            $aln .= substr( $sequence, $pos, 1 );
        }

        if ( $aln !~ /^[.N]+$/smx ) {
            $readPats{$aln} += $alignments{$sequence};
        }
    }

    store( \%readPats, $patFile );
}
