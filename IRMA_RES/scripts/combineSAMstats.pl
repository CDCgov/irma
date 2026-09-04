#!/usr/bin/env perl
# combineSAMstats.pl
#
# Sam Shepard - 2020
#
# Description: combine assembly statistics for iterative final assembly generation

use POSIX;
use Getopt::Long;
use English qw(-no_match_vars);

Getopt::Long::Configure('no_ignore_case');

GetOptions(
            'name|N=s'                      => \$name,
            'insertion-threshold|I=f'       => \$insertionThreshold,
            'deletion-threshold|D=f'        => \$deletionThreshold,
            'insertion-depth-threshold|i=i' => \$insertionDepthThreshold,
            'deletion-depth-threshold|d=i'  => \$deletionDepthThreshold,
            'alternative-threshold|A=f'     => \$alternativeThreshold,
            'alternative-count|C=i'         => \$alternativeCount,
            'store-stats|S=s'               => \$storeStats,
            'mark-deletions|M'              => \$markDeletions,
);

if ( scalar(@ARGV) < 2 ) {
    $message = "Usage:\t$PROGRAM_NAME [options] <REF> <STAT1> <...>\n";
    $message .= "\t\t-N|--name <STR>\t\t\t\tName of consensus sequence.\n";
    $message .=
      "\t\t-I|--insertion-threshold <#>\t\tInsertion frequency where consensus is altered. Default = 0.15 or 15%.\n";
    $message .= "\t\t-D|--deletion-threshold <#>\t\tDeletion frequency where consensus is altered. Default = 0.75 or 75%.\n";
    $message .=
"\t\t-i|--insertion-depth-threshold <#>\tInsertion coverage depth where consensus is edit (given frequency). Default = 1.\n";
    $message .=
"\t\t-d|--deletion-depth-threshold <#>\tDeletion coverage depth where consensus is edited (given frequency). Default = 1.\n";
    $message .=
      "\t\t-A|--alternative-threshold <#>\t\tFrequency where alternative reference allele is changed. Default is off.\n";
    $message .=
      "\t\t-C|--alternative-count <#>\t\tVariant count where alternative reference allele is changed. Default is off.\n";
    $message .= "\t\t-S|--store-stats <FILE>\t\t\tSave aggregate stats to a .sto file.\n";
    $message .=
      "\t\t-M|--mark-deletions\t\t\tOutput '-' for deletions in the consensus instead ommitting the deleted states.\n";
    die( $message . "\n" );

    # TODO: revisit minimum dropout edge feature in Rust port
}

open( REF, '<', $ARGV[0] ) or die("$PROGRAM_NAME ERROR: cannot open REF $ARGV[0] for reading.\n");
local $RS = ">";
while ( $record = <REF> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/smx, $record );
    $REF_NAME = shift(@lines);
    $REF_SEQ  = join( q{}, @lines );
    if ( length($REF_SEQ) < 1 ) {
        next;
    }
    $N = length($REF_SEQ);
    last;
}
close(REF);
if ( !defined($N) ) { die("$PROGRAM_NAME ERROR: no reference found in $ARGV[0].\n"); }

# Insert if >= thresholds
# Given as >= T AND >= C
if ( !defined($insertionThreshold) )      { $insertionThreshold      = 0.15; }
if ( !defined($insertionDepthThreshold) ) { $insertionDepthThreshold = 1; }

# Delete if >= thresholds
# Given as < T OR < C implies >= T AND >= C
if ( !defined($deletionThreshold) )      { $deletionThreshold      = 0.75; }
if ( !defined($deletionDepthThreshold) ) { $deletionDepthThreshold = 1; }

# Alternative if >= T AND C
# Given as < T OR < C implies >= T AND >= C
if ( !defined($alternativeThreshold) ) { $alternativeThreshold = 2; }
if ( !defined($alternativeCount) )     { $alternativeCount     = LONG_MAX; }

# Flags
$markDeletions = defined($markDeletions) ? 1 : 0;
$storeStats    = defined($storeStats)    ? 1 : 0;

# Aggregate data
@agg_counts       = ();
@agg_quals        = ();
%agg_inserts      = ();
%agg_insert_quals = ();

local $RS = "\n";
for my $i ( 1 .. $#ARGV ) {
    open( my $STAT, '<', $ARGV[$i] ) or die("$PROGRAM_NAME ERROR: cannot open STAT $ARGV[$i] for reading.\n");
    while ( my $line = <$STAT> ) {
        chomp($line);
        my ( $type, $pos, $allele, $partition_count, $partition_enc_qualities ) = split( "\t", $line );
        if ( $type eq 'M' ) {
            $agg_counts[$pos]{$allele} += $partition_count;
            $agg_quals[$pos]{$allele}  += $partition_enc_qualities;
        } elsif ( $type eq 'I' ) {
            $agg_inserts{$pos}{$allele}      += $partition_count;
            $agg_insert_quals{$pos}{$allele} += $partition_enc_qualities;
        } else {
            die("$PROGRAM_NAME ERROR: unrecognized record type '$type' in $ARGV[$i].\n");
        }
    }
    close($STAT);
}

my @alpha      = split( q{}, 'AaCcGgTtUuRrYySsWwKkMmBbDdHhVvNn-.' );
my %base_order = ();
@base_order{@alpha} = 0 .. $#alpha;

my @totals = ();
for my $p ( 0 .. ( $N - 1 ) ) {

    # Different from call.pl
    # includes all counts for sake of deletion editing
    my $total = 0;
    foreach my $count ( values( %{ $agg_counts[$p] } ) ) {
        $total += $count;
    }
    $totals[$p] = $total;
}

my ( $header,    $header2 )     = ( q{}, q{} );
my ( $consensus, $alternative ) = ( q{}, q{} );
if ($name) {
    $header  = '>' . $name . "\n";
    $header2 = $header;
} else {
    $header  = ">consensus\n";
    $header2 = ">alternative\n";
}

my ( @sorted_alleles, $cons );
for my $p ( 0 .. ( $N - 1 ) ) {
    @sorted_alleles = sort {
        $agg_counts[$p]{$b}                        <=> $agg_counts[$p]{$a}                      # desc count
          or $agg_quals[$p]{$b}                    <=> $agg_quals[$p]{$a}                       # desc quality
          or ( $base_order{$a} // scalar(@alpha) ) <=> ( $base_order{$b} // scalar(@alpha) )    # asc alpha order
          or $a cmp $b                                                                          # ASCII byte order
      }
      keys( %{ $agg_counts[$p] } );
    $cons = $sorted_alleles[0] // q{};

    if ( $cons ne '-' ) {
        $consensus .= $cons;

        # alternative non-gap
        if ( scalar @sorted_alleles > 1 && $sorted_alleles[1] ne '-' ) {
            $altCount = $agg_counts[$p]{ $sorted_alleles[1] };
            $altFreq  = $altCount / $totals[$p];

            if ( $altFreq < $alternativeThreshold || $altCount < $alternativeCount ) {
                $alternative .= $cons;
            } else {
                $alternative .= $sorted_alleles[1];
            }
        } else {
            $alternative .= $cons;
        }
    } else {

        # Plurality consensus is '-'
        $freq = $agg_counts[$p]{$cons} / $totals[$p];

        # Ignore deletion if below threshold and there exists another allele
        if ( ( $agg_counts[$p]{$cons} < $deletionDepthThreshold || $freq < $deletionThreshold )
             && scalar(@sorted_alleles) > 1 ) {
            $consensus .= $sorted_alleles[1];

            # alternative non-gap
            if ( scalar @sorted_alleles > 2 ) {
                $altCount = $agg_counts[$p]{ $sorted_alleles[2] };
                $altFreq  = $altCount / $totals[$p];

                if ( $altFreq < $alternativeThreshold || $altCount < $alternativeCount ) {
                    $alternative .= $sorted_alleles[1];
                } else {
                    $alternative .= $sorted_alleles[2];
                }
            } else {
                $alternative .= $sorted_alleles[1];
            }

            # >= Thresholds for deletion OR just the deletion allele is found
            # Skip unless we are to mark the deletion
        } elsif ($markDeletions) {

            # Consensus is a deletion.
            # Note: the alternative consensus cannot differ in length, so the alternative allele is not evaluated.
            $consensus   .= $cons;
            $alternative .= $cons;
        }
    }

    if ( defined( $agg_inserts{$p} ) ) {
        @sortedIns =
          sort {
            $agg_inserts{$p}{$b} <=> $agg_inserts{$p}{$a}                                                   # desc count
              or ( $agg_insert_quals{$p}{$b} / length $b ) <=> ( $agg_insert_quals{$p}{$a} / length $a )    # desc quality
              or $a cmp $b                                                                                  # string order
          }
          keys( %{ $agg_inserts{$p} } );
        if ( $p < ( $N - 1 ) ) {
            $avgTotal = int( ( $totals[$p] + $totals[$p + 1] ) / 2 );
        } else {
            $avgTotal = $totals[$p];
        }

        $freq = $agg_inserts{$p}{ $sortedIns[0] } / $avgTotal;
        if ( $freq >= $insertionThreshold && $agg_inserts{$p}{ $sortedIns[0] } >= $insertionDepthThreshold ) {
            $consensus   .= lc( $sortedIns[0] );
            $alternative .= lc( $sortedIns[0] );
        }
    }
}

# print the consensus sequences
print $header, $consensus, "\n";
if ( $alternative ne q{} && $alternative ne $consensus ) {
    print $header2, $alternative, "\n";
}
