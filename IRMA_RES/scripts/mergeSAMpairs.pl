#!/usr/bin/env perl
# Filename:         mergeSAMpairs
# Description:      Reference-based read-pair merging for paired-end alignments.
#                   Detects error and corrects paired overlaps parsimoniously.
#
# Date dedicated:   2023-06-21
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

use 5.16.1;
use strict;
use warnings;
use Getopt::Long;
use English qw(-no_match_vars);
use Carp qw(croak);

use constant { false => 0, true => 1 };

my ( $storeStats, $bowtieFormat ) = ( false, false );
GetOptions( 'store-stats|S'   => \$storeStats,
            'bowtie-format|B' => \$bowtieFormat );

if ( scalar(@ARGV) != 3 ) {
    my $message = "Usage:\n\tperl $PROGRAM_NAME <ref> <sam> <prefix>\n";
    $message .= "\t\t-S|--use-storable\tStore statistics object..\n";
    die( $message . "\n" );
}

# FUNCTIONS #
sub condenseCigar($) {
    my $cig   = $_[0];
    my $cigar = q{};
    my $state = q{};
    while ( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/g ) {
        $state = $1;
        $cigar .= length($state);
        $cigar .= substr( $state, 0, 1 );
    }
    return $cigar;
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
#############

my $REisBase = qr/[ATCG]/smx;

local $RS = ">";
my $REF_LEN  = 0;
my $REF_NAME = q{};
my @REF_SEQ  = ();
my $REF;
open( $REF, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
while ( my $fasta_record = <$REF> ) {
    chomp($fasta_record);
    my @lines = split( /\r\n|\n|\r/, $fasta_record );
    $REF_NAME = shift(@lines);
    my $seq = join( q{}, @lines );
    if ( length($seq) < 1 ) {
        next;
    }
    @REF_SEQ = split( q{}, uc($seq) );
    $REF_LEN = length($seq);
    last;
}
close($REF) or croak("Cannot close file: $OS_ERROR\n");

local $RS = "\n";
my $SAM;
open( $SAM, '<', $ARGV[1] ) or die("Cannot open $ARGV[1] for reading.\n");
my @sam = <$SAM>;
chomp(@sam);
close($SAM) or croak("Cannot close file: $OS_ERROR\n");

sub sra_prefix {
    my $s = shift;
    $s = substr( $s, 0, 4 );
    return $s eq '@SRR' || $s eq '@DRR' || $s eq '@ERR';
}

sub get_molID_side {
    my $s       = shift // q{};
    my $default = shift // '0';
    my $pre     = substr( $s, 0, 3 );

    my ( $id, $read );

    if ( index( $s, q{ } ) != -1 ) {
        my @remainder = ();

        ( $id, @remainder ) = split q{ }, $s;
        if ( !( $pre eq 'SRR' || $pre eq 'DRR' || $pre eq 'ERR' ) || index( $id, '.' ) == -1 ) {

            # Illumina format
            ($read) = split ':', ( $remainder[0] // q{} );
        } else {

            # SRA format
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $read = $fields[2];
                $id   = join( '.', @fields[0 .. 1] );
            }
        }

    } elsif ( index( $s, '/' ) != -1 ) {

        # Legacy Illumina
        ( $id, $read ) = split '/', $s;
    } else {

        ## ASSUME SRA with phony Illumina
        if ( ( $pre eq 'SRR' || $pre eq 'DRR' || $pre eq 'ERR' ) && index( $s, '.' ) != -1 ) {
            my @remainder = ();
            ( $id, @remainder ) = split '_', $s;
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $read = $fields[2];
                $id   = join( '.', @fields[0 .. 1] );
            }
        } else {

            # IRMA Illumina output
            my @fields = split ':', $s;
            if ( scalar @fields > 6 ) {
                my $umi;
                ( $umi, $read ) = split '_', $fields[6];
                $id = join( ':', ( @fields[0 .. 5], $umi // q{} ) );
            }
        }
    }

    if ( !defined $read || $read lt '0' || $read gt '3' ) {
        $read = $default;
    }

    return ( $id // $s, $read // $default );
}

sub make_merged_qname {
    my $s = shift // q{};
    my $id;
    my $pre = substr( $s, 0, 3 );

    if ( index( $s, q{ } ) != -1 ) {
        my @remainder = ();
        ( $id, @remainder ) = split q{ }, $s;
        if ( !( $pre eq 'SRR' || $pre eq 'DRR' || $pre eq 'ERR' ) || index( $id, '.' ) == -1 ) {

            # Illumina format
            my @fields = split ':', ( $remainder[0] // q{} );
            $fields[0]    = '3';
            $remainder[0] = join( ':', @fields );

        } else {

            # SRA format
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $fields[2] = '3';
                $id = join( '.', @fields );
            } else {
                $id .= '.3';
            }
        }
        $s = join( q{ }, ( $id, @remainder ) );
    } elsif ( index( $s, '/' ) != -1 ) {

        # Legacy Illumina
        ($id) = split '/', $s;
        $s = $id . '/3';
    } else {

        ## ASSUME SRA with phony Illumina
        if ( ( $pre eq 'SRR' || $pre eq 'DRR' || $pre eq 'ERR' ) && index( $s, '.' ) != -1 ) {
            my @remainder = ();
            ( $id, @remainder ) = split '_', $s;
            my @fields = split '\.', $id;
            if ( scalar @fields == 3 ) {
                $fields[2] = '3';
                $id = join( '.', @fields );
            } else {
                $id .= '.3';
            }
            $s = join( '_', ( $id, @remainder ) );
        } else {

            # IRMA Illumina output
            my @fields = split ':', $s;
            if ( scalar @fields > 6 ) {
                my $umi;
                ($umi) = split '_', $fields[6];
                $fields[6] = ( $umi // q{} ) . '_3';
                $s = join( ':', @fields );
            }
        }
    }

    return $s;
}

my %pairs      = ();
my %insByIndex = ();
foreach my $K ( 0 .. $#sam ) {
    if ( substr( $sam[$K], 0, 1 ) eq '@' ) {
        next;
    }

    my ( $qname, $flag, $rn, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual ) = split( "\t", $sam[$K] );

    my ( $qMolID, $qSide ) = get_molID_side( $qname, '0' );

    # TO CONSIDER: switching to an order-based approach in the future for
    # everything
    if ($bowtieFormat) {
        if ( defined $pairs{$qMolID} ) {
            $qSide = 2;
        } else {
            $qSide = 1;
        }
    }

    if ( $REF_NAME eq $rn ) {
        my @NTs    = split( q{}, uc($seq) );
        my @QCs    = split( q{}, $qual );
        my @Qint   = unpack( "c* i*", $qual );
        my @cigars = split( q{}, $cigar );
        my $rpos   = $pos - 1;
        my $qpos   = 0;
        my ( $aln, $qAln ) = ( q{}, q{} );

        while ( $cigar =~ /(\d+)([MIDNSHP])/g ) {
            my $inc = $1;
            my $op  = $2;
            if ( $op eq 'M' ) {
                for ( 1 .. $inc ) {
                    $qAln .= $QCs[$qpos];
                    $aln  .= $NTs[$qpos];
                    $qpos++;
                    $rpos++;
                }
            } elsif ( $op eq 'D' ) {
                $qAln .= q{ } x $inc;
                $aln  .= '-' x $inc;
                for ( 1 .. $inc ) {
                    $rpos++;
                }
            } elsif ( $op eq 'I' ) {
                $insByIndex{$K}{ $rpos - 1 } = [substr( $seq, $qpos, $inc ), substr( $qual, $qpos, $inc )];
                $qpos += $inc;
            } elsif ( $op eq 'S' ) {
                $qpos += $inc;
            } elsif ( $op eq 'N' ) {
                $aln  .= 'N' x $inc;
                $qAln .= q{ } x $inc;
                $rpos += $inc;
            } elsif ( $op eq 'H' ) {
                next;
            } else {
                die("Extended CIGAR ($op) not yet supported.\n");
            }
        }

        $pairs{$qMolID}{$qSide} = [$aln, $qAln, $K, ( $pos - 1 ), ( $rpos - 1 ), $qname, $mapq];
    }
}

my %observations = ();
my $OSAM;
open( $OSAM, '>', $ARGV[2] . '.sam' ) or die("Cannot open $ARGV[2].sam\n");
foreach my $line (@sam) {
    if ( substr( $line, 0, 1 ) eq '@' ) {
        print $OSAM $line, "\n";
    } else {
        last;
    }
}
my ( $dmv, $obs, $fmv, $tmv ) = ( 0, 0, 0, 0 );
my ( $insObs, $insErr ) = ( 0, 0 );
foreach my $mID ( keys(%pairs) ) {
    my @mPairs = keys( %{ $pairs{$mID} } );
    if ( scalar(@mPairs) == 2 ) {
        my @a1 = @{ $pairs{$mID}{'1'} };
        my @a2 = @{ $pairs{$mID}{'2'} };
        my ( $s1, $e1 ) = ( $a1[3], $a1[4] );
        my ( $s2, $e2 ) = ( $a2[3], $a2[4] );

        my $start = min( $s1, $s2 );
        my $end   = max( $e1, $e2 );

        my $mSeq   = q{};
        my $cigars = q{};
        my $qSeq   = q{};

        my $K1     = $a1[2];
        my $K2     = $a2[2];
        my @bases1 = split( q{}, $a1[0] );
        my @bases2 = split( q{}, $a2[0] );
        my @quals1 = unpack( "c* i*", $a1[1] );
        my @quals2 = unpack( "c* i*", $a2[1] );

        my $FB1 = sub {
            my $i = $_[0];
            if ( $i < $s1 || $i > $e1 ) {
                return '.';
            } else {
                if ( ( $i - $s1 ) > $#bases1 ) { die( "Bad: " . ( $i - $s1 ) . " > " . $#bases1 . "\n" ); }
                return $bases1[$i - $s1];
            }
        };

        my $FQ1 = sub {
            my $i = $_[0];
            if ( $i < $s1 || $i > $e1 ) {
                return q{ };
            } else {
                return $quals1[$i - $s1];
            }
        };

        my $FB2 = sub {
            my $i = $_[0];
            if ( $i < $s2 || $i > $e2 ) {
                return '.';
            } else {
                return $bases2[$i - $s2];
            }
        };

        my $FQ2 = sub {
            my $i = $_[0];
            if ( $i < $s2 || $i > $e2 ) {
                return q{ };
            } else {
                return $quals2[$i - $s2];
            }
        };

        foreach my $i ( $start .. $end ) {
            my $x  = $FB1->($i);
            my $y  = $FB2->($i);
            my $qx = $FQ1->($i);
            my $qy = $FQ2->($i);
            my $r  = $REF_SEQ[$i];

            if ( $x ne '.' && $y ne '.' ) {
                $obs++;
                if ( $x eq $y ) {
                    if ( $x eq '-' ) {
                        $tmv++;
                        $cigars .= 'D';
                    } else {
                        if ( $x ne $r ) { $tmv++; }
                        $cigars .= 'M';
                        $mSeq   .= $x;
                        $qSeq   .= chr( max( $qx, $qy ) );
                    }
                } elsif ( $x eq $r ) {
                    $fmv++;
                    if ( $y eq '-' ) { $dmv++; }
                    $mSeq   .= $x;
                    $qSeq   .= chr($qx);
                    $cigars .= 'M';
                } elsif ( $y eq $r ) {
                    $fmv++;
                    if ( $x eq '-' ) { $dmv++; }
                    $mSeq   .= $y;
                    $qSeq   .= chr($qy);
                    $cigars .= 'M';
                } else {
                    $fmv++;
                    if ( $x =~ $REisBase && $y !~ $REisBase ) {
                        $cigars .= 'M';
                        $mSeq   .= $x;
                        $qSeq   .= chr($qx);
                        if ( $y eq '-' ) { $dmv++; }
                    } elsif ( $x !~ $REisBase && $y =~ $REisBase ) {
                        $cigars .= 'M';
                        $mSeq   .= $y;
                        $qSeq   .= chr($qy);
                        if ( $x eq '-' ) { $dmv++; }
                    } elsif ( $qx > ( $qy + 4 ) ) {
                        $cigars .= 'M';
                        $mSeq   .= $x;
                        $qSeq   .= chr($qx);
                        if ( $y eq '-' ) { $dmv++; }
                    } elsif ( $qy > ( $qx + 4 ) ) {
                        $cigars .= 'M';
                        $mSeq   .= $y;
                        $qSeq   .= chr($qy);
                        if ( $x eq '-' ) { $dmv++; }
                    } else {
                        $cigars .= 'M';
                        $mSeq   .= 'N';
                        $qSeq   .= chr( int( avg( $qx, $qy ) ) );
                    }
                }
            } elsif ( $x eq '.' && $y ne '.' ) {
                if ( $y eq '-' ) {
                    $cigars .= 'D';
                } else {
                    $cigars .= 'M';
                    $mSeq   .= $y;
                    $qSeq   .= chr($qy);
                }
            } elsif ( $x ne '.' && $y eq '.' ) {
                if ( $x eq '-' ) {
                    $cigars .= 'D';
                } else {
                    $cigars .= 'M';
                    $mSeq   .= $x;
                    $qSeq   .= chr($qx);
                }
            } else {
                $cigars .= 'N';
            }

            if ( defined $insByIndex{$K1}{$i} && defined $insByIndex{$K2}{$i} ) {
                my $ins1 = lc( $insByIndex{$K1}{$i}[0] );
                my $ins2 = lc( $insByIndex{$K2}{$i}[0] );
                $insObs++;
                if ( $ins1 eq $ins2 ) {
                    $mSeq .= $ins1;
                    my @qIns1   = split( q{}, $insByIndex{$K1}{$i}[1] );
                    my @qIns2   = split( q{}, $insByIndex{$K2}{$i}[1] );
                    my $qSeqNew = q{};
                    foreach my $qIndex ( 0 .. ( length($ins1) - 1 ) ) {
                        $qSeqNew .= chr( max( ord( $qIns1[$qIndex] ), ord( $qIns2[$qIndex] ) ) );
                    }
                    $qSeq   .= $qSeqNew;
                    $cigars .= 'I' x length($ins1);
                } elsif ( $ins2 =~ /$ins1/ ) {

                    # 1 in 2
                    $mSeq   .= $ins1;
                    $qSeq   .= $insByIndex{$K1}{$i}[1];
                    $cigars .= 'I' x length($ins1);
                    $insErr++;
                } elsif ( $ins1 =~ /$ins2/ ) {

                    # 2 in 1
                    $mSeq   .= $ins2;
                    $qSeq   .= $insByIndex{$K2}{$i}[1];
                    $cigars .= 'I' x length($ins2);
                    $insErr++;
                } else {
                    $insErr++;
                }
            } elsif ( defined $insByIndex{$K1}{$i} ) {
                my $w = q{};
                if ( $i != $end ) {
                    $w = $FB2->( $i + 1 );
                } else {
                    $w = '.';
                }

                # TO-DO: can ssw permit hanging insertions?
                if ( $y ne '.' && $w ne '.' ) {
                    $insObs++;
                    $insErr++;
                } else {
                    my $ins1 = lc( $insByIndex{$K1}{$i}[0] );

                    $mSeq   .= $ins1;
                    $qSeq   .= $insByIndex{$K1}{$i}[1];
                    $cigars .= 'I' x length($ins1);
                }
            } elsif ( defined $insByIndex{$K2}{$i} ) {
                my $v = q{};
                if ( $i != $end ) {
                    $v = $FB1->( $i + 1 );
                } else {
                    $v = '.';
                }

                if ( $x ne '.' && $v ne '.' ) {
                    $insObs++;
                    $insErr++;
                } else {
                    my $ins2 = lc( $insByIndex{$K2}{$i}[0] );

                    $mSeq   .= $ins2;
                    $qSeq   .= $insByIndex{$K2}{$i}[1];
                    $cigars .= 'I' x length($ins2);
                }
            }
        }

        my $qname = $a1[5];
        my $mapq  = int( avg( $a1[6], $a2[6] ) );

        if ( !$bowtieFormat ) {
            $qname = make_merged_qname($qname);
        }
        print $OSAM $qname, "\t", '0', "\t", $REF_NAME, "\t", ( $start + 1 ), "\t", $mapq;
        print $OSAM "\t", condenseCigar($cigars), "\t*\t0\t0\t", $mSeq, "\t", $qSeq, "\n";
    } else {
        my $K = $pairs{$mID}{ $mPairs[0] }[2];
        print $OSAM $sam[$K], "\n";
    }
}

close($OSAM) or croak("Cannot close file: $OS_ERROR\n");
if ($storeStats) {
    open( my $STATS, '>', $ARGV[2] . '.stats' ) or die("Cannot open $ARGV[2].stats: $OS_ERROR\n");

    print $STATS
      "${REF_NAME}\tobs\t${obs}\n",
      "${REF_NAME}\ttmv\t${tmv}\n",
      "${REF_NAME}\tfmv\t${fmv}\n",
      "${REF_NAME}\tdmv\t${dmv}\n",
      "${REF_NAME}\tinsObs\t${insObs}\n",
      "${REF_NAME}\tinsErr\t${insErr}\n";

    close($STATS) or croak("Cannot close file: $OS_ERROR\n");
}
