#!/usr/bin/env perl
# samStats.pl
#
# Samuel S. Shepard - 2014
#
# Description: tabulate SAM assembly statistics

use Getopt::Long;
GetOptions( 'ignore-annotation|G' => \$ignoreAnnotation, 'silence-complex-indel|S' => \$silenceBadIndels );
if ( scalar(@ARGV) != 3 ) {
    die("Usage:\n\t$0 <REF> <SAM> <OUT>\n");
}

open( REF, '<', $ARGV[0] ) or die("$0 ERROR: cannot open REF $ARGV[0] for reading.\n");
$/ = ">";
while ( $record = <REF> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $REF_NAME = shift(@lines);
    $REF_SEQ  = join( '', @lines );
    if ( length($REF_SEQ) < 1 ) {
        next;
    }
    $N = length($REF_SEQ);
    last;
}
close(REF);
if ( !defined($N) ) { die("$0 ERROR: no reference found in $ARGV[0].\n"); }

if ( $ignoreAnnotation && $REF_NAME =~ /^([^{]+)\{[^}]*}/ ) {
    $REF_NAME = $1;
}

$silenceBadIndels = defined($silenceBadIndels) ? 1 : 0;

open( SAM, '<', $ARGV[1] ) or die("$0 ERROR: cannot open SAM $ARGV[1] for reading.\n");
$/ = "\n";

my ( @counts, @enc_qualities, %insert_counts, %insert_enc_qualities ) = ();

while ( $line = <SAM> ) {
    chomp($line);
    if ( substr( $line, 0, 1 ) eq '@' ) {
        next;
    }

    ( $qname, $flag, $rname, $pos, $mapq, $cigar, $mrnm, $mpos, $isize, $seq, $qual ) = split( "\t", $line );
    if ( $cigar eq '*' ) { next; }

    if ( $silenceBadIndels && $cigar =~ /\d\d+[DI]\d+M+/ ) {
        if ( $cigar =~ /^\d+M(\d+[DI]\d+M){4,}+$/ ) {
            $cigar =~ s/(\d+)D/$1N/g;
            $cigar =~ s/(\d+)I/$1S/g;
        }
    }

    if ( $ignoreAnnotation && $rname =~ /^([^{]+)\{[^}]*}/ ) { $rname = $1; }

    if ( $REF_NAME eq $rname ) {
        $seq  = uc($seq);
        $rpos = $pos - 1;
        $qpos = 0;

        while ( $cigar =~ /(\d+)([MIDNSHP])/g ) {
            $inc = $1;
            $op  = $2;
            if ( $op eq 'M' ) {
                for ( 1 .. $inc ) {
                    my $query_allele = substr( $seq, $qpos, 1 );
                    $counts[$rpos]{$query_allele}++;
                    $enc_qualities[$rpos]{$query_allele} += ord( substr( $qual, $qpos, 1 ) );
                    $qpos++;
                    $rpos++;
                }
            } elsif ( $op eq 'D' ) {
                for ( 1 .. $inc ) {
                    $counts[$rpos]{'-'}++;
                    $enc_qualities[$rpos]{'-'} = 0;
                    $rpos++;
                }
            } elsif ( $op eq 'I' ) {
                $insert = lc( substr( $seq, $qpos, $inc ) );
                $insert_counts{ $rpos - 1 }{$insert}++;

                # Sum the encoded insertion qualities into a 32-bit integer.
                # Note that 2^32 >> 126 * expected virus lengths. The value is
                # normalized downstream if needed
                $insert_enc_qualities{ $rpos - 1 }{$insert} += unpack( "%32C*", substr( $qual, $qpos, $inc ) );
                $qpos += $inc;
            } elsif ( $op eq 'N' ) {
                $rpos += $inc;
                next;
            } elsif ( $op eq 'S' ) {
                $qpos += $inc;
            } elsif ( $op eq 'H' ) {
                next;
            } else {
                die("Extended CIGAR ($op) not yet supported.\n");
            }
        }
    }
}
close(SAM);

my $OUT;
open( $OUT, '>', $ARGV[2] ) or die("$0 ERROR: cannot open OUT $ARGV[2] for writing.\n");
for my $rpos ( 0 .. $#counts ) {
    next if ( !defined $counts[$rpos] );
    foreach my $allele ( keys( %{ $counts[$rpos] } ) ) {
        print $OUT join( "\t", 'M', $rpos, $allele, $counts[$rpos]{$allele}, $enc_qualities[$rpos]{$allele} ), "\n";
    }
}

foreach my $rpos_upstream ( keys(%insert_counts) ) {
    foreach my $insert ( keys( %{ $insert_counts{$rpos_upstream} } ) ) {
        print $OUT join( "\t",
                         'I', $rpos_upstream, $insert,
                         $insert_counts{$rpos_upstream}{$insert},
                         $insert_enc_qualities{$rpos_upstream}{$insert} ),
          "\n";
    }
}
close($OUT);

