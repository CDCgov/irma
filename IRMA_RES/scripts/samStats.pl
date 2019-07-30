#!/usr/bin/env perl
# samStats.pl
# Sam Shepard - 9.2014

use Storable;
use Getopt::Long;
GetOptions(	'ignore-annotation|G' => \$ignoreAnnotation, 'silence-complex-indel|S' => \$silenceBadIndels );
if ( scalar(@ARGV) != 3 ) {
	die("Usage:\n\t$0 <REF> <SAM> <OUT>\n");
}

open(REF,'<',$ARGV[0]) or die("$0 ERROR: cannot open REF $ARGV[0] for reading.\n");
$/ = ">";
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$REF_SEQ = join('',@lines);
	if ( length($REF_SEQ) < 1 ) {
		next;
	}
	$N = length($REF_SEQ);
	last;
}
close(REF);
if ( !defined($N) ) { die("$0 ERROR: no reference found in $ARGV[0].\n"); }

if ( $ignoreAnnotation && $REF_NAME =~ /^([^{]+)\{[^}]*}/  ) {
	$REF_NAME = $1;
}

$silenceBadIndels = defined($silenceBadIndels) ? 1 : 0;

open(SAM,'<',$ARGV[1]) or die("$0 ERROR: cannot open SAM $ARGV[1] for reading.\n");
$/ = "\n"; @table = ();
while($line=<SAM>) {
	chomp($line);
	if ( substr($line,0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	if ( $cigar eq '*' ) { next; } 
		
	if ( $silenceBadIndels && $cigar =~ /\d\d+[DI]\d+M+/ ) {
		if ( $cigar =~ /^\d+M(\d+[DI]\d+M){4,}+$/ ) {
			$cigar =~ s/(\d+)D/$1N/g; $cigar =~ s/(\d+)I/$1S/g;
		}
	}

	if ( $ignoreAnnotation && $rname =~ /^([^{]+)\{[^}]*}/  ) { $rname = $1; }

	if ( $REF_NAME eq $rname ) {
		$seq = uc($seq);
		$rpos=$pos-1;
		$qpos=0;
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				for( 1..$inc ) {
					$table[0][$rpos]{substr($seq,$qpos,1)}++;
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				for( 1..$inc ) {
					$table[0][$rpos]{'-'}++;
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insert = lc(substr($seq,$qpos,$inc));
				$table[1]{$rpos-1}{$insert}++;
				$qpos += $inc;
			} elsif ( $op eq 'N' ) {
				$rpos += $inc;
				next;
			} elsif( $op eq 'S' ) {
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
store(\@table, $ARGV[2]);
