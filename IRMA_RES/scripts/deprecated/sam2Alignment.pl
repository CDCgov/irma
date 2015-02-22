#!/usr/bin/env perl
# Sam Shepard -- print sam as fasta format
# 3.2014

use Storable;
use Getopt::Long;
GetOptions(
		'use-storable|S' => \$useStorable
	);

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
	$message .= "\t\t-O|--use-storable\tUse storable objects rather than FASTA.\n";
	die($message."\n");
}

$REisBase = qr/[ATCG]/;
$REgetMolID = qr/(.+?)[_ ]([12]):.+/;

%seqByRname = ();
$/ = ">"; $REF_NAME = $REF_SEQ = $REF_N = '';
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$REF_SEQ = [split('',uc($seq))];
	$REF_N = length($seq);
}
close(REF);

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
@sam = <SAM>; chomp(@sam);
close(SAM);

if ( $useStorable ) {
	%pairs = %insByIndex = ();
	for($K=0;$K<scalar(@sam);$K++) {
		if ( substr($sam[$K],0,1) eq '@' ) {
			next;
		}

		($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
		if ( $qname =~ $REgetMolID ) {
			$qMolID = $1;
			$qSide = $2;
		}

		if ( $REF_NAME eq $rn ) {
			@NTs = split('',uc($seq));
			@QCs = split('',$qual);
			@Qint = unpack("c* i*",$qual);
			@cigars = split('',$cigar);
			$rpos=$pos-1;
			$qpos=0;
			
			if ( $rpos > 0 ) {
				$aln = '.' x $rpos; 
				$qAln = ' ' x $rpos;
			} else {
				$aln = '';
				$qAln = '';
			}
			
			while($cigar =~ /(\d+)([MIDNSHP])/g ) {
				$inc=$1;
				$op=$2;
				if ( $op eq 'M' ) {
					for(1..$inc) {
						$qAln .= $QCs[$qpos];
						$aln .= $NTs[$qpos];
						$qpos++; $rpos++;
					}
				} elsif ( $op eq 'D' ) {
					$qAln .= ' ' x $inc;
					$aln .= '-' x $inc;
					for(1..$inc) {
						$rpos++;
					}
				} elsif( $op eq 'I' ) {
					$insByIndex{$K}{$rpos} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
					$qpos += $inc;
				} elsif( $op eq 'S' ) {
					$qpos += $inc;
				} elsif( $op eq 'N' ) {
					$aln .= 'N' x $inc;
					$qAln .= ' ' x $inc;
					$rpos += $inc;
				} elsif ( $op eq 'H' ) {
					next;
				} else {
					die("Extended CIGAR ($op) not yet supported.\n");
				}
			}
			$aln .= '.' x (($REF_N)-$rpos);
			$qAln .= ' ' x (($REF_N)-$rpos);
			$pairs{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
		}
	}
	store(\%pairs, $ARGV[2].'.aln');
	store(\%insByIndex, $ARGV[2].'.ins');
} else {
	open(FASTA,'>',$ARGV[2].'.fasta') or die("Cannot open $ARGV[2].fasta for writing.\n");
	for($K=0;$K<scalar(@sam);$K++) {
		if ( substr($sam[$K],0,1) eq '@' ) {
			next;
		}

		($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
		if ( $qname =~ $REgetMolID ) {
			$qMolID = $1;
			$qSide = $2;
		}

		if ( $REF_NAME eq $rn ) {
			@NTs = split('',uc($seq));
			@QCs = split('',$qual);
			@Qint = unpack("c* i*",$qual);
			@cigars = split('',$cigar);
			$rpos=$pos-1;
			$qpos=0;
			
			if ( $rpos > 0 ) {
				$aln = '.' x $rpos; 
				$qAln = ' ' x $rpos;
			} else {
				$aln = '';
				$qAln = '';
			}
			
			while($cigar =~ /(\d+)([MIDNSHP])/g ) {
				$inc=$1;
				$op=$2;
				if ( $op eq 'M' ) {
					for(1..$inc) {
						$qAln .= $QCs[$qpos];
						$aln .= $NTs[$qpos];
						$qpos++; $rpos++;
					}
				} elsif ( $op eq 'D' ) {
					$qAln .= ' ' x $inc;
					$aln .= '-' x $inc;
					for(1..$inc) {
						$rpos++;
					}
				} elsif( $op eq 'I' ) {
					$qpos += $inc;
				} elsif( $op eq 'S' ) {
					$qpos += $inc;
				} elsif( $op eq 'N' ) {
					$aln .= 'N' x $inc;
					$qAln .= ' ' x $inc;
					$rpos += $inc;
				} elsif ( $op eq 'H' ) {
					next;
				} else {
					die("Extended CIGAR ($op) not yet supported.\n");
				}
			}
			$aln .= '.' x (($REF_N)-$rpos);
			$qAln .= ' ' x (($REF_N)-$rpos);
			print FASTA '>',$qname,"\n",$aln,"\n";
		}
	}
}
