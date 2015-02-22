#!/usr/bin/env perl
# Sam Shepard - post-processing quality alignments
# 3.2014

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam>\n";
	die($message."\n");
}

%seqByRname = ();
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
$/ = ">";
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$rname = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$seqByRname{$rname} = [split('',uc($seq))];
	$N = $lenByRname{$rname} = length($seq);
}
close(REF);


open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
$/ = "\n";
@sam = <SAM>; chomp(@sam);
%tableByRname = (); %insByRname = ();
%qTableByRname = ();
$totalScore = 0;
%readsByRname = ();
foreach $line ( @sam ) {
	chomp($line);
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	$readsByRname{$rname}++;
	if ( defined($seqByRname{$rname}) ) {
		$N = $lenByRname{$rname};
		@qBases = split('',uc($seq));
		@cigars = split('',$cigar);
		$rpos=$pos-1;
		$qpos=0;
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				while($inc > 0 ) {
					$tableByRname{$rname}[$rpos]{$qBases[$qpos]}++;
					$qpos++; $rpos++; $inc--;
				}
			} elsif ( $op eq 'D' ) {
				while($inc > 0 ) {
					$tableByRname{$rname}[$rpos]{'-'}++;
					$rpos++; $inc--;
				}
			} elsif( $op eq 'I' ) {
				$insert = lc(substr($seq,$qpos,$inc));
				$insByRname{$rname}{$rpos}{$insert}++;
				$qpos += $inc;
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
@counts = ();
foreach $rname ( keys(%tableByRname) ) {
	@counts = ();
	$N = $lenByRname{$rname};
	$avg = 0;
	for($p=0;$p<$N;$p++) {
		$total = 0;
		@bases = keys(%{$tableByRname{$rname}[$p]});
		foreach $base ( @bases ) {
			if ( $base =~ /[atcgATCG]/ ) {
				$counts[$p] += $tableByRname{$rname}[$p]{$base};
			}
		}
		$avg += $counts[$p];
	}
	$avg /= $N;
	print $readsByRname{$rname},"\t",$N,"\t",$rname,"\t",$avg,"\n";
}


