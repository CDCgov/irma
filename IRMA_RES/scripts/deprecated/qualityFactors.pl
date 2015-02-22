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
%qByRname = ();
$totalScore = 0;
print "Allele\tQuality\tQpos\tRpos\tStrand\tPrefix\n";
foreach $line ( @sam ) {
	chomp($line);
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	@Q = unpack("c* i*",$qual);
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
					print $qBases[$qpos],"\t",($Q[$qpos]-33),"\t",($qpos+1),"\t",($rpos+1),"\t";
					print $flag,"\t";
					if ( $qpos > 0 ) {
						print $qBases[$qpos-1],"\n";
					} else {
						print "*\n";
					}
					$qpos++; $rpos++; $inc--;
				}
			} elsif ( $op eq 'D' ) {
				while($inc > 0 ) {
					$rpos++; $inc--;
				}
			} elsif( $op eq 'I' ) {
				$qpos += $inc;
			} else {
				die("Extended CIGAR not yet supported.\n");
			}
		}
	}
}
close(SAM);
