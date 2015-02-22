#!/usr/bin/env perl
# Sam Shepard - post-processing quality alignments
# 3.2014

use Getopt::Long;
GetOptions(
		'min-freq|F=f' => \$minFreq,
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality
	);

if ( !defined($minCount) ) {
	$minCount = 2;
} elsif ( $minCount < 0 ) {
	$minCount = 0;
}

if ( !defined($minQuality) ) {
	$minQuality = 30;
} elsif($minQuality < 0 ) {
	$minQuality = 0;
}

if ( !defined($minFreq) ) {
	$minFreq = 0.005;	
} elsif ( $minFreq < 0 || $minFreq > 100 ) {
	die("ERROR: frequency must be [0,1] or [0,100].\n");
} elsif ( $minFreq >= 1 ) {
	$minFreq /= 100;
}

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
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
					$tableByRname{$rname}[$rpos]{$qBases[$qpos]}++;
					$qByRname{$rname}[$rpos]{$qBases[$qpos]} += $Q[$qpos];
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
			} else {
				die("Extended CIGAR not yet supported.\n");
			}
		}
	}
}
close(SAM);

# Statistics
open(VARS,'>',$ARGV[2].'.txt') or die("ERROR: cannot open $ARGV[2] for writing.\n");
print VARS 'Reference_Name',"\t",'Position',"\t",'Allele',"\t",'Count',"\t",'Total',"\t",'Frequency',"\t",'Average_Quality',"\n";
%totals = ();
foreach $rname ( keys(%tableByRname) ) {
	$N = $lenByRname{$rname};
	for($p=0;$p<$N;$p++) {
		$total = 0;
		@sorted = sort { $tableByRname{$rname}[$p]{$b} <=> $tableByRname{$rname}[$p]{$a} } keys(%{$tableByRname{$rname}[$p]});
		foreach $base ( @sorted ) {
			if ( $base ne 'N' ) {
				$count = $tableByRname{$rname}[$p]{$base};
				$total += $count;
			}
		}
		$totals{$rname}[$p] = $total;
	
		foreach($i=1; $i<scalar(@sorted);$i++) {
			$base = $sorted[$i];
			$count = $tableByRname{$rname}[$p]{$base};
			$freq = $count / $total;
			if ( $base ne '-' ) {
				$quality = ($qByRname{$rname}[$p]{$base} - $count*33)/$count;
			} else {
				$quality = $minQuality;
			}

			if ( $base ne 'N' && $base ne '-' ) {
				if ( $freq < $minFreq || $count < $minCount || $quality < $minQuality ) {
					$maskAllele{$rname}[$p]{$base} = 1;
				}
				print VARS $rname,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality,"\n";
			}
		}
	}
}
close(VARS);


open(SAM,'>',$ARGV[2].'.sam') or die("ERROR: cannot open $ARGV[2].sam for writing.\n");
$numLines = scalar(@sam);
print SAM $sam[0],"\n";
print SAM $sam[1],"\n";
for($l=2;$l<$numLines;$l++) {
	$line = $sam[$l];
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = @fields = split("\t",$line);
	if ( defined($seqByRname{$rname}) ) {
		$N = $lenByRname{$rname};
		@qBases = split('',uc($seq));
		@qualities = unpack("c* i*",$qual);
		@cigars = split('',$cigar);
		$rpos=$pos-1;
		$qpos=0;
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				while($inc > 0 ) {
					$base = $qBases[$qpos];
					if ( $maskAllele{$rname}[$rpos]{$base} ) {
						$qBases[$qpos] = 'N';
					}
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
			} else {
				die("Extended CIGAR not yet supported.\n");
			}
		}

		$seq = join('',@qBases);
		$fields[9] = $seq;
		print SAM join("\t",@fields),"\n";
	}
}
close(SAM);
