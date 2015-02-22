#!/usr/bin/env perl
# Sam Shepard - consensus builder with indel support
# 3.2014

use Getopt::Long;
GetOptions(	'insertion-threshold|I=f' => \$insertionThreshold
	);

if ( !defined($insertionThreshold) ) {
	$insertionThreshold = 0.15;
}

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <consensus_out>\n";
	$message .= "\t\t-I|--insertion-threshold <#>\tInsertion frequency where consensus is altered. Default = .15 or 15%.";
	die($message."\n");
}

if ( !defined($suffix) ) {
	$suffix = '.ref';
} else {
	$suffix = '.'.$suffix;
}

%seqByRname = ();
open(REF,'<',$ARGV[0]) or die("Cannot open REF $ARGV[0] for reading.\n");
$/ = ">";
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$rname = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$seqByRname{$rname} = uc($seq);
	$N = $lenByRname{$rname} = length($seq);
}
close(REF);


open(SAM,'<',$ARGV[1]) or die("Cannot open SAM $ARGV[1] for reading.\n");
$/ = "\n"; %tableByRname = (); %insByRname = ();
while($line=<SAM>) {
	chomp($line);
	if ( substr($line,0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	if ( defined($seqByRname{$rname}) ) {
		$N = $lenByRname{$rname};
		$seq = uc($seq);
		$rpos=$pos-1;
		$qpos=0;
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				while($inc > 0 ) {
					$tableByRname{$rname}[$rpos]{substr($seq,$qpos,1)}++;
					$qpos++; $rpos++; $inc--;
				}
			} elsif ( $op eq 'D' ) {
				while($inc > 0 ) {
					$tableByRname{$rname}[$rpos]{'-'}++;
					$rpos++; $inc--;
				}
			} elsif( $op eq 'I' ) {
				$insert = lc(substr($seq,$qpos,$inc));
				$insByRname{$rname}{$rpos-1}{$insert}++;
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

%totals = (); %con = ();
foreach $rname ( keys(%tableByRname) ) {
	$N = $lenByRname{$rname};

	for $p ( 0 .. ($N-1) ) {
		$total = 0;
		@sorted = sort { $tableByRname{$rname}[$p]{$b} <=> $tableByRname{$rname}[$p]{$a} } keys(%{$tableByRname{$rname}[$p]});
		$con{$rname}[$p] = $sorted[0];
		foreach $base ( @sorted ) {
			$count = $tableByRname{$rname}[$p]{$base};
			$total += $count;
		}
		$totals{$rname}[$p] = $total;
#		print STDERR ($p+1),"\t",$total,"\n";
	}
}

open(CONS,'>',$ARGV[2]) or die("Cannot open $ARGV[2] for writing.\n");
foreach $rname ( keys(%tableByRname) ) {
	$N = $lenByRname{$rname};
	print CONS '>',$rname,"\n";
	for $p ( 0 .. ($N-1) ) {
		$base = $con{$rname}[$p];
		if ( $base ne '-') {
			if (  $tableByRname{$rname}[$p]{$base} > 0 ) {
				print CONS $base;
			}
		} else {
			$freq = $tableByRname{$rname}[$p]{$cons[$p]} / $totals{$rname}[$p];
			@alleles = keys(%{$tableByRname{$rname}[$p]});
			if ( $freq < $deletionThreshold && scalar(@alleles) > 1 ) {
				@sortedAlleles = sort { $tableByRname{$rname}[$p]{$b} <=> $tableByRname{$rname}[$p]{$a} } @alleles;
				if ( $tableByRname{$rname}[$p]{$sortedAlleles[1]} > 0 ) {
					print $sortedAlleles[1];
				}
			}
		}
		
		if ( defined($insByRname{$rname}{$p}) ) {
			@sortedIns = sort { $insByRname{$rname}{$p}{$b} <=> $insByRname{$rname}{$p}{$a} } keys(%{$insByRname{$rname}{$p}});
			if ( $p < ($N-1) ) {
				$avgTotal = int(($totals{$rname}[$p] + $totals{$rname}[$p+1])/2);
			} else {
				$avgTotal = $totals{$rname}[$p];
			}

			$freq = $insByRname{$rname}{$p}{$sortedIns[0]} / $avgTotal;
			if ( $freq >= $insertionThreshold ) {
				print CONS lc($sortedIns[0]);
			}
		}
	}
	print CONS "\n";
}
close(CONS);
