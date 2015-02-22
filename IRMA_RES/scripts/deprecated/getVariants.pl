#!/usr/bin/env perl
# Sam Shepard - variant caller
# 3.2014

use Getopt::Long;
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'print-all-vars|P' => \$printAllVars,
		'merged-counts|M=s' => \$mergedCounts
	);

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
	$message .= "\t\t-G|--no-gap-allele\t\t\tDo not count gaps alleles as variants.\n";
	$message .= "\t\t-F|--min-freq <FLT>\t\t\tMinimum frequency for a variant to be processed. Default = 0.01.\n";
	$message .= "\t\t-C|--min-count <INT>\t\t\tMinimum count of variant. Default = 2.\n";
	$message .= "\t\t-Q|--min-quality <INT>\t\t\tMinimum average variant quality, preprocesses data. Default = 20.\n";
	$message .= "\t\t-T|--min-total-col-coverage <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t-P|--print-all-vars\t\t\tPrint all variants.\n";
	$message .= "\t\t-M|--merged-counts <FILE>\t\tGet counts from merge operation.\n";
	die($message."\n");
}

# FUNCTIONS #
sub calcProb($$) {
	my $w = $_[0];
	my $e = 1/(10**($_[1]/10));
	if ( $e > $w ) {
		return 0;
	} else {
		return (($w-$e)/$w);
	}
}

sub lgg($) {
	return log($_[0])/log(10);
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
	return (($_[0]+$_[1])/2);
}
##############


if ( !defined($noGap) ) {
	$noGap = 0;
}

if ( !defined($minCount) ) {
	$minCount = 2;
} elsif ( $minCount < 0 ) {
	$minCount = 0;
}

if ( !defined($minFreq) ) {
	$minFreq = 0.005;
} elsif( $minFreq < 0 ) {
	$minFreq = 0;
}

if ( !defined($minQuality) ) {
	$minQuality = 20;
} elsif($minQuality < 0 ) {
	$minQuality = 0;
}

if ( !defined($minTotal) || $minTotal < 0 ) {
	$minTotal = 2;
}


if ( defined($mergedCounts) ) {
	open(IN,'<',$mergedCounts) or die("Cannot open $mergedCounts for reading.\n");
	$/ = "\n";
	%meCounts = ();
	while($line=<IN>) {
		chomp($line);
		($rn,$key,$count) = split("\t",$line);
		$meCounts{$rn}{$key} = $count;
	}
	close(IN);
	
	%priors = ();
	foreach $rn ( keys(%meCounts) ) {
		if ( $meCounts{$rn}{'Observations'} > 0 && defined($meCounts{$rn}{'TrueMinor'}) && defined($meCounts{$rn}{'FalseMinor'}) ) {
			$sum = $meCounts{$rn}{'TrueMinor'}+$meCounts{$rn}{'FalseMinor'};
			$priors{$rn}{'tmv'} = $meCounts{$rn}{'TrueMinor'}/$sum;
			$priors{$rn}{'fmv'} = $meCounts{$rn}{'FalseMinor'}/$sum;
		} else {
			print STDERR "WARNING: bad priors, setting to 0.5\n";
			$priors{$rn}{'tmv'} = $priors{$rn}{'fmv'} = 0.5;
		}
	}
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
#$totalScore = 0;
%qByRname = %indexByRname = %tableByRname = %insByRname = ();
foreach $line (@sam) {
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	($qname,$flag,$rname,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	@Q = unpack("c* i*",$qual);
#	$totalScore += $mapq;

	if ( defined($seqByRname{$rname}) ) {
		$N = $lenByRname{$rname};
		if ( !defined($indexByRname{$rname}) ) {
			$indexByRname{$rname} = 0;
		} else {
			$indexByRname{$rname}++;
		}
		$headerByRname{$rname}[$indexByRname{$rname}] = $qname;
		@qBases = split('',uc($seq));
		@cigars = split('',$cigar);
		$rpos=$pos-1;
		$qpos=0;
		
		if ( $rpos > 0 ) {
			$aln = '.' x $rpos; 
		} else {
			$aln = '';
		}
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				while($inc > 0 ) {
					$aln .= $qBases[$qpos];
					$tableByRname{$rname}[$rpos]{$qBases[$qpos]}++;
					$qByRname{$rname}[$rpos]{$qBases[$qpos]} += $Q[$qpos];
					$qpos++; $rpos++; $inc--;
				}
			} elsif ( $op eq 'D' ) {
				while($inc > 0 ) {
					$aln .= '-';
					$tableByRname{$rname}[$rpos]{'-'}++;
					$rpos++; $inc--;
				}
			} elsif( $op eq 'I' ) {
				$insert = lc(substr($seq,$qpos,$inc));
				$insByRname{$rname}{$rpos}{$insert}++;
				#TO-DO: support insertions as well as quality averaging
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} elsif ( $op eq 'N' ) {
				$aln .= '.' x $inc;
				$qAln .= ' ' x $inc;
				$rpos += $inc;
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}

		$aln .= '.' x (($N)-$rpos);
		$sequenceByRname{$rname}[$indexByRname{$rname}] = $aln;
	}
}
close(SAM);

open(VARS,'>',$ARGV[2].'.txt') or die("ERROR: cannot open $ARGV[2] for writing.\n");
print VARS 'Reference_Name',"\t",'Position',"\t",'Allele',"\t",'Count',"\t",'Total',"\t",'Frequency',"\t",'Average_Quality',"\t",'ConfidenceNotMacErr',"\n";
%totals = (); %variants = ();
%refBases = ();
%refProbs = ();
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

		$refBases{$rname}[$p] = $sorted[0];
		$refProbs{$rname}[$p] = $tableByRname{$rname}[$p]{$sorted[0]} / $total;
	
		foreach($i=1; $i<scalar(@sorted);$i++) {
			$base = $sorted[$i];
			$count = $tableByRname{$rname}[$p]{$base};
			$freq = $count / $total;
			if ( $base ne '-' ) {
				$quality = ($qByRname{$rname}[$p]{$base} - $count*33)/$count;
			} else {
				$quality = $minQuality;
			}

			if ( $base ne 'N' && $freq >= $minFreq && $count >= $minCount && !($noGap && $base eq '-') && $quality >= $minQuality  && $total >= $minTotal ) {
				$variants{$rname}[$p]{$base} = $freq;
				$varPos{$rname}{$p} = 1;
				print VARS $rname,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality,"\t",calcProb($freq,$quality),"\n";
			} elsif ( $printAllVars && $base ne 'N' ) {
				print VARS $rname,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality,"\t",calcProb($freq,$quality),"\n";
			}
		}
	}
}
close(VARS);

open(EXP,'>',$ARGV[2]."_EXPENRD.sqm") or die("Cannot open $ARGV[2]_EXPENRD.sqm for writing.\n");
open(JAC,'>',$ARGV[2]."_JACCARD.sqm") or die("Cannot open $ARGV[2]_JACCARD.sqm for writing.\n");
open(MUT,'>',$ARGV[2]."_MUTUALD.sqm") or die("Cannot open $ARGV[2]_MUTUALD.sqm for writing.\n");
@refEnr = ();
@refs = keys(%headerByRname);
foreach $ref (@refs) {
	%readPats = ();
#	@headers = @{$headerByRname{$ref}};
	@sequences = @{$sequenceByRname{$ref}};
	@vars = sort { $a <=> $b } keys(%{$varPos{$ref}});

	for($s=0;$s<scalar(@sequences);$s++) {
		$aln = '';
		foreach $pos ( @vars ) {
			$aln .= substr($sequences[$s],$pos,1);
		}

		if ( $aln !~ /^[.N]+$/ ) {
			$readPats{$aln}++;
		}
	}

	$V = scalar(@vars);
	for($i=0;$i<$V;$i++) {
		$v1 = $vars[$i];
		for($j=0;$j<$V;$j++) {
			if ( $i == $j ) {
				$refEnr[$i][$i] = 1;
				next;
			}

			$v2 = $vars[$j];
			$b1 = $refBases{$ref}[$v1];
			$b2 = $refBases{$ref}[$v2];
			$Fb1 = $refProbs{$ref}[$v1];
			$Fb2 = $refProbs{$ref}[$v2];
			$total = $Fb1b2 = $Eb1 = $Eb2 = 0;
			foreach $pat ( keys(%readPats) ) {
				$s1 = substr($pat,$i,1);
				$s2 = substr($pat,$j,1);
				if ( $s1 eq '.' || $s2 eq '.' ) { next; }
				
				$X = $readPats{$pat};
				$total += $X;
				if ( $s1 eq $b1 && $s2 eq $b2 ) 	{ $Fb1b2 += $X;	}	
				if ( $s1 eq $b1 ) 			{ $Eb1 += $X	}
				if ( $s2 eq $b2 ) 			{ $Eb2 += $X	}
			}
			
			if ( $total == 0 ) {
				$enr = 0;
			} else {
				$Fb1b2	/= $total;
				$Eb1	/= $total;
				$Eb2	/= $total;
				$enr = $Fb1b2 / (max($Eb1,$Fb1) * max($Eb2,$Fb2));
			}
			$refEnr[$i][$j] = $enr; 
		}
	}

	for($i=0;$i<$V;$i++) {
		$v1 = $vars[$i];
		@B1 = keys(%{$variants{$ref}[$v1]});
		foreach $b1 (@B1) {
			print MUT ($v1+1),$b1;
			print JAC ($v1+1),$b1;
			print EXP ($v1+1),$b1;
			for($j=0;$j<$V;$j++) {
				$v2 = $vars[$j];
				@B2 = keys(%{$variants{$ref}[$v2]});
				foreach $b2 (@B2) {
					if ( $i == $j ) {
						if (  $b1 eq $b2 ) {
							print MUT "\t0";
							print JAC "\t0";
							print EXP "\t0";
							next;
						} else {
							print MUT "\t1";
							print JAC "\t1";
							print EXP "\t1";
							next;
						}
					}
					
					$Fb1 = $variants{$ref}[$v1]{$b1};
					$Fb2 = $variants{$ref}[$v2]{$b2};
					if ( $Fb1 == 0 || $Fb2 == 0 ) {
						print MUT "\t1";
						print JAC "\t1";
						print EXP "\t1";
						next;
					}

					$total = $Fb1b2 = $Eb1 = $Eb2 = 0;
					foreach $pat ( keys(%readPats) ) {
						$s1 = substr($pat,$i,1);
						$s2 = substr($pat,$j,1);
						if ( $s1 eq '.' || $s2 eq '.' ) { next; }
						
						$X = $readPats{$pat};
						$total += $X;
						if ( $s1 eq $b1 && $s2 eq $b2 ) 	{ $Fb1b2 += $X;	}	
						if ( $s1 eq $b1 ) 			{ $Eb1 += $X	}
						if ( $s2 eq $b2 ) 			{ $Eb2 += $X	}
					}

					if ( $total != 0 ) {
						$Fb1b2 /= $total;
						$Eb1 /= $total;
						$Eb2 /= $total;
						if ( $Fb1b2 == 0 || $Eb1 == 0 || $Eb2 == 0 ) {
							print MUT "\t1";
							print JAC "\t1";
							print EXP "\t1";
							next;
						}
					} else {
						print MUT "\t1";
						print EXP "\t1";
						print JAC "\t1";
						next;
					}

					$mn1 = min($Eb1,$Fb1);
					$mn2 = min($Eb2,$Fb2);
					$mnA = min($mn2,$mn1);
					$mx1 = max($Eb1,$Fb1);
					$mx2 = max($Eb2,$Fb2);
					$mutd = 1 - $Fb1b2**2/($mx1*$mx2);
					$jacc = 1 - $Fb1b2/($mx1+$mx2-$Fb1b2);
					$expd = 1 - (($Fb1b2*$mnA)/($mx1*$mx2));
					print JAC "\t$mutd";
					print MUT "\t$jacc";
					print EXP "\t$expd";
				}
			}
			print MUT "\n";
			print JAC "\n";
			print EXP "\n";
		}
	}
}
close(MUT);
close(JAC);
close(EXP);
#print $ARGV[1],"\t",$totalScore,"\n";

