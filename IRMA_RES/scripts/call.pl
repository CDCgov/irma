#!/usr/bin/env perl
# Sam Shepard -- call and do phasing -- 9.2014

use Storable;
use Getopt::Long;
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-insertion-freq|I=f' => \$minFreqIns, 
		'min-deletion-freq|D=f' => \$minFreqDel, 
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'print-all-sites|P' => \$printAllAlleles,
		'conf-not-mac-err|M=f' => \$minConf,
		'sig-level|S=f' => \$sigLevel,
		'paired-error|E=s' => \$pairedStats,
		'auto-min-freq|A' => \$autoFreq,
		'call-table|B=s' => \$callTable
	);

if ( scalar(@ARGV) < 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <prefix> <aln.sto> <...>\n";
	$message .= "\t\t-G|--no-gap-allele\t\t\tDo not count gaps alleles as variants.\n";
	$message .= "\t\t-F|--min-freq <FLT>\t\t\tMinimum frequency for a variant to be processed. Default = 0.01.\n";
	$message .= "\t\t-C|--min-count <INT>\t\t\tMinimum count of variant. Default = 2.\n";
	$message .= "\t\t-Q|--min-quality <INT>\t\t\tMinimum average variant quality, preprocesses data. Default = 20.\n";
	$message .= "\t\t-T|--min-total-col-coverage <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t-P|--print-all-vars\t\t\tPrint all variants.\n";
	$message .= "\t\t-M|--conf-not-mac-err <FLT>\t\tConfidence not machine error allowable minimum. Default = 0.5\n";
	$message .= "\t\t-S|--sig-level <FLT>\t\t\tSignificance test (90, 95, 99, 99.9) variant is not machine error.\n";
	$message .= "\t\t-E|--paired-error <FILE>\t\tFile with paired error estimates.\n";
	$message .= "\t\t-B|--call-table <FILE>\t\t\tCall table specified in file, format:<POS>[TAB]<ALLELE>\n";
	$message .= "\t\t-A|--auto-min-freq\t\t\tAutomatically find minimum frequency heuristic.\n";
	die($message."\n");
}

# FUNCTIONS #
sub calcProb($$) {
	my $w = $_[0];
	my $e = $_[1];
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

sub toIndicesZero($) {
	my $aln = $_[0];
	my $coords = '';
	my $first = '';
	my $index = 0;
	my $length = 0;

	while( $aln =~ m/(\.+|[^.]+)/g ) {
		($length,$first)= (length($1),substr($1,0,1));
		if ( $first ne '.' ) {
			$coords .= $index.',';
			$index += $length;
			$coords .= ($index-1).';';
		} else {
			$index += $length;
		}
	}
	chop($coords);
	return $coords;	
}
#############

if ( defined($sigLevel) ) {
	$takeSig = 1;
	if ( $sigLevel >= 1 ) {
		$sigLevel /= 100;
	}

	if ( $sigLevel >= .999 ) {
		$kappa = 3.090232;
	} elsif ( $sigLevel >= .99 ) {
		$kappa = 2.326348;
	} elsif ( $sigLevel >= .95 ) {
		$kappa = 1.644854;
	} elsif ( $sigLevel >= .90 ) {
		$kappa = 1.281552;
	} else {
		$kappa = 3.090232;
	}

	### second order correction ###
	$kappa2 = $kappa ** 2;
	$eta = $kappa2/3 + 1/6;
	$gamma1 = $kappa2*(13/18) + 17/18;
	$gamma2 = $kappa2*(1/18)+7/36;
	###############################
	
	sub UB($$) {
		my $p = $_[0];
		my $N = $_[1];

		if ( $N <= 0 ) {
			print STDERR "Unexpected error: $N coverage depth.\n";
			return 0;
		}
		
		if ( $p == 1 ) {
			return 1;
		}
		#my $u= $p/(1-$p); 	# negative binomial
		#my $u = $p;		# binomial

		# Let b=-1, so V= u - u^2
		my $V= $p - $p**2;
		# And N + 2*eta
		my $u2= ($p*$N + $eta)/($N + 2*$eta);

#		print STDERR $p,"\t",$N,"\t",$V,"\t",$gamma2,"\t",$gamma1,"\t",$conFreq,"\t",$u2,"\n";
		my $inRoot = $V+($gamma2-$gamma1*$V)/$N;

		# N < gamma1 - 4*gamma2
		# N < kappa^2/2 + 31/18 
		if ( $inRoot < 0 ) {
			return 1;
			#return max(min($u2,1),0);
		} else {
			my $UB = $u2 + $kappa*sqrt($inRoot)/sqrt($N);
			return max(min($UB,1),0);
		}
	}
} else {
	$takeSig = 0;
}

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

if ( !defined($minFreqIns) ) {
	$minFreqIns = $minFreq;
} elsif ( $minFreqIns < 0 ) {
	$minFreqIns = 0;
}

if ( !defined($minFreqDel) ) {
	$minFreqDel = $minFreq;
} elsif ( $minFreqDel < 0 ) {
	$minFreqDel = 0;
}

if ( !defined($minConf) ) {
	$minConf = 0.5;
} elsif( $minConf < 0 ) {
	$minConf = 0;
}

if ( !defined($minQuality) ) {
	$minQuality = 20;
} elsif($minQuality < 0 ) {
	$minQuality = 0;
}

if ( !defined($minTotal) || $minTotal < 0 ) {
	$minTotal = 2;
}

%variants = (); $doCallTable = 0;
if ( defined($callTable) ) {
	$/ = "\n";
	open(CALLTBL,'<',$callTable) or die("Cannot open $callTable for reading.\n");
	while($line=<CALLTBL>) {
		chomp($line);
		($p,$base) = split("\t",$line);
		$variants{($p-1)}{uc($base)} = 0;
	}
	close(CALLTBL);

	if ( scalar(keys(%variants)) > 0 ) {
		$doCallTable = 1;
	}
}

# Consider implementing multiple references
$/ = ">";
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$REF_SEQ = join('',@lines);
	if ( length($REF_SEQ) < 1 ) {
		next;
	}
	$REF_LEN = length($REF_SEQ);
	last;
}
close(REF);
if ( !defined($REF_LEN) ) { die("No reference found.\n"); }

my ($DE, $PE, $IE, $is_paired) = (0,0,0,0);
if ( defined($pairedStats) ) {
	$/ = "\n";
	%pStats = ();
	open(PSF,'<',$pairedStats) or die("Cannot open $pairedStats for reading.\n");
	while($line=<PSF>) {
		chomp($line);
		($rn,$type,$value) = split("\t",$line);
		$pStats{$rn}{$type} = $value;
	}
	close(PSF);	
	$DE = $pStats{$REF_NAME}{'MinimumDeletionErrorRate'};
	$PE = $pStats{$REF_NAME}{'ExpectedErrorRate'};
	$IE = $pStats{$REF_NAME}{'MinimumInsertionErrorRate'};
	$is_paired = 1;
}

%icTable = ();
%iqTable = ();
%cTable = ();
%qTable = ();
%alignments = ();
@data = ();
for($i=2;$i<scalar(@ARGV);$i++) {
	@data = @{retrieve($ARGV[$i])};
	# combine alignments
	foreach $aln ( keys(%{$data[5]}) ) {
		$alignments{$aln} += $data[5]{$aln};
	}

	# combine allele and quality counts
	for $p (0..($REF_LEN-1)) {
		foreach $allele ( keys(%{$data[0][$p]}) ) {
			$cTable[$p]{$allele} += $data[0][$p]{$allele};
			$qTable[$p]{$allele} += $data[2][$p]{$allele};
		}
	}
	
	# handle insertion data
	foreach $p ( keys(%{$data[1]}) ) {
		foreach $insert ( keys(%{$data[1]{$p}}) ) {
			$icTable{$p}{$insert} += $data[1]{$p}{$insert};
			$iqTable{$p}{$insert} += $data[3]{$p}{$insert};
		}
	}

	# handle deletion data
	foreach $p ( keys(%{$data[4]}) ) {
		foreach $inc ( keys(%{$data[4]{$p}}) ) {
			$dcTable{$p}{$inc} += $data[4]{$p}{$inc};
		}
	}
}

$prefix = $ARGV[1];
open(COVG,'>',$prefix.'-coverage.txt') or die("Cannot open $prefix-coverage.txt for writing.\n");
open(CONS,'>',$prefix.'.fasta') or die("Cannot open $prefix.fasta for writing.\n");
if ( $printAllAlleles || $doCallTable ) {
	open(ALLA,'>',$prefix.'-allAlleles.txt') or die("ERROR: cannot open $prefix-allAlleles.txt for writing.\n");
	print ALLA 'Reference_Name',"\t",'Position',"\t";
	print ALLA 'Allele',"\t",'Count',"\t",'Total',"\t",'Frequency',"\t";
	print ALLA 'Average_Quality',"\t",'ConfidenceNotMacErr';
	print ALLA "\t",'PairedUB',"\t",'QualityUB',"\t",'Allele_Type',"\n";
}
open(VARS,'>',$prefix.'-variants.txt') or die("ERROR: cannot open $prefix-variants.txt for writing.\n");
print VARS 'Reference_Name',"\t",'Position',"\t",'Total';
print VARS "\t",'Consensus_Allele',"\t",'Minority_Allele';
print VARS "\t",'Consensus_Count',"\t",'Minority_Count';
print VARS "\t",'Consensus_Frequency',"\t",'Minority_Frequency';
print VARS "\t",'Consensus_Average_Quality',"\t",'Minority_Average_Quality';
print VARS "\t",'ConfidenceNotMacErr',"\t",'PairedUB',"\t",'QualityUB',"\n";

print COVG "Reference_Name\tPosition\tCoverage Depth\tConsensus\tDeletions\tAmbiguous\tConsensus_Count\tConsensus_Average_Quality\n";
print CONS '>',$REF_NAME,"\n";

$hFreq = 0; @alphabet = split('','ACGT-');

#TO-DO consider refactoring
%totals = (); $consensusSeq = '';
for($p=0;$p<$REF_LEN;$p++) {
	$consensus = '.';
	$conCount = -1;
	$total = 0;
	@bases = keys(%{$cTable[$p]});
	$nAlleles = scalar(@bases);

	if ( $nAlleles == 1 ) {
		$consensus = $bases[0];
		$total = $conCount = $cTable[$p]{$consensus};	
	} else {
		foreach $base ( @bases ) {
			if ( $cTable[$p]{$base} > $conCount && $base ne '-' ) {
				$conCount = $cTable[$p]{$base};
				$consensus = $base;
			}
			$total += $cTable[$p]{$base};
		}
	}

	# Account for ambiguous (do not count as part of coverage)
	$total -= $cTable[$p]{'N'};
	$totals[$p] = $total;

	print CONS $consensus;
	$consensusSeq .= $consensus;

	if ( $total != 0 ) { $conFreq = $conCount / $total; } else { $conFreq = 0; }
	if ( $conCount != 0 ) { $conQuality = ($qTable[$p]{$consensus} - $conCount*33)/$conCount; } else { $conQuality = $minQuality; }

	if ( defined($cTable[$p]{'-'}) ) {
		print COVG $REF_NAME,"\t",($p+1),"\t",($total-$cTable[$p]{'-'}),"\t",$consensus,"\t",$cTable[$p]{'-'};
	} else {
		print COVG $REF_NAME,"\t",($p+1),"\t",$total,"\t",$consensus,"\t",0;
	}

	if ( ! defined($cTable[$p]{'N'}) ) {
		print COVG "\t",0;
	} else {
		print COVG "\t",$cTable[$p]{'N'};
	}
	print COVG "\t",$conCount,"\t",$conQuality,"\n";

	if ( $doCallTable ) {
		if ( !defined($variants{$p}) ) { next; }
		foreach $base ( @alphabet ) {
			if ( defined($cTable[$p]{$base}) ) {
				$count = $cTable[$p]{$base};
				$freq = $count / $total;

				if ( $base eq '-' ) {
					$qualityUB = $quality = $confidence = 'NA';
					$pairedUB = UB($DE,$total);
				} else {
					$quality = ($qTable[$p]{$base} - $count*33)/$count;
					$ee = 1/(10**($quality/10));
					$confidence = calcProb($freq,$ee);
					$pairedUB = UB($PE,$total);
					$qualityUB = UB($ee,$total);
				}
			} else {
				$freq = $count = 0;
				$qualityUB = $quality = $confidence = 'NA';
				if ( $total == 0 ) {
					$pairedUB = 'NA';
				} elsif ( $base eq '-' ) {
					$pairedUB = UB($DE,$total);
				} else {
					$pairedUB = UB($PE,$total);
				}
			}

			if ( $base eq $consensus ) { $baseType = 'Consensus'; } else { $baseType = 'Minority'; }
			print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
			print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",$baseType,"\n";

			if ( defined($variants{$p}{$base}) ) {
				$variants{$p}{$base} = $freq;
				$varLine{$p}{$base} = $REF_NAME."\t".($p+1)."\t".$total."\t";
				$varLine{$p}{$base} .= $consensus."\t".$base."\t".$conCount."\t".$count."\t";
				$varLine{$p}{$base} .= $conFreq."\t".$freq."\t".$conQuality."\t".$quality."\t";
				$varLine{$p}{$base} .= $confidence."\t".$pairedUB."\t".$qualityUB."\n";
			}
		}
	} else {
		foreach $base ( @bases ) {
			# majority allele
			if ( $base eq $consensus ) {
				if ( $printAllAlleles ) {
					# Please revisit
					if ( $base eq 'N' ) {
						$ee = 1/(10**($conQuality/10));
						$confidence = calcProb($conFreq,$ee);
						$quality = $conQuality;
						$pairedUB = UB($PE,$conCount);
						$qualityUB = UB($ee,$conCount);
						$total = $conCount;
					} elsif ( $base eq '-' ) {
						$quality = 'NA';
						$confidence = 'NA';
						$pairedUB = UB($DE,$total);
						$qualityUB = 0;
						$ee = 0;
					} else {
						$ee = 1/(10**($conQuality/10));
						$confidence = calcProb($conFreq,$ee);
						$quality = $conQuality;
						$pairedUB = UB($PE,$total);
						$qualityUB = UB($ee,$total);
					}
					print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$conCount,"\t",$total,"\t",$conFreq,"\t",$quality;
					print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Consensus',"\n";
				}
			} else {
				$count = $cTable[$p]{$base};
				if ( $count == 0 || $base eq 'N' ) { 
					next;
				}
				$freq = $count / $total;
				if ( $base ne '-'  && $base ne 'N' ) {
					$quality = ($qTable[$p]{$base} - $count*33)/$count;
				} else {
					$quality = $minQuality;
				}

				# minor allele
				if ( $base ne 'N' ) {	
					# valid variant
					if ( !($noGap && $base eq '-') && $freq >= $minFreq && $count >= $minCount && $quality >= $minQuality  && $total >= $minTotal ) {
						if ( $base eq '-' ) {
							$confidence = 'NA';
							$quality = 'NA';
							$pairedUB = UB($DE,$total);
							$qualityUB = 0;
							$ee = 0;
						} else {
							$ee = 1/(10**($quality/10));
							$confidence = calcProb($freq,$ee);
							$pairedUB = UB($PE,$total);
							$qualityUB = UB($ee,$total);
						}

						if ( $freq <= $ee && $freq > $hFreq ) { $hFreq = $freq; }
						if ( $printAllAlleles ) {
							print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
							print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
						}

						if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) {
							next;
						}

						$variants{$p}{$base} = $freq;
						$varLine{$p}{$base} = $REF_NAME."\t".($p+1)."\t".$total."\t";
						$varLine{$p}{$base} .= $consensus."\t".$base."\t".$conCount."\t".$count."\t";
						$varLine{$p}{$base} .= $conFreq."\t".$freq."\t".$conQuality."\t".$quality."\t";
						$varLine{$p}{$base} .= $confidence."\t".$pairedUB."\t".$qualityUB."\n";
					# any variant
					} elsif ( $printAllAlleles ) {
						if ( $base eq '-' ) {
							$quality = 'NA';
							$confidence = 'NA';
							$pairedUB = UB($DE,$total);
							$qualityUB = 0;
							$ee = 0;
						} else {
							$ee = 1/(10**($quality/10));
							$confidence = calcProb($freq,$ee);
							$pairedUB = UB($PE,$total);
							$qualityUB = UB($ee,$total);
						}
						if ( $freq <= $ee && $freq > $hFreq ) { $hFreq = $freq; }
						print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
						print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
					}
				}
			}
		}
	}
}
print CONS "\n";
close(CONS); close(COVG);

foreach $p ( sort { $a <=> $b } keys(%varLine) ) {
	foreach $base ( sort { $varLine{$p}{$a} cmp $varLine{$p}{$b} } keys(%{$varLine{$p}}) ) {
		if ( $doCallTable ) {
			print VARS $varLine{$p}{$base};
		} else {
			if ( $autoFreq ) { 
				if ( $variants{$p}{$base} > $hFreq ) {
					print VARS $varLine{$p}{$base};
				} else {
					delete($variants{$p}{$base});	
				}
			} else {
				print VARS $varLine{$p}{$base};
			}
		}
	}
}
close(VARS);
close(ALLA);

%coordSupport = %coordList = ();
foreach $aln ( keys(%alignments) ) {
	$coordList{toIndicesZero($aln)} += $alignments{$aln};
}

foreach $listOfCoords ( keys(%coordList) ) {
	@coords = split(';',$listOfCoords);
	foreach $coord ( @coords ) {
		($start,$stop) = split(',',$coord);
		$coordSupport{$start}{$stop} += $coordList{$listOfCoords};
	}
}

%coordStops = ();
@coordStarts = sort { $a <=> $b } keys(%coordSupport);
foreach $start ( @coordStarts ) {
	$coordStops{$start} = [ sort { $b <=> $a } keys(%{$coordSupport{$start}}) ];
}

open(INSV,'>',$prefix.'-insertions.txt') or die("ERROR: cannot open $prefix-insertions.txt for writing.\n");
	print INSV "Reference_Name\tUpstream_Position\t";
	print INSV "Insert\tContext\tCalled\tCount\tTotal\tFrequency\t";
	print INSV 'Average_Quality',"\t",'ConfidenceNotMacErr';
	print INSV "\t",'PairedUB',"\t",'QualityUB',"\n";
foreach $p ( sort {$a <=> $b } keys(%icTable) ) {
	$pp = $p + 1;
	$total=0;
	foreach $start ( @coordStarts ) {
		if ( $start <= $p ) {
			foreach $stop ( @{$coordStops{$start}} ) {
				if ( $pp <= $stop ) {
					$total += $coordSupport{$start}{$stop};
				} else {
					last;
				}
			}
		} else {
			last;
		}
	}
	
	foreach $insert ( sort {$a cmp $b} keys(%{$icTable{$p}}) ) {
		$count = $icTable{$p}{$insert};
		if ( $count < $minCount ) { next; }

		$called = "TRUE";
		if ( $count > 0 ) {
			$quality = $iqTable{$p}{$insert} / $count;
		} else {
			$quality = 0;
		}

		if ( $quality < $minQuality ) { $called = "FALSE"; }

		if ( $total > 0 ) {
			$freq = $count / $total;
		} else {
			$freq = 0;
		}

		if ( $freq < $minFreqIns || $total < $minTotal ) { $called = "FALSE"; }
		$EE = 1/(10**($quality/10));
		$confidence = calcProb($freq,$EE);
		$pairedUB = UB($IE,$total);
		$qualityUB = UB($EE,$total);

		if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) { $called = "FALSE"; }
		$left = $right = '';
		if ( $p < 5 ) {
			$left = substr($consensusSeq,0,$p+1);
		} else {
			$left = substr($consensusSeq,$p-4,5);
		}

		if ( $p > ($REF_LEN-6) ) {
			$right = substr($consensusSeq,$pp,$REF_LEN-$pp);
		} else {
			$right = substr($consensusSeq,$pp,5);
		}

		print INSV $REF_NAME,"\t",($p+1),"\t",uc($insert),"\t",lc($left),uc($insert),lc($right);
		print INSV "\t",$called,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
		print INSV "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\n";
	}
}
close(INSV);


open(DELV,'>',$prefix.'-deletions.txt') or die("ERROR: cannot open $prefix-deletions.txt for writing.\n");
	print DELV "Reference_Name\tUpstream_Position\t";
	print DELV "Length\tContext\tCalled\tCount\tTotal\tFrequency\tPairedUB\n";
foreach $p ( sort { $a <=> $b } keys(%dcTable) ) {
	foreach $inc ( sort { $a <=> $b } keys(%{$dcTable{$p}}) ) {
		$count = $dcTable{$p}{$inc};
		if ( $count < $minCount ) { next; }
		$total=0; $pp=$p+$inc+1;
		$called = "TRUE";

		# get depth
		foreach $start ( @coordStarts ) {
			if ( $start <= $p ) {
				foreach $stop ( @{$coordStops{$start}} ) {
					if ( $pp <= $stop ) {
						$total += $coordSupport{$start}{$stop};
					} else {
						last;
					}
				}
			} else {
				last;
			}
		}

		if ( $total > 0 ) {
			$freq = $count / $total;
		} else {
			$freq = 0;
		}

		if ( $freq < $minFreqDel || $total < $minTotal ) { $called = "FALSE"; }

		$pairedUB = UB($DE,$total);
		if ( $freq <= $pairedUB ) { $called = "FALSE"; }

		$left = $right = '';
		if ( $p < 5 ) {
			$left = substr($consensusSeq,0,$p+1);
		} else {
			$left = substr($consensusSeq,$p-4,5);
		}

		if ( $p > ($REF_LEN-6-$inc) ) {
			$right = substr($consensusSeq,$pp,$REF_LEN-$pp);
		} else {
			$right = substr($consensusSeq,$pp,5);
		}
	
		$mid = '-' x $inc;
		print DELV $REF_NAME,"\t",($p+1),"\t",$inc,"\t",$left,$mid,$right;
		print DELV "\t",$called,"\t",$count,"\t",$total,"\t",$freq,"\t",$pairedUB,"\n";
	}
}
close(DELV);

$variantCount = 0;
foreach my $variantPosition ( keys(%variants) ) {
	$variantCount += scalar(%{$variants{$variantPosition}});
}
if ( $variantCount > 1 ) {
	$varFile = $prefix.'-vars.sto';
	$patFile = $prefix.'-pats.sto';
	store(\%variants,$varFile);

	%readPats = ();
	@vars = sort { $a <=> $b } keys(%variants);
	foreach $sequence ( keys(%alignments) ) {
		$aln = '';
		foreach $pos ( @vars ) {
			$aln .= substr($sequence,$pos,1);
		}

		if ( $aln !~ /^[.N]+$/ ) {
			$readPats{$aln} += $alignments{$sequence};
		}
	}

	store(\%readPats,$patFile);
}
