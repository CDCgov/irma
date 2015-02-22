#!/usr/bin/env perl
# Sam Shepard -- merge SAM, call variants, create consensus, do coverage
# 3.2014

use Getopt::Long;
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'print-all-sites|P' => \$printAllAlleles,
		'name|N' => \$name,
		'conf-not-mac-err|M=f' => \$minConf,
		'sig-level|S=f' => \$sigLevel
	);

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
	$message .= "\t\t-G|--no-gap-allele\t\t\tDo not count gaps alleles as variants.\n";
	$message .= "\t\t-F|--min-freq <FLT>\t\t\tMinimum frequency for a variant to be processed. Default = 0.01.\n";
	$message .= "\t\t-C|--min-count <INT>\t\t\tMinimum count of variant. Default = 2.\n";
	$message .= "\t\t-Q|--min-quality <INT>\t\t\tMinimum average variant quality, preprocesses data. Default = 20.\n";
	$message .= "\t\t-T|--min-total-col-coverage <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t-P|--print-all-vars\t\t\tPrint all variants.\n";
	$message .= "\t\t-M|--conf-not-mac-err <FLT>\t\tConfidence not machine error allowable minimum. Default = 0.5\n";
	$message .= "\t\t-S|--sig-level <FLT>\t\t\tSignificance test (90, 95, 99, 99.9) variant is not machine error.\n";
	die($message."\n");
}

# FUNCTIONS #
sub condenseCigar($) {
	my $cig = $_[0];
	my $cigar = '';
	my $state = '';
	while( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+)/g ) {
		$state = $1;
		$cigar .= length($state);
		$cigar .= substr($state,0,1);
	}
	return $cigar;
}

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

	### NEGATIVE BINOMIAL ###
	$kappa2 = $kappa ** 2;
	$eta = $kappa2/3 + 1/6;
	$gamma1 = $kappa2*(13/18) + 17/18;
	$gamma2 = $kappa2*(1/18)+7/36;
	########################
	
	sub UB($$) {
		my $p = $_[0];
		my $N = $_[1];

		if ( $N == 0 ) {
			print STDERR "Unexpected error: 0 coverage depth.\n";
			return 0;
		}

		my $u= $p/(1-$p);
		# Let b=1, so V= u + u^2
		my $V= $u + $u**2;
		# And N + -2*eta
		my $u2= ($p*$N + $eta)/($N - 2*$eta);
		my $UB = $u2 + $kappa*sqrt($V+($gamma1*$V+$gamma2)/$N)*(1/sqrt($N));
		return max(min($UB,1),0);
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

$REisBase = qr/[ATCG]/;
$REgetMolID = qr/^(.+?)[_ ]([12]):.+/;

%seqByRname = ();
$/ = ">";
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$rn = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$seqByRname{$rn} = [split('',uc($seq))];
	$N = $lenByRname{$rn} = length($seq);
}
close(REF);

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
@sam = <SAM>; chomp(@sam);
close(SAM);

%pairs = %insByIndex = %tableByRname = ();
for($K=0;$K<scalar(@sam);$K++) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
	if ( $qname =~ $REgetMolID ) {
		$qMolID = $1;
		$qSide = $2;
	}

	if ( defined($seqByRname{$rn}) ) {
		$N = $lenByRname{$rn};
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
					$qByRname{$rn}[$rpos]{$NTs[$qpos]} += $Qint[$qpos];
					$tableByRname{$rn}[$rpos]{$NTs[$qpos]}++;
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$qAln .= ' ' x $inc;
				$aln .= '-' x $inc;
				for(1..$inc) {
					$tableByRname{$rn}[$rpos]{'-'}++;
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insByIndex{$rn}{$K}{$rpos} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
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
		$aln .= '.' x (($N)-$rpos);
		$qAln .= ' ' x (($N)-$rpos);
		$pairs{$rn}{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
	}
}

%alignmentByRname = %EEbyRef = %DEbyRef = ();
open(SAM,'>',$ARGV[2].'.sam') or die("Cannot open $ARGV[2].sam for writing.\n");
open(OUT,'>',$ARGV[2].'-pairingStats.txt') or die("Cannot open $ARGV[2]-merged.txt\n");
foreach $line ( @sam ) {
	if ( substr($line,0,1) eq '@' ) {
		print SAM $line,"\n";
	} else {
		last;
	}
}
foreach $rn ( keys(%pairs) ) {
	$dmv = $obs = $fmv = $hObs = $tmv = $TMJ = 0;
	@seqRef = @{$seqByRname{$rn}};
	$N = scalar(@seqRef);
	$IND = 0;
	foreach $mID ( keys(%{$pairs{$rn}}) ) {
		@mPairs = keys(%{$pairs{$rn}{$mID}});
		if ( scalar(@mPairs) == 2 ) {
			@a1 = @{$pairs{$rn}{$mID}{'1'}};			
			@a2 = @{$pairs{$rn}{$mID}{'2'}};
			($s1,$e1) = ($a1[3],$a1[4]);
			($s2,$e2) = ($a2[3],$a2[4]);
	
			$start = min($s1,$s2);
			$end = max($e1,$e2);

			$left = $right = $mid = '';
			$left = '.' x ($start);
			$right = '.' x ($N-$end-1);

			$mSeq = '';
			$cigars = '';
			$qSeq = '';
			
			$K1 = $a1[2];
			$K2 = $a2[2];
			@bases1 = split('',$a1[0]);
			@bases2 = split('',$a2[0]);
			@quals1 = unpack("c* i*",$a1[1]);
			@quals2 = unpack("c* i*",$a2[1]);

			for($i=$start;$i<=$end;$i++) {
				$x = $bases1[$i];
				$y = $bases2[$i];
				$qx = $quals1[$i];
				$qy = $quals2[$i];
				$r = $seqRef[$i];

				if ( defined($insByIndex{$rn}{$K1}{$i}) && defined($insByIndex{$rn}{$K2}{$i}) ) {
					$ins1 = lc($insByIndex{$rn}{$K1}{$i}[0]);
					$ins2 = lc($insByIndex{$rn}{$K2}{$i}[0]);
					if ( $ins1 eq $ins2 ) {
						$mSeq .= $ins1;
						@qIns1 = split('',$insByIndex{$rn}{$K1}{$i}[1]);	
						@qIns2 = split('',$insByIndex{$rn}{$K2}{$i}[1]);
						for($qIndex=0;$qIndex<length($ins1);$qIndex++ ) {
							$qSeq .= chr(max(ord($qIns1[$qIndex]),ord($qIns2[$qIndex])));
						}
						$cigars .= 'I' x length($ins1);
					} elsif ( $ins2 =~ /$ins1/ ) {
						# 1 in 2
						$mSeq .= $ins1;
						$qSeq .= $insByIndex{$rn}{$K1}{$i}[1];
						$cigars .= 'I' x length($ins1);
					} elsif ( $ins1 =~ /$ins2/ ) {
						# 2 in 1
						$mSeq .= $ins2;
						$qSeq .= $insByIndex{$rn}{$K2}{$i}[1];
						$cigars .= 'I' x length($ins2);
					}
				}

				if ( $x ne '.' && $y ne '.' ) {
					$obs++;
					if ( $x eq $y ) {
						if ( $x eq '-' ) {
							$tmv++;
							$cigars .= 'D';
							$mid .= '-';
							$tableByRname{$rn}[$i]{'-'}--;
						} else {
							if ( $x ne $r ) { $tmv++; }
							$cigars .= 'M';
							$mSeq .= $x;	
							$mid .= $x;
							$qSeq .= chr(max($qx,$qy));
							$qByRname{$rn}[$i]{$x} -= min($qx,$qy);
							$tableByRname{$rn}[$i]{$x}--;
						}
					} elsif( $x eq $r ) {
						$fmv++;
						if ( $y eq '-' ) { $dmv++; }
						$mSeq .= $x;
						$mid .= $x;
						$qSeq .= chr($qx);
						$cigars .= 'M';
						$tableByRname{$rn}[$i]{$y}--;
						$qByRname{$rn}[$i]{$y} -= $qy;
					} elsif( $y eq $r ) {
						$fmv++;
						if ( $x eq '-' ) { $dmv++; }
						$mSeq .= $y;
						$mid .= $y;
						$qSeq .= chr($qy);
						$cigars .= 'M';
						$tableByRname{$rn}[$i]{$x}--;
						$qByRname{$rn}[$i]{$x} -= $qx;
					} else {
						$fmv++;
						if ( $x =~ $REisBase && $y !~ $REisBase ) {
							$cigars .= 'M';
							$mSeq .= $x;
							$mid .= $x;
							$qSeq .= chr($qx);
							$tableByRname{$rn}[$i]{$y}--;
							$qByRname{$rn}[$i]{$y} -= $qy;
							if ( $y eq '-' ) { $dmv++; }
						} elsif ( $x !~ $REisBase && $y =~ $REisBase ) {
							$cigars .= 'M';
							$mSeq .= $y;
							$mid .= $y;
							$qSeq .= chr($qy);
							$tableByRname{$rn}[$i]{$x}--;
							$qByRname{$rn}[$i]{$x} -= $qx;
							if ( $x eq '-' ) { $dmv++; }
						} elsif ( $qx > ($qy+4) ) {
							$cigars .= 'M';
							$mSeq .= $x;
							$mid .= $x;
							$qSeq .= chr($qx);
							$tableByRname{$rn}[$i]{$y}--;
							$qByRname{$rn}[$i]{$y} -= $qy;
							if ( $y eq '-' ) { $dmv++; }
						} elsif ( $qy > ($qx+4) ) {
							$cigars .= 'M';
							$mSeq .= $y;
							$mid .= $y;
							$qSeq .= chr($qy);
							$tableByRname{$rn}[$i]{$x}--;
							$qByRname{$rn}[$i]{$x} -= $qx;
							if ( $x eq '-' ) { $dmv++; }
						} else {
							$cigars .= 'M';
							$mSeq .= 'N';
							$mid .= 'N';
							$qSeq .= chr(int(avg($qx,$qy)));
							$tableByRname{$rn}[$i]{$x}--;
							$tableByRname{$rn}[$i]{$y}--;
							$tableByRname{$rn}[$i]{'N'}++;
							$qByRname{$rn}[$i]{$x} -= $qx;
							$qByRname{$rn}[$i]{$y} -= $qy;
							$qByRname{$rn}[$i]{'N'} += int(avg($qx,$qy));
						}
					}
				} elsif ( $x eq '.' && $y ne '.' ) {
					if ( $y eq '-' ) {
						$mid .= '-';
						$cigars .= 'D';
					} else {
						$cigars .= 'M';	
						$mSeq .= $y;
						$mid .= $y;
						$qSeq .= chr($qy);
					}
				} elsif( $x ne '.' && $y eq '.' ) {
					if ( $x eq '-' ) {
						$cigars .= 'D';
						$mid .= '-';
					} else {
						$cigars .= 'M';	
						$mSeq .= $x;
						$mid .= $x;
						$qSeq .= chr($qx);
					}
				} else {
					$cigars .= 'N';
					$mid .= '.';
				}
			}
			
			$qname = $a1[5];
			$mapq = int(avg($a1[6],$a2[6]));
			$qname =~ s/(.+?[_ ])[12](:.+)/${1}3${2}/;
			print SAM $qname,"\t",'0',"\t",$rn,"\t",($start+1),"\t",$mapq;
			print SAM "\t",condenseCigar($cigars),"\t*\t0\t0\t",$mSeq,"\t",$qSeq,"\n";
			$alignmentByRname{$rn}[$IND] = $left.$mid.$right;
		} else {
			$K = $pairs{$rn}{$mID}{$mPairs[0]}[2];
			print SAM $sam[$K],"\n";
			$alignmentByRname{$rn}[$IND] = $pairs{$rn}{$mID}{$mPairs[0]}[0];
		}
		$IND++;
	}

	$TMJ = $obs - $fmv - $tmv;
	$hObs = $TMJ + $tmv;

	if ( $obs > 0 ) {
		print OUT $rn,"\tObservations\t",$obs,"\n";
		print OUT $rn,"\tExpectedErrorRate\t",($fmv/$obs),"\n";
		print OUT $rn,"\tMinimumExpectedVariation\t",($tmv/$obs),"\n";
		print OUT $rn,"\tMinimumDeletionErrorRate\t",($dmv/$obs),"\n";
		$EEbyRef{$rn} = ($fmv/$obs);
		$DEbyRef{$rn} = ($dmv/$obs);
	} else {
		print OUT $rn,"\tObservations\t",$obs,"\n";
		print OUT $rn,"\tExpectedErrorRate\t0\n";
		print OUT $rn,"\tMinimumExpectedVariation\t0\n";
		print OUT $rn,"\tMinimumDeletionErrorRate\t0\n";
		$EEbyRef{$rn} = 0;
		$DEbyRef{$rn} = 0;
	}
}
close(SAM);
close(OUT);


%variants = ();
open(COVG,'>',$ARGV[2].'-coverage.txt') or die("Cannot open $ARGV[2]-coverage.txt for writing.\n");
open(CONS,'>',$ARGV[2].'.fasta') or die("Cannot open $ARGV[2].fasta for writing.\n");
if ( $printAllAlleles ) {
	open(ALLA,'>',$ARGV[2].'-allAlleles.txt') or die("ERROR: cannot open $ARGV[2]-allAlleles.txt for writing.\n");
	print ALLA 'Reference_Name',"\t",'Position',"\t";
	print ALLA 'Allele',"\t",'Count',"\t",'Total',"\t",'Frequency',"\t";
	print ALLA 'Average_Quality',"\t",'ConfidenceNotMacErr';
	print ALLA "\t",'PairedUB',"\t",'QualityUB',"\t",'Allele_Type',"\n";
}
open(VARS,'>',$ARGV[2].'-variants.txt') or die("ERROR: cannot open $ARGV[2]-variants.txt for writing.\n");
print VARS 'Reference_Name',"\t",'Position',"\t",'Total';
print VARS "\t",'Major Allele',"\t",'Minor Allele';
print VARS "\t",'Major Count',"\t",'Minor Count';
print VARS "\t",'Major Frequency',"\t",'Minor Frequency';
print VARS "\t",'Major_Average_Quality',"\t",'Minor_Average_Quality';
print VARS "\t",'ConfidenceNotMacErr',"\t",'PairedUB',"\t",'QualityUB',"\n";

print COVG "Position\tCoverage\tConsensus\tAmbiguous\n";
foreach $rn ( keys(%tableByRname) ) {
	$N = $lenByRname{$rn};
	print CONS '>',$rn,"\n";

	#TO-DO consider refactoring
	for($p=0;$p<$N;$p++) {
		$consensus = '.';
		$conCount = -1;
		$total = 0;
		@bases = keys(%{$tableByRname{$rn}[$p]});
		$nAlleles = scalar(@bases);

		if ( $nAlleles == 1 ) {
			$consensus = $bases[0];
			$total = $conCount = $tableByRname{$rn}[$p]{$consensus};	
		} else {
			foreach $base ( @bases ) {
				if ( $tableByRname{$rn}[$p]{$base} > $conCount && $base ne '-' ) {
					$conCount = $tableByRname{$rn}[$p]{$base};
					$consensus = $base;
				}
				$total += $tableByRname{$rn}[$p]{$base};
			}
		}

		# Account for ambiguous (do not count as part of coverage)
		$total -= $tableByRname{$rn}[$p]{'N'};

		if ( ! defined($tableByRname{$rn}[$p]{'N'}) ) {
			print COVG ($p+1),"\t",$total,"\t",$consensus,"\t",0,"\n";
		} else {
			print COVG ($p+1),"\t",$total,"\t",$consensus,"\t",$tableByRname{$rn}[$p]{'N'},"\n";
		}
		print CONS $consensus;

		if ( $total != 0 ) { 
			$conFreq = $conCount / $total;
		} else {
			$conFreq = 0;
		}

		if ( $conCount != 0 ) {
			$conQuality = ($qByRname{$rn}[$p]{$consensus} - $conCount*33)/$conCount;
		} else {
			$conQuality = $minQuality;
		}

		foreach $base ( @bases ) {
			# majority allele
			if ( $base eq $consensus ) {
				if ( $printAllAlleles ) {
					# Please revisit
					if ( $base eq 'N' ) {
						$EE = 1/(10**($conQuality/10));
						$confidence = calcProb($conFreq,$EE);
						$quality = $conQuality;
						$pairedUB = UB($EEbyRef{$rn},$conCount);
						$qualityUB = UB($EE,$conCount);
						$total = $conCount;
					} elsif ( $base eq '-' ) {
						$quality = 'NA';
						$confidence = 'NA';
						$pairedUB = UB($DEbyRef{$rn},$total);
						$qualityUB = 0;
					} else {
						$EE = 1/(10**($conQuality/10));
						$confidence = calcProb($conFreq,$EE);
						$quality = $conQuality;
						$pairedUB = UB($EEbyRef{$rn},$total);
						$qualityUB = UB($EE,$total);
					}
					print ALLA $rn,"\t",($p+1),"\t",$base,"\t",$conCount,"\t",$total,"\t",$conFreq,"\t",$quality;
					print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Majority',"\n";
				}
			} else {
				$count = $tableByRname{$rn}[$p]{$base};
				if ( $count == 0 || $base eq 'N' ) { 
					next;
				}
				$freq = $count / $total;
				if ( $base ne '-'  && $base ne 'N' ) {
					$quality = ($qByRname{$rn}[$p]{$base} - $count*33)/$count;
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
							$pairedUB = UB($DEbyRef{$rn},$total);
							$qualityUB = 0;
						} else {
							$EE = 1/(10**($quality/10));
							$confidence = calcProb($freq,$EE);
							$pairedUB = UB($EEbyRef{$rn},$total);
							$qualityUB = UB($EE,$total);
						}

						if ( $printAllAlleles ) {
							print ALLA $rn,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
							print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
						}

						if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) {
							next;
						}

						$variants{$rn}[$p]{$base} = $freq;
						$varPos{$rn}{$p} = 1;

						print VARS $rn,"\t",($p+1),"\t",$total,"\t";
						print VARS $consensus,"\t",$base,"\t",$conCount,"\t",$count,"\t";
						print VARS $conFreq,"\t",$freq,"\t",$conQuality,"\t",$quality,"\t";
						print VARS $confidence,"\t",$pairedUB,"\t",$qualityUB,"\n";
					# any variant
					} elsif ( $printAllAlleles ) {
						if ( $base eq '-' ) {
							$quality = 'NA';
							$confidence = 'NA';
							$pairedUB = UB($DEbyRef{$rn},$total);
							$qualityUB = 0;
						} else {
							$EE = 1/(10**($quality/10));
							$confidence = calcProb($freq,$EE);
							$pairedUB = UB($EEbyRef{$rn},$total);
							$qualityUB = UB($EE,$total);
						}
						print ALLA $rn,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
						print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
					}
				}
			}
		}
	}
	print CONS "\n";
}
close(CONS);
close(COVG);
close(VARS);
close(ALL);


@varKeys = keys(%variants);
if ( scalar(@varKeys) > 0 ) {
	$writtenToFile = 0;
	open(EXP,'>',$ARGV[2]."-EXPENRD.sqm") or die("Cannot open $ARGV[2]_EXPENRD.sqm for writing.\n");
	open(JAC,'>',$ARGV[2]."-JACCARD.sqm") or die("Cannot open $ARGV[2]_JACCARD.sqm for writing.\n");
	open(MUT,'>',$ARGV[2]."-MUTUALD.sqm") or die("Cannot open $ARGV[2]_MUTUALD.sqm for writing.\n");
	foreach $ref ( @varKeys ) {
		%readPats = ();
		@sequences = @{$alignmentByRname{$ref}};
		@vars = sort { $a <=> $b } keys(%{$varPos{$ref}});

		if ( scalar(@vars) > 1 ) {
			$writtenToFile = 1;
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

							$mn1 = min($Eb1,$Fb1); $mn2 = min($Eb2,$Fb2); $mnA = min($mn2,$mn1);
							$mx1 = max($Eb1,$Fb1); $mx2 = max($Eb2,$Fb2);
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
	}
	close(MUT);
	close(JAC);
	close(EXP);
	if ( ! $writtenToFile ) {
		unlink($ARGV[2]."-EXPENRD.sqm") or die("Cannot delete $ARGV[2]_EXPENRD.sqm\n");
		unlink($ARGV[2]."-JACCARD.sqm") or die("Cannot delete $ARGV[2]_JACCARD.sqm\n");
		unlink($ARGV[2]."-MUTUALD.sqm") or die("Cannot delete $ARGV[2]_MUTUALD.sqm\n");
	}
}
