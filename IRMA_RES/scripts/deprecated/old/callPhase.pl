#!/usr/bin/env perl
# Sam Shepard -- call and do phasing -- 9.2014

use Storable;
use Getopt::Long;
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'print-all-sites|P' => \$printAllAlleles,
		'name|N' => \$name,
		'conf-not-mac-err|M=f' => \$minConf,
		'sig-level|S=f' => \$sigLevel,
		'paired-error|E=s' => \$pairedStats,
		'auto-min-freq|A' => \$autoFreq
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
	$message .= "\t\t-A|--auto-min-freq\t\t\tAutomatically find minimum frequency heuristic.\n";
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
	$REF_NAME = shift(@lines);
	$REF_SEQ = join('',@lines);
	if ( length($REF_SEQ) < 1 ) {
		next;
	}
	$REF_LEN = length($REF_SEQ);
}
close(REF);
if ( !defined($REF_LEN) ) { die("No reference found.\n"); }

$DE = $PE = 0;
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
}

%cTable = ();
%qTable = ();
%alignments = ();
@data = ();
for($i=2;$i<scalar(@ARGV);$i++) {
	@data = @{retrieve($ARGV[$i])};
	# combine alignments
	foreach $aln ( keys(%{$data[4]}) ) {
		$alignments{$aln} += $data[4]{$aln};
	}

	# combine allele and quality counts
	for $p (0..($REF_LEN-1)) {
		foreach $allele ( keys(%{$data[0][$p]}) ) {
			$cTable[$p]{$allele} += $data[0][$p]{$allele};
			$qTable[$p]{$allele} += $data[2][$p]{$allele};
		}
	}
}

%variants = (); $prefix = $ARGV[1];
open(COVG,'>',$prefix.'-coverage.txt') or die("Cannot open $prefix-coverage.txt for writing.\n");
open(CONS,'>',$prefix.'.fasta') or die("Cannot open $prefix.fasta for writing.\n");
if ( $printAllAlleles ) {
	open(ALLA,'>',$prefix.'-allAlleles.txt') or die("ERROR: cannot open $prefix-allAlleles.txt for writing.\n");
	print ALLA 'Reference_Name',"\t",'Position',"\t";
	print ALLA 'Allele',"\t",'Count',"\t",'Total',"\t",'Frequency',"\t";
	print ALLA 'Average_Quality',"\t",'ConfidenceNotMacErr';
	print ALLA "\t",'PairedUB',"\t",'QualityUB',"\t",'Allele_Type',"\n";
}
open(VARS,'>',$prefix.'-variants.txt') or die("ERROR: cannot open $prefix-variants.txt for writing.\n");
print VARS 'Reference_Name',"\t",'Position',"\t",'Total';
print VARS "\t",'Major Allele',"\t",'Minor Allele';
print VARS "\t",'Major Count',"\t",'Minor Count';
print VARS "\t",'Major Frequency',"\t",'Minor Frequency';
print VARS "\t",'Major_Average_Quality',"\t",'Minor_Average_Quality';
print VARS "\t",'ConfidenceNotMacErr',"\t",'PairedUB',"\t",'QualityUB',"\n";

print COVG "Position\tCoverage\tConsensus\tAmbiguous\n";
print CONS '>',$REF_NAME,"\n";

#TO-DO consider refactoring
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

	if ( ! defined($cTable[$p]{'N'}) ) {
		print COVG ($p+1),"\t",$total,"\t",$consensus,"\t",0,"\n";
	} else {
		print COVG ($p+1),"\t",$total,"\t",$consensus,"\t",$cTable[$p]{'N'},"\n";
	}
	print CONS $consensus;

	if ( $total != 0 ) { 
		$conFreq = $conCount / $total;
	} else {
		$conFreq = 0;
	}

	if ( $conCount != 0 ) {
		$conQuality = ($qTable[$p]{$consensus} - $conCount*33)/$conCount;
	} else {
		$conQuality = $minQuality;
	}

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
				} else {
					$ee = 1/(10**($conQuality/10));
					$confidence = calcProb($conFreq,$ee);
					$quality = $conQuality;
					$pairedUB = UB($PE,$total);
					$qualityUB = UB($ee,$total);
				}
				print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$conCount,"\t",$total,"\t",$conFreq,"\t",$quality;
				print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Majority',"\n";
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
					} else {
						$ee = 1/(10**($quality/10));
						$confidence = calcProb($freq,$ee);
						$pairedUB = UB($PE,$total);
						$qualityUB = UB($ee,$total);
					}

					if ( $printAllAlleles ) {
						print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
						print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
					}

					if ( $confidence < $minConf || $freq <= $pairedUB || $freq <= $qualityUB ) {
						next;
					}

					$variants[$p]{$base} = $freq;
					$varPos{$p} = 1;

					print VARS $REF_NAME,"\t",($p+1),"\t",$total,"\t";
					print VARS $consensus,"\t",$base,"\t",$conCount,"\t",$count,"\t";
					print VARS $conFreq,"\t",$freq,"\t",$conQuality,"\t",$quality,"\t";
					print VARS $confidence,"\t",$pairedUB,"\t",$qualityUB,"\n";
				# any variant
				} elsif ( $printAllAlleles ) {
					if ( $base eq '-' ) {
						$quality = 'NA';
						$confidence = 'NA';
						$pairedUB = UB($DE,$total);
						$qualityUB = 0;
					} else {
						$ee = 1/(10**($quality/10));
						$confidence = calcProb($freq,$ee);
						$pairedUB = UB($PE,$total);
						$qualityUB = UB($ee,$total);
					}
					print ALLA $REF_NAME,"\t",($p+1),"\t",$base,"\t",$count,"\t",$total,"\t",$freq,"\t",$quality;
					print ALLA "\t",$confidence,"\t",$pairedUB,"\t",$qualityUB,"\t",'Minority',"\n";
				}
			}
		}
	}
}
print CONS "\n";
close(CONS);
close(COVG);
close(VARS);
close(ALL);

if ( scalar(@variants) > 1 ) {
	open(EXP,'>',$prefix."-EXPENRD.sqm") or die("Cannot open ${prefix}_EXPENRD.sqm for writing.\n");
	open(JAC,'>',$prefix."-JACCARD.sqm") or die("Cannot open ${prefix}_JACCARD.sqm for writing.\n");
	open(MUT,'>',$prefix."-MUTUALD.sqm") or die("Cannot open ${prefix}_MUTUALD.sqm for writing.\n");
	%readPats = ();
	@vars = sort { $a <=> $b } keys(%varPos);
	foreach $sequence ( keys(%alignments) ) {
		$aln = '';
		foreach $pos ( @vars ) {
			$aln .= substr($sequence,$pos,1);
		}

		if ( $aln !~ /^[.N]+$/ ) {
			$readPats{$aln} += $alignments{$sequence};
		}
	}

	$V = scalar(@vars);
	for($i=0;$i<$V;$i++) {
		$v1 = $vars[$i];
		@B1 = keys(%{$variants[$v1]});
		foreach $b1 (@B1) {
			print MUT ($v1+1),$b1;
			print JAC ($v1+1),$b1;
			print EXP ($v1+1),$b1;
			for($j=0;$j<$V;$j++) {
				$v2 = $vars[$j];
				@B2 = keys(%{$variants[$v2]});
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
					
					$Fb1 = $variants[$v1]{$b1};
					$Fb2 = $variants[$v2]{$b2};
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
	close(MUT);
	close(JAC);
	close(EXP);
}
