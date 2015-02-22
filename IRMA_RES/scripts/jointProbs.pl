#!/usr/bin/env perl
# Sam Shepard -- merge SAM, call variants, create consensus, do coverage
# 3.2014
# V2

use Getopt::Long;
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'print-all-alleles|P' => \$printAllAlleles,
		'name|N' => \$name,
		'conf-not-mac-err|M=f' => \$minConf,
		'sig-level|S=f' => \$sigLevel,
		'site-list|L=s' => \$siteListStr
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam>\n";
	$message .= "\t\t-G|--no-gap-allele\t\t\tDo not count gaps alleles as variants.\n";
	$message .= "\t\t-F|--min-freq <FLT>\t\t\tMinimum frequency for a variant to be processed. Default = 0.01.\n";
	$message .= "\t\t-C|--min-count <INT>\t\t\tMinimum count of variant. Default = 2.\n";
	$message .= "\t\t-Q|--min-quality <INT>\t\t\tMinimum average variant quality, preprocesses data. Default = 20.\n";
	$message .= "\t\t-T|--min-total-col-coverage <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t-P|--print-all-alleles\t\t\tPrint all alleles.\n";
	$message .= "\t\t-M|--conf-not-mac-err <FLT>\t\tConfidence not machine error allowable minimum. Default = 0.5\n";
	$message .= "\t\t-S|--sig-level <FLT>\t\t\tSignificance test (90, 95, 99, 99.9) variant is not machine error.\n";
	$message .= "\t\t-E|--extra-stats <FILE>\t\t\tPairing stats file.\n";
	$message .= "\t\t-L|--site-list <#,#,#..>\t\t\tList of sites to output.\n";
	die($message."\n");
}


if ( !defined($siteListStr) ) {
	die("Must define -L <#,#,..>\n");
} else {
	@siteList = split(',',$siteListStr);
	for($i=0;$i< scalar(@siteList);$i++) {
		$siteList[$i]--;
	}
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

	if ( 1 ) {
		### BINOMIAL ###
		$kappa2 = $kappa ** 2;
		$b_star = -1;
		$eta = $kappa2/3 + 1/6;
		$gamma1 = $b_star*($kappa2*(13/18) + 17/18);
		$gamma2 = $kappa2*(1/18)+7/36;
		########################
		sub UB($$) {
			my $p = $_[0];
			my $N = $_[1];

			if ( $N == 0 ) {
				print STDERR "Unexpected error: 0 coverage depth.\n";
				return 0;
			}

			my $u= $p;
			# Let b=1, so V= u + u^2
			my $V= $u + $u**2;
			# And N + -2*eta
			my $u2= ($p*$N + $eta)/($N - 2*$eta*$b_star);
			my $UB = $u2 + $kappa*sqrt($V+($gamma1*$V+$gamma2)/$N)*(1/sqrt($N));
			return max(min($UB,1),0);
		}
	} else {
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
$REgetMolID = qr/(.+?)[_ ]([12]):.+/;

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

%insByIndex = %tableByRname = %alignedSequences = %alnIndex = ();
for($K=0;$K<scalar(@sam);$K++) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);

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
		} else {
			$aln = '';
		}
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1; $op=$2;
			if ( $op eq 'M' ) {
				for(1..$inc) {
					$aln .= $NTs[$qpos];
					$qByRname{$rn}[$rpos]{$NTs[$qpos]} += $Qint[$qpos];
					$tableByRname{$rn}[$rpos]{$NTs[$qpos]}++;
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$aln .= '-' x $inc;
				for(1..$inc) {
					$tableByRname{$rn}[$rpos]{'-'}++;
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insByIndex{$rn}{$K}{$rpos-1} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif( $op eq 'N' ) {
				$aln .= 'N' x $inc;
				$rpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}
		$aln .= '.' x (($N)-$rpos);
		$alignedSequences{$rn}[$alnIndex{$rn}] = $aln;
		$alnIndex{$rn}++;
	}
}
print STDERR "ALN\n";

foreach $rn ( keys(%tableByRname) ) {
	if ( !defined($pStats{$rn}{$sEE}) ) {
		$pStats{$rn}{$sEE} = 0;
	}

	if ( !defined($pStats{$rn}{$sDE}) ) {
		$pStats{$rn}{$sDE} = 0;
	}

	$N = $lenByRname{$rn};
	foreach $p ( @siteList ) {
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

		foreach $base ( @bases ) {
			$count = $tableByRname{$rn}[$p]{$base};
			if ( $count == 0 || $base eq 'N' || $base eq '-' ) { 
				next;
			}
			$freq = $count / $total;

			$variants{$rn}[$p]{$base} = $freq;
			$varPos{$rn}{$p} = 1;
		}
	}
}

print STDERR "VARS\n";

@varKeys = keys(%variants);
if ( scalar(@varKeys) > 0 ) {
	$writtenToFile = 0;
	foreach $ref ( @varKeys ) {
		%readPats = ();
		@sequences = @{$alignedSequences{$rn}}; #@{$alignmentByRname{$ref}};
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
					$tmp = $v1 + 1;
					print $tmp,$b1;
					for($j=0;$j<$V;$j++) {
						$v2 = $vars[$j];
						@B2 = keys(%{$variants{$ref}[$v2]});
						foreach $b2 (@B2) {
							if ( $i == $j ) {
								if (  $b1 eq $b2 ) {
									print "\t1";
									next;
								} else {
									print "\t0";
									next;
								}
							}
							
							$Fb1 = $variants{$ref}[$v1]{$b1};
							$Fb2 = $variants{$ref}[$v2]{$b2};
							if ( $Fb1 == 0 || $Fb2 == 0 ) {
								print "\t0";
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
#								if ( $s1 eq $b1 ) 			{ $Eb1 += $X	}
#								if ( $s2 eq $b2 ) 			{ $Eb2 += $X	}
							}

							
							if ( $total != 0 ) {
								$Fb1b2 /= $total;
								print "\t",$Fb1b2;
							} else {
								print "\t0";
							}
						}
					}
					print "\n";
				}
			}
		}
	}
}
