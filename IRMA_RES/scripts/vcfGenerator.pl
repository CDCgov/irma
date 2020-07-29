#!/usr/bin/env perl
# Samuel Shepard - 12.2014

use Getopt::Long;
print '##fileformat=VCFv4.2',"\n";
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
print sprintf('##fileDate=%4d%02d%02d',$year+1900,$mon+1,$mday),"\n";
print '##source=',$0,' ',join(' ',@ARGV),"\n";
GetOptions(	'no-gap-allele|G' => \$noGap,
		'min-freq|F=f' => \$minFreq,
		'min-insertion-freq|I=f' => \$minFreqIns, 
		'min-deletion-freq|D=f' => \$minFreqDel, 
		'min-count|C=i' => \$minCount,
		'min-quality|Q=i' => \$minQuality,
		'min-total-col-coverage|T=i' => \$minTotal,
		'name|N' => \$name,
		'conf-not-mac-err|M=f' => \$minConf,
		'sig-level|S=f' => \$sigLevel,
		'paired-error|E=s' => \$pairedStats,
		'auto-min-freq|A' => \$autoFreq,
		'rewrite-var-table|V=s' => \$rewriteVariants
	);

print "##reference=$ARGV[0]\n";

if ( !defined($name) ) {
	$name = '.';
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
	$minTotal = 100;
}

if ( scalar(@ARGV) < 2 || scalar(@ARGV) > 4 ) {
	$message = "Usage:\n\tperl $0 <reference.fasta> <all-alleles> <insertions> <deletions>\n";
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


if ( defined($sigLevel) ) {
	$takeSig = 1;
	if ( $sigLevel >= 1 ) {
		$sigLevel /= 100;
	}

	if ( $sigLevel >= .999 ) {
		$sigLevel = '99.9';
	} elsif ( $sigLevel >= .99 ) {
		$sigLevel = '99';
	} elsif ( $sigLevel >= .95 ) {
		$sigLevel = '95';
	} elsif ( $sigLevel >= .90 ) {
		$sigLevel = '90';
	} else {
		$sigLevel = '99.9';
	}

} else {
	$takeSig = 0;
}

print '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',"\n";
print '##INFO=<ID=DPI,Number=.,Type=Integer,Description="Total Insertion Depth">',"\n";
print '##INFO=<ID=DPD,Number=.,Type=Integer,Description="Total Deletion Depth">',"\n";
print '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">',"\n";
print '##INFO=<ID=AQ,Number=A,Type=Float,Description="Average Allele Quality">',"\n";

if ( $takeSig ) {
	print '##INFO=<ID=QUB,Number=A,Type=Float,Description="Second-order corrected ',$sigLevel,'% binomial confidence upper bound on quality error estimate">',"\n";
	print '##INFO=<ID=PUB,Number=.,Type=Float,Description="Second-order corrected ',$sigLevel,'% binomial confidence upper bound on paired error estimate">',"\n";
}

if ( $minTotal > 0 ) {
	print sprintf('##FILTER=<ID=dp%d,Description="Minimum Total Coverage Depth %d">',$minTotal,$minTotal),"\n";
}

if ( $minCount > 0 ) {
	print sprintf('##FILTER=<ID=cnt%d,Description="Minimum Minority Allele Count %d">',$minCount,$minCount),"\n";
}

if ( $minQuality > 0 ) {
	print sprintf('##FILTER=<ID=q%d,Description="Minimum Minority Allele Average Quality %d">',$minQuality,$minQuality),"\n";
}


if ( $minFreqIns > 0 ) {
	print sprintf('##FILTER=<ID=fins%0f,Description="Minimum Insertion Frequency %0f">',$minFreqIns,$minFreqIns),"\n";
}

if ( $minFreqDel > 0 ) {
	print sprintf('##FILTER=<ID=fdel%0f,Description="Minimum Deletion Frequency %0f">',$minFreqDel,$minFreqDel),"\n";
}

if ( $minConf > 0 ) {
	print sprintf('##FILTER=<ID=conf%0f,Description="Minimum Confidence (Proportion not Machine Error) %0f">',$minConf,$minConf),"\n";
}

if ( $takeSig ) {
	print '##FILTER=<ID=sigP',$sigLevel,',Description="Second-order corrected ',$sigLevel,'% binomial confidence upper bound on paired error estimate">',"\n";
	print '##FILTER=<ID=sigQ',$sigLevel,',Description="Second-order corrected ',$sigLevel,'% binomial confidence upper bound on quality error estimate">',"\n";
}

%majorityTable = ();
%minorityTable = ();
%deletionTable = ();
%insertionTable = ();


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
$/ = "\n";


if ( defined($rewriteVariants) ) {
	open(VARS,'>',$rewriteVariants ) or die("ERROR: cannot open $rewriteVariants for writing.\n");
	print VARS 'Reference_Name',"\t",'Position',"\t",'Total';
	print VARS "\t",'Major Allele',"\t",'Minor Allele';
	print VARS "\t",'Major Count',"\t",'Minor Count';
	print VARS "\t",'Major Frequency',"\t",'Minor Frequency';
	print VARS "\t",'Major_Average_Quality',"\t",'Minor_Average_Quality';
	print VARS "\t",'ConfidenceNotMacErr',"\t",'PairedUB',"\t",'QualityUB',"\n";
}

$heurFreq = 0;
open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
$header = <IN>;
while($line=<IN>) {
	chomp($line);
	@fields = split("\t",$line);
	#Reference_Name Position Allele Count Total Frequency Average_Quality ConfidenceNotMacErr PairedUB QualityUB Allele_Type	
	($Reference_Name,$Position,$Allele,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB,$Allele_Type) = @fields;
	
	if ( defined($autoFreq) ) {
		if ( $ConfidenceNotMacErr == 0 && $Allele ne '-' ) {
			if ( $Frequency > $heurFreq ) {
				$heurFreq = $Frequency;
			}
		}
	}

	if ( $Allele_Type eq 'Consensus' || $Allele_Type eq 'Majority' || $Allele_Type eq 'Plurality' ) {
		$majorityTable{$Position} = [@fields];
	} else {
		if ( $Allele ne '-' ) {
			$minorityTable{$Position}{$Allele} = [@fields];
		}
	}
}
close(IN);


if ( defined($autoFreq) ) {
	if ( $heurFreq > $minFreq ) {
		$minFreq = $heurFreq;
		print sprintf('##FILTER=<ID=fsnv%0f,Description="(AUTO) Minority Allele Frequency Lower Bound %0f">',$minFreq,$minFreq),"\n";
	} else {
		print sprintf('##FILTER=<ID=fsnv%0f,Description="Minimum Minority Allele Frequency %0f">',$minFreq,$minFreq),"\n";
	}
} elsif ( $minFreq > 0 ) {
	print sprintf('##FILTER=<ID=fsnv%0f,Description="Minimum Minority Allele Frequency %0f">',$minFreq,$minFreq),"\n";
}

open(IN,'<',$ARGV[2]) or die("Cannot open $ARGV[2] for reading.\n");
$header = <IN>;
while($line=<IN>) {
	chomp($line);
	@fields = split("\t",$line);	
	# Reference_Name Upstream_Position Insert Context Called Count Total Frequency Average_Quality ConfidenceNotMacErr PairedUB QualityUB
	($Reference_Name,$Upstream_Position,$Insert,$Context,$Called,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB) = @fields;
	$insertionTable{$Upstream_Position}{$Insert} = [@fields];
}
close(IN);


open(IN,'<',$ARGV[3]) or die("Cannot open $ARGV[3] for reading.\n");
$header = <IN>;
while($line=<IN>) {
	chomp($line);
	@fields = split("\t",$line);	
	# Reference_Name Upstream_Position Length Context Called Count Total Frequency PairedUB
	($Reference_Name,$Upstream_Position,$Length,$Context,$Called,$Count,$Total,$Frequency,$PairedUB) = @fields;
	$deletionTable{$Upstream_Position}{$Length} = [@fields];
}
close(IN);


%deleted = ();
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO","\n";
foreach $p ( sort { $a <=> $b } keys(%minorityTable) ) {
	%calledSNV = %calledINS = %calledDEL = ();
	$REF = $majorityTable{$p}[2];
	$AQ = $QUAL = $QUB = $PUB = $AF = $ALT = '';

	foreach $l ( keys(%{$deletionTable{$p}}) ) {
		($Reference_Name,$Upstream_Position,$Length,$Context,$Called,$Count,$Total,$Frequency,$PairedUB) = @{$deletionTable{$p}{$l}};
		if ( $minTotal > $Total ) { next; }
		if ( $minCount > $Count ) { next; }
		if ( $minFreqDel > $Frequency ) { next; }
		if ( $takeSig && ($Frequency <= $PairedUB ) ) { next; }
		$calledDEL{$l} = 1;
	}


	$REF_EXT = '';
	@dels = (); $DPd = '';
	@dels = sort { $b <=> $a } keys(%calledDEL);
	$maxDel = $dels[0];
	$delCalls = scalar(@dels);
	if ( $delCalls > 0 ) {
		$ALT = $REF.substr($REF_SEQ,$p+$dels[0],$maxDel-$dels[0]);
		$AF =  sprintf('AF=%.10f',$deletionTable{$p}{$dels[0]}[7]);
		if ( $takeSig ) {
			$PUB = sprintf('PUB=%.10f',$deletionTable{$p}{$dels[0]}[8]);	
		}
		$start = 1;
		$AQ = 'AQ=.';
		$DPd = ';DPD='.$deletionTable{$p}{$dels[0]}[6];

		for($i=1;$i<$delCalls;$i++) {
			$DPd = ','.$deletionTable{$p}{$dels[$i]}[6];
		}

		for($i=$start;$i<$delCalls;$i++) {
			$AQ .= ',.';
			$ALT .= ','.$REF.substr($REF_SEQ,$p+$dels[$i],$maxDel-$dels[$i]);
			$AF .=  sprintf(',%.10f',$deletionTable{$p}{$dels[$i]}[7]);
			if ( $takeSig ) {
				$PUB .= sprintf(',%.10f',$deletionTable{$p}{$dels[$i]}[10]);	
			}
		}
		$REF_EXT = substr($REF_SEQ,$p,$maxDel);
	}

	foreach $b ( keys(%{$minorityTable{$p}} ) ) {
		($Reference_Name,$Position,$Allele,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB,$Allele_Type) = @{$minorityTable{$p}{$b}};
		if ( $minTotal > $Total ) { next; }
		if ( $minCount > $Count ) { next; }
		if ( $minQuality > $Average_Quality ) { next; }
		if ( defined($autoFreq) ) {
			if ( $minFreq > $Frequency ) { next; }
		} else {
			if ( $minFreq >= $Frequency ) { next; }
		}
		
		if ( $minConf > $ConfidenceNotMacErr ) { next; }
		if ( $takeSig && ($Frequency <= $PairedUB || $Frequency <= $QualityUB ) ) { next; }
		$calledSNV{$b} = 1;
		if ( defined($rewriteVariants) ) {	
			($cReference_Name,$cPosition,$cAllele,$cCount,$cTotal,$cFrequency,$cAverage_Quality,
				$cConfidenceNotMacErr,$cPairedUB,$cQualityUB,$cAllele_Type) = @{$majorityTable{$p}};
			print VARS $Reference_Name,"\t",$p,"\t",$Total,"\t";
			print VARS $cAllele,"\t",$Allele,"\t",$cCount,"\t","$Count","\t";
			print VARS $cFrequency,"\t",$Frequency,"\t",$cAverage_Quality,"\t",$Average_Quality,"\t";
			print VARS $ConfidenceNotMacErr,"\t",$PairedUB,"\t",$QualityUB,"\n";
		}
	}

	@bases = ();
	@bases = sort { $a cmp $b } keys(%calledSNV);
	$snvCalls = scalar(@bases);
	if ( $snvCalls > 0 ) {

		if ( $delCalls == 0 ) {
			$AF = sprintf('AF=%.10f',$minorityTable{$p}{$bases[0]}[5]);
			$ALT = $bases[0].$REF_EXT;
			$AQ = sprintf('AQ=%.2f',$minorityTable{$p}{$bases[0]}[6]);
			if ( $takeSig ) {
				$PUB = sprintf('PUB=%.10f',$minorityTable{$p}{$bases[0]}[8]);
			}
		} else {
			$AF .= sprintf(',%.10f',$minorityTable{$p}{$bases[0]}[5]);
			$ALT .= ','.$bases[0].$REF_EXT;
			$AQ .= sprintf(',%.2f',$minorityTable{$p}{$bases[0]}[6]);
			if ( $takeSig ) {
				$PUB = sprintf(',%.10f',$minorityTable{$p}{$bases[0]}[8]);
			}
		}
		if ( $takeSig ) {
			$QUB = sprintf('QUB=%.10f',$minorityTable{$p}{$bases[0]}[9]);
		}

		for($i=1;$i<$snvCalls;$i++) {
			$ALT .= ','.$bases[$i].$REF_EXT;
			$AF .= sprintf(',%.10f',$minorityTable{$p}{$bases[$i]}[5]);
			$AQ .= sprintf(',%.2f',$minorityTable{$p}{$bases[$i]}[6]);
			if ( $takeSig ) {
				$PUB .= sprintf(',%.10f',$minorityTable{$p}{$bases[$i]}[8]);
				$QUB .= sprintf(',%.10f',$minorityTable{$p}{$bases[$i]}[9]);
			}
		}
	}

	foreach $ins ( keys(%{$insertionTable{$p}}) ) {
		($Reference_Name,$Upstream_Position,$Insert,$Context,$Called,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB) = @{$insertionTable{$p}{$ins}};
		if ( $minTotal > $Total ) { next; }
		if ( $minCount > $Count ) { next; }
		if ( $minQuality > $Average_Quality ) { next; }

		if ( $minFreqIns > $Frequency ) { next; }
		if ( $minConf > $ConfidenceNotMacErr ) { next; }
		if ( $takeSig && ($Frequency <= $PairedUB || $Frequency <= $QualityUB) ) { next; }
		$calledINS{$ins} = 1;
	}
	@inserts = (); $DPi = '';
	@inserts = sort { $a cmp $b } keys(%calledINS);
	$insCalls = scalar(@inserts);
	if ( $insCalls > 0 ) {
		if ( $snvCalls == 0 && $delCalls == 0 ) {
			$ALT = $REF.$inserts[0].$REF_EXT;
			$AQ = sprintf('AQ=%.2f',$insertionTable{$p}{$inserts[0]}[8]);
			$AF =  sprintf('AF=%.10f',$insertionTable{$p}{$inserts[0]}[7]);
			if ( $takeSig ) {
				$PUB = sprintf('PUB=%.10f',$insertionTable{$p}{$inserts[0]}[10]);	
				$QUB = sprintf('QUB=%.10f',$insertionTable{$p}{$inserts[0]}[11]);
			}
			$start = 1;	
		} else {
			$start = 0;			
		}

		$DPi = ';DPI='.$insertionTable{$p}{$inserts[0]}[6];
		for($i=1;$i<$insCalls;$i++) {
			$DPi = ','.$insertionTable{$p}{$inserts[$i]}[6];
		}

		for($i=$start;$i<$insCalls;$i++) {
			$ALT .= ','.$REF.$inserts[$i].$REF_EXT;
			$AQ .= sprintf(',%.2f',$insertionTable{$p}{$inserts[$i]}[8]);
			$AF .=  sprintf(',%.10f',$insertionTable{$p}{$inserts[$i]}[7]);
			if ( $takeSig ) {
				$PUB .= sprintf(',%.10f',$insertionTable{$p}{$inserts[$i]}[10]);	
				$QUB .= sprintf(',%.10f',$insertionTable{$p}{$inserts[$i]}[11]);
			}
		}
	}

	$totalCalls = $snvCalls + $insCalls + $delCalls;
	if ( $totalCalls > 0 ) {
		$QUAL = '.';
		$DP = 'DP='.$majorityTable{$p}[4];
		print $majorityTable{$p}[0],"\t",$majorityTable{$p}[1],"\t$name\t";
		print $REF,$REF_EXT,"\t",$ALT,"\t",$QUAL,"\t",'PASS';
		print "\t",$DP,$DPd,$DPi,';',$AF,';',$AQ;
		if ( $takeSig ) { 
			print ';',$PUB,';',$QUB;
		}
		print "\n";
	}
}

if ( defined($rewriteVariants) ) {
	close(VARS);
}
