#!/usr/bin/env perl
# amendConsensus.pl - Samuel Shepard
# v2 - 12/16/2014

use Getopt::Long;
Getopt::Long::Configure('no_ignore_case');
GetOptions(	
       		'count|C=i' => \$minCount,
		'min-total-depth|T=i' => \$minTotal,
		'freq|F=f' => \$minFreq,
		'name|N=s' => \$name,
		'seg|S=s' => \$convertSeg,
		'prefix|P=s' => \$prefix,
		'fa-header-suffix|H' => \$faHeaderSuffix,
		'min-del-freq|D=f' => \$minFreqDel,
		'min-ins-freq|I=f' => \$minFreqIns,
		'deletion-file|d=s'=> \$delFile,
		'insertion-file|i=s'=> \$insFile,
		'rewrite-coverage|c=s' => \$covgRewrite,
		'replace-not-ambiguate|R' => \$replaceNotEncode,
		'belong-to-phase|B=i' => \$replaceWithPhase,
		'min-consensus-support=i' => \$minConsensusSupport,
		'replace-coverage-file' => \$replaceCoverageFile
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <reference.fasta> <variants.txt> [options]\n";
	$message .= "\t\t-C|--count <INT>\t\tMinimum total count needed (coverage).\n";
	$message .= "\t\t-F|--freq <FLT>\t\t\tMinimum needed minor variant frequency.\n";
	$message .= "\t\t-d|--deletion-file <STR>\tIRMA generated deletion table.\n";
	$message .= "\t\t-D|--min-del-freq <FLT>\t\tMinimum deletion frequency.\n";
	$message .= "\t\t-i|--insertion-file <STR>\tIRMA generated insertion table.\n";
	$message .= "\t\t-I|--min-ins-freq <FLT>\t\tMinimum insertion frequency.\n";
	$message .= "\t\t-N|--name <STR>\t\t\tName of header and file.\n";
	$message .= "\t\t-H|--fa-header-suffix\t\tAdd fasta header as suffix to <name>.\n";
	$message .= "\t\t-S|--seg\t\t\tConvert the seg name to number and add to ISA.\n";
	$message .= "\t\t-P|--prefix <STR>\t\tOutput prefix for file.\n";
	$message .= "\t\t-T|--min-total-depth <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t--min-consensus-support <INT>\tMinimum total depth per site to call plurality, otherwise 'N'.\n";
	$message .= "\t\t--replace-coverage-file\t\tImplied by min-consensus-support, replaces file exactly.\n";
	die($message."\n");
}

# DEFAULT THRESHOLDS
if ( !defined($minFreq) ) 		{ $minFreq = 0.25; 		}
if ( !defined($minCount) ) 		{ $minCount = 2; 		}
if ( !defined($prefix) ) 		{ $prefix = '.'; 		}
if ( defined($replaceWithPhase) )	{ $replaceNotEncode = 1;	}

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

if ( !defined($minTotal) || $minTotal < 0 ) {
	$minTotal = 100;
}

if ( !defined($minConsensusSupport) || $minConsensusSupport < 0 ) {
	$minConsensusSupport = 1;
} else {
	$replaceCoverageFile = 1;
}

# mappings of nucleotides
%M = (
	'AT' => 'W',
	'CG' => 'S',
	'AC' => 'M',
	'GT' => 'K',
	'AG' => 'R',
	'CT' => 'Y',
	'CGT' => 'B',
	'AGT' => 'D',
	'ACT' => 'H',
	'ACG' => 'V',
	'ACGT' => 'N'
);

# PROCESS fasta data
open(IN, '<', $ARGV[0]) or die("$0 ERROR: cannot open $ARGV[0].\n");
$/ = ">"; $i = 0; %count = (); $L = 0;
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$faHeader = shift(@lines);
	$sequence = uc(join('',@lines));

	if ( length($sequence) == 0 ) {
		next;
	}
	@seq = split('',$sequence);
}
close(IN);


# GET called variant data
open(IN,'<',$ARGV[1]) or die("$0 ERROR: cannot open $ARGV[1].\n");
$/ = "\n"; $header = <IN>;
@variants = <IN>; chomp(@variants);
%validPos = %freqByAllele = ();
foreach $line ( @variants ) {
	($ref,$pos,$total,$majAllele,$allele,$majCount,$count,$majFreq,$freq,$majAQ,$aq,$con,$pairedUB,$qualityUB,$phase) = split("\t",$line);
	$allele = uc($allele);
	# necessary for a called variant to be used for mixed base call
	
	if ( defined($replaceNotEncode) && defined($replaceWithPhase) && $replaceWithPhase == $phase ) {
			$validPos{$pos} = lc($allele);
	} else {	
		if ( $count >= $minCount && $freq >= $minFreq && $total >= $minTotal ) {
			# get all minor alleles greater than thresholds
			$validPos{$pos} .= $allele;
		}
	}
}
close(IN);

# if custom header name
if ( defined($name) ) {
	$outHdr = $name;
} else {
	$outHdr = $faHeader;
}

# if fasta header has a protein name, convert to flu segment number
if ( defined($convertSeg) ) {
	@pairs = split(',',$convertSeg);
	foreach $pair ( @pairs ) {
		($prot,$numbering) = split(':',$pair);
		if ( $faHeader =~ /$prot/ ) {
			$outHdr .= '_'.$numbering;
			last;
		}
	}
}

# Add fasta header as suffix
if ( defined($faHeaderSuffix) && defined($name) ) {
	$outHdr .= '-'.$faHeader;
} 

# encode variants that are valid
foreach $pos ( keys(%validPos) ) {
	$p = $pos-1;				# base zero

	if ( !defined($replaceNotEncode) ) {
		$nts = $seq[$p] . $validPos{$pos};	# major + minor
		$seq[$p] = encode($nts);
	} else {
		$seq[$p] = $validPos{$pos};		# set to minor
	}
}

if ( defined($delFile) ) {
	open(IN,'<',$delFile) or die("$0 ERROR: cannot open $delFile for reading.\n");
	$header = <IN>;
	while($line=<IN>) {
		chomp($line);
		($Reference_Name,$Upstream_Position,$Length,$Context,$Called,$Count,$Total,$Frequency,$PairedUB) = split("\t",$line);
		if ( $Count >= $minCount && $Frequency >= $minFreqDel && $Total >= $minTotal  ) {
			for $p ( $Upstream_Position .. ($Upstream_Position + $Length - 1) ) {
				# deletions merely stack
				$seq[$p] = '-';
			}
		}
	}
	close(IN);
}

%insertions = ();
if ( defined($insFile) ) {
	open(IN,'<',$insFile) or die("$0 ERROR: cannot open $insFile for reading.\n");
	$header = <IN>;
	while($line=<IN>) {
		chomp($line);
		($Reference_Name,$Upstream_Position,$Insert,$Context,$Called,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB) = split("\t",$line);
		$p = $Upstream_Position - 1;
		if ( $Count >= $minCount && $Frequency >= $minFreqIns && $Total >= $minTotal ) {
			# most frequent insertion wins
			if ( !defined($insertions{$p}) ) {
				$insertions{$p} = [$Insert,$Frequency,$Total];
			} elsif ( $Frequency > $insertions{$p}[1] ) {
				$insertions{$p} = [$Insert,$Frequency,$Total];
			}
		}
	}
	close(IN);
}

if ( $covgRewrite ) {
	$/ = "\n";
	open(IN,'<',$covgRewrite) or die("$0 ERROR: cannot open $covgRewrite.\n");
	$header = <IN>; chomp($header);
	@coverages = <IN>; chomp(@coverages); close(IN);

	my $coverage_filename = $prefix.'/'.$outHdr.'-coverage.txt';
	if ( defined($replaceCoverageFile) ) {
		$coverage_filename = $covgRewrite;
	}
	open(OUT,'>',$coverage_filename) or die("$0 ERROR: cannot open $outHdr-coverage.txt.\n");
	print OUT $header,"\n";

	$iPos = 1; $iDepth = 2; $iCon = 3; $cursor = 1; $offset = 0;
	foreach $line ( @coverages ) {
		@fields = split("\t",$line);
	}

	foreach $line ( @coverages ) {
		@fields = split("\t",$line);

		if ( $fields[$iCon] eq '.' ) {
			next;
		} else {
			$p = $fields[$iPos]-1;
		}

		if ( $fields[$iDepth] < $minConsensusSupport ) {
			if ( 0 <= $p && $p <= $#seq ) {
				$seq[$p] = 'N';
				$fields[$iCon] = 'N';
			}
		}

		if ( $fields[$iCon] ne '-' && $seq[$p] ne '-' ) {
			$fields[$iPos] = $cursor;
			print OUT join("\t",@fields),"\n";
			$cursor++;
		}

		# Add insertions		
		if ( defined($insertions{$p}) ) {
			@insertedBases = split('', $insertions{$p}[0] );
			$total = $insertions{$p}[2];
			foreach $insert ( @insertedBases ) {
				print OUT $fields[0],"\t",($cursor++),"\t",$total,"\t",$insert,"\tNA\tNA\n";
			}		
		}
	}
	close(OUT);	
}

# amended consensus is written to a file
open(OUT,'>',$prefix.'/'.$outHdr.'.fa') or die("$0 ERROR: cannot open $outHdr.fa for writing.\n");
print OUT '>',$outHdr,"\n";
for $p (0..$#seq) {
	if ( $seq[$p] ne '-' && $seq[$p] ne '.' ) {
		print OUT $seq[$p];
	}

	if ( defined($insertions{$p}) ) {
		print OUT lc($insertions{$p}[0]);
	}
}
print OUT "\n";
close(OUT);

# function ENCODE
# Accepts string of alleles.
# Returns mapped base call.
sub encode($) {
	my $nts = $_[0];
	if ( length($nts) == 1 ) {
		return $nts;
	} elsif( $nts =~ /-/ ) {
		return '?';
	} elsif( $nts =~ /[^ATCG]/ ) {
		return 'N';
	} else {
		return $M{ join('', sort( split('',$nts) ) )};
	}
}
