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
		'seg|S' => \$convertSeg,
		'prefix|P=s' => \$prefix,
		'min-del-freq|D=f' => \$minFreqDel,
		'min-ins-freq|I=f' => \$minFreqIns,
		'deletion-file|d=s'=> \$delFile,
		'insertion-file|i=s'=> \$insFile,
		'rewrite-coverage|c=s' => \$covgRewrite
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <reference.fasta> <variants.txt> [options]\n";
	$message .= "\t\t-C|--count <INT>\t\tMinimum total count needed (coverage).\n";
	$message .= "\t\t-F|--freq <FLT>\t\t\tMinimum needed minor variant frequency.\n";
	$message .= "\t\t-d|--deletion-file <STR>\tIRMA generated deletion table.\n";
	$message .= "\t\t-D|--min-del-freq <FLT>\t\tMinimum deletion frequency.\n";
	$message .= "\t\t-i|--insertion-file <STR>\tIRMA generated insertion table.\n";
	$message .= "\t\t-I|--min-ins-freq <FLT>\t\tMinimum insertion frequency.\n";
	$message .= "\t\t-N|--name <STR>\t\t\tName of header.\n";
	$message .= "\t\t-S|--seg\t\t\tConvert the seg name to number and add to ISA.\n";
	$message .= "\t\t-P|--prefix <STR>\t\tOutput prefix for file.\n";
	$message .= "\t\t-T|--min-total-depth <INT>\tMinimum non-ambiguous column coverage. Default = 2.\n";
	die($message."\n");
}

# DEFAULT THRESHOLDS
if ( !defined($minFreq) ) {
	$minFreq = 0.25;
}

if ( !defined($minCount) ) {
	$minCount = 2;
}

if ( !defined($prefix) ) {
	$prefix = '.';
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

if ( !defined($minTotal) || $minTotal < 0 ) {
	$minTotal = 100;
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

%segMap = 	(
		'B_PB1' => 1,
		'B_PB2' => 2,
		'A_PB2' => 1,
		'A_PB1' => 2,
		'PA' => 3,
		'HA' => 4,
		'NP' => 5,
		'NA' => 6,
		'M' => 7,
		'NS' => 8
		);

# PROCESS fasta data
open(IN, '<', $ARGV[0]) or die("Cannot open $ARGV[0].\n");
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
open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1].\n");
$/ = "\n"; $header = <IN>;
@variants = <IN>; chomp(@variants);
%validPos = %freqByAllele = ();
foreach $line ( @variants ) {
	($ref,$pos,$total,$majAllele,$allele,$majCount,$count,$majFreq,$freq,$majAQ,$aq,$con,$pairedUB,$qualityUB) = split("\t",$line);
	$allele = uc($allele);
	# necessary for a called variant to be used for mixed base call
	if ( $count >= $minCount && $freq >= $minFreq && $total >= $minTotal  ) {
		# get all minor alleles greater than thresholds
		$validPos{$pos} .= $allele;
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
	foreach	$prot ( keys(%segMap) ) {
		if ( $faHeader =~ /$prot/ ) {
			$outHdr .= '_'.$segMap{$prot};
			last;
		}
	}
}

# encode variants that are valid
foreach $pos ( keys(%validPos) ) {
	$p = $pos-1;				# base zero
	$nts = $seq[$p] . $validPos{$pos};	# major + minor

	$seq[$p] = encode($nts);
# 	TRIPLET debugging
#	if ( length($nts) > 2 ) {
#		print STDERR $pos,"\t",$inHdr,"\n";
#	}
}

if ( defined($delFile) ) {
	open(IN,'<',$delFile) or die("Cannot open $delFile for reading.\n");
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
	open(IN,'<',$insFile) or die("Cannot open $insFile for reading.\n");
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

# amended consensus is written to a file
open(OUT,'>',$prefix.'/'.$outHdr.'.fa') or die("Cannot open $outHdr.fa for writing.\n");
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



if ( $covgRewrite ) {
	$/ = "\n";
	open(IN,'<',$covgRewrite) or die("Cannot open $covgRewrite.\n");
	$header = <IN>; chomp($header);
	@coverages = <IN>; chomp(@coverages); close(IN);

	open(OUT,'>',$prefix.'/'.$outHdr.'-coverage.txt') or die("Cannot open $outHdr-coverage.txt.\n");
	print OUT $header,"\n";

	$iPos = 1; $iCon = 3; $cursor = 1; $offset = 0;
	foreach $line ( @coverages ) {
		@fields = split("\t",$line);

		if ( $fields[$iCon] eq '.' ) {
			next;
		} else {
			$p = $fields[$iPos]-1;
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
