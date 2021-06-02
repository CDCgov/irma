#!/usr/bin/env perl
# amendConsensus.pl - Samuel Shepard
# 2021-02

use strict;
use warnings;

use Getopt::Long;

my ($minCount,$minTotal,$minFreq,$name,$convertSeg,$prefix,$faHeaderSuffix,$minFreqDel,$minFreqIns,$delFile,$insFile,$covgRewrite,$replaceNotEncode,$replaceWithPhase,$minConsensusSupport,$minConsensusQuality,$outputCoverageFile,$printFinalName,$a2mReference);

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
		'min-consensus-quality=f' => \$minConsensusQuality,
		'output-coverage-file=s' => \$outputCoverageFile,
		'print-name' => \$printFinalName,
		'a2m-reference' => \$a2mReference
	);

if ( scalar(@ARGV) != 2 ) {
	my $message = "Usage:\n\tperl $0 <reference.fasta> <variants.txt> [options]\n";
	$message .= "\t\t-C|--count <INT>\t\t\tMinimum total count needed (coverage).\n";
	$message .= "\t\t-F|--freq <FLT>\t\t\t\tMinimum needed minor variant frequency.\n";
	$message .= "\t\t-d|--deletion-file <STR>\t\tIRMA generated deletion table.\n";
	$message .= "\t\t-D|--min-del-freq <FLT>\t\t\tMinimum deletion frequency.\n";
	$message .= "\t\t-i|--insertion-file <STR>\t\tIRMA generated insertion table.\n";
	$message .= "\t\t-I|--min-ins-freq <FLT>\t\t\tMinimum insertion frequency.\n";
	$message .= "\t\t-N|--name <STR>\t\t\t\tName of header and file.\n";
	$message .= "\t\t-H|--fa-header-suffix\t\t\tAdd fasta header as suffix to <name>.\n";
	$message .= "\t\t-S|--seg\t\t\t\tConvert the seg name to number and add to ISA.\n";
	$message .= "\t\t-P|--prefix <STR>\t\t\tOutput prefix for file.\n";
	$message .= "\t\t-T|--min-total-depth <INT>\t\tMinimum non-ambiguous column coverage. Default = 2.\n";
	$message .= "\t\t-c|--rewrite-coverage <FILE>\t\tRe-write the coverage file.\n";
	$message .= "\t\t-o|--output-coverage-file <FILE>\tOutput coverage file, otherwise replaces file exactly.\n";
	$message .= "\t\t--min-consensus-support <INT>\t\tMinimum plurality count to call consensus, otherwise 'N'.\n";
	$message .= "\t\t--min-consensus-quality <FLT>\t\tMinimum plurality average quality to call consensus, otherwise 'N'.\n";
	$message .= "\t\t--a2m-reference\t\t\t\tInterpret reference as an align to module fasta file.\n";
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
}

if ( !defined($minConsensusQuality) || $minConsensusQuality < 0 ) {
	$minConsensusQuality = 0;
}

# mappings of nucleotides
my %M = (
	'A'=>'A',
	'C'=>'C',
	'G'=>'G',
	'T'=>'T',
	'N'=>'N',
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
		$nts =  join('', sort( split('',$nts) ) );
		$nts =~ tr/ACGT//s;
		return $M{ $nts };
	}
}

# PROCESS fasta data
open(IN, '<', $ARGV[0]) or die("$0 ERROR: cannot open $ARGV[0].\n");
$/ = ">";
my %count = ();
my ($i,$L) = (0,0);
my @seq;
my $faHeader = '';
while(my $record = <IN> ) {
	chomp($record);
	my @lines = split(/\r\n|\n|\r/, $record);
	$faHeader = shift(@lines);
		
	my $sequence;
	if ( defined($a2mReference) ) {
		$sequence = join('',@lines);
	} else {
		$sequence = uc(join('',@lines));
	}

	if ( length($sequence) == 0 ) {
		next;
	} else {
		@seq = split('',$sequence);
		$L = length($sequence);
	}
}
close(IN);

my @a2mMap = ();	# Coordinate Map for Plurality to A2M
my @refMap = ();	# Coordinate Map for A2M to A2M Reference
my $CM; 		# Coordinate Mapper Function: Plurality to A2M
my $ext;		# File extension
if ( defined($a2mReference) ) {
	# Map 0-based plurality to 0-based A2M
	my $pos = 0;
	for my $idx ( 0 .. $#seq ) {
		if ( $seq[$idx] ne '-' && $seq[$idx] ne '.' ) {
			$a2mMap[$pos] = $idx;
			$pos++;			
		}
	}

	# Map 0-based A2M to 1-based reference position
	$pos = 0;
	for my $idx ( 0 .. $#seq ) {
		if ( $seq[$idx] !~ /[a-z.]/ ) {
			$pos++;
		}
 
		$refMap[$idx] = $pos;
	}

	# Maps table position (plurality) to index of a2m (not necessarily HMM position)
	$CM = sub {
		if ( $_[0] < 0 ) {
			print STDERR "Coordinate $_[0] out of bounds.\n";
			return $a2mMap[0];;
		} elsif ( $_[0] > $#a2mMap ) {
			print STDERR "Coordinate $_[0] out of bounds.\n";
			return $a2mMap[$#a2mMap];
		} else {
			return $a2mMap[ $_[0] ];
		}
	};

	$ext = '.a2m';
} else {
	$CM = sub { return $_[0] };
	$ext = '.fa';
}

# if custom header name
my $outHdr = 'unknown';
if ( defined($name) ) {
	$outHdr = $name;
} else {
	$outHdr = $faHeader;
}

# if fasta header has a protein name, convert to flu segment number
if ( defined($convertSeg) ) {
	my @pairs = split(',',$convertSeg);
	foreach my  $pair ( @pairs ) {
		my ($prot,$numbering) = split(':',$pair);
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

# exit early
if ( defined($printFinalName) ) {
	print STDOUT $outHdr;
	exit(0);
}


# GET called variant data
open(IN,'<',$ARGV[1]) or die("$0 ERROR: cannot open $ARGV[1].\n");
$/ = "\n";
my $header = <IN>;
my @variants = <IN>; chomp(@variants);
my %validPos = ();
my  %freqByAllele = ();
foreach my $line ( @variants ) {
	my ($ref,$pos,$total,$majAllele,$allele,$majCount,$count,$majFreq,$freq,$majAQ,$aq,$con,$pairedUB,$qualityUB,$phase) = split("\t",$line);
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


# encode variants that are valid
foreach my $pos ( keys(%validPos) ) {
	my $p = $CM->($pos-1);				# base zero

	if ( !defined($replaceNotEncode) ) {
		my $nts = $seq[$p] . $validPos{$pos};	# major + minor
		$seq[$p] = encode($nts);
	} else {
		$seq[$p] = $validPos{$pos};		# set to minor
	}
}

if ( defined($delFile) ) {
	open(IN,'<',$delFile) or die("$0 ERROR: cannot open $delFile for reading.\n");
	$header = <IN>;
	while(my $line=<IN>) {
		chomp($line);
		my ($Reference_Name,$Upstream_Position,$Length,$Context,$Called,$Count,$Total,$Frequency,$PairedUB) = split("\t",$line);
		if ( $Count >= $minCount && $Frequency >= $minFreqDel && $Total >= $minTotal  ) {
			for my $p ( $Upstream_Position .. ($Upstream_Position + $Length - 1) ) {
				# deletions merely stack
				$seq[$CM->($p)] = '-';
			}
		}
	}
	close(IN);
}

my %insertions = ();
if ( defined($insFile) ) {
	open(IN,'<',$insFile) or die("$0 ERROR: cannot open $insFile for reading.\n");
	my $header = <IN>;
	while(my $line=<IN>) {
		chomp($line);
		my ($Reference_Name,$Upstream_Position,$Insert,$Context,$Called,$Count,$Total,$Frequency,$Average_Quality,$ConfidenceNotMacErr,$PairedUB,$QualityUB) = split("\t",$line);
		my $p = $CM->($Upstream_Position - 1);
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
	my $header = <IN>; chomp($header);
	my @coverages = <IN>; chomp(@coverages); close(IN);

	my $coverage_filename = $prefix.'/'.$outHdr.'-coverage.txt';
	if ( ! defined($outputCoverageFile) ) {
		$coverage_filename = $covgRewrite;
	} else {
		#$prefix.'/'.$outHdr.'-coverage.txt';
		$coverage_filename = $outputCoverageFile;
	}
	open(OUT,'>',$coverage_filename) or die("$0 ERROR: cannot open $outHdr-coverage.txt.\n");

	my $iPos = 1;
	my $iDepth = 6;
	my $iCon = 3;
	my $iQuality = 7;

	if ( defined($a2mReference) ) {
		# M is match, D is deletion, I is insertion, X is missing, and P is padded
		print OUT $header,"\tHMM_Position\tAlignment_State\n";
		my @fields = split("\t",$coverages[0]);
		my $del_suffix = "NA\t0\t-\t0\t0\t0\t0\t";
		my $miss_suffix = "NA\t0\t.\t0\t0\t0\t0\n";
		my $gene = $fields[0]."\t";

		my %cTable = ();
		my $state = '';
		my $last = 0;
		foreach my $line ( @coverages ) {
			my @fields = split("\t",$line);
			if ( $fields[$iCon] eq '.' ) {
				next;
			}
	
			my $p = $CM->($fields[$iPos] - 1);	
			if ( $seq[$p] =~ /[a-z.]/ ) {
				$state = 'I';
				$fields[$iCon] = lc($fields[$iCon]);
			} else {
				$last = $p;
				$state = 'M';
			}
				
			if ( $fields[$iDepth] < $minConsensusSupport || $fields[$iQuality] < $minConsensusQuality ) {
				if ( 0 <= $p && $p <= $#seq ) {
					if ( $state eq 'M' ) {
						$seq[$p] = 'N';
						$fields[$iCon] = 'N';
					} else {
						$seq[$p] = 'n';
						$fields[$iCon] = 'n';
					}
				}
			}

			if ( $fields[$iCon] ne '-' && $seq[$p] ne '-' ) {
				$cTable{$p} = join("\t",(@fields,$refMap[$p],$state));
			}
		}

		for my $p ( 0 .. $#seq ) {
			my $pp = $refMap[$p];
			if ( $seq[$p] eq '-') {
				print OUT $gene,$del_suffix,$pp,"\tD\n";	
			} elsif ( $seq[$p] eq '.' ) {
				print OUT $gene,$miss_suffix,$pp,"\tX\n";
			} elsif ( defined($cTable{$p}) ) {
				print OUT $cTable{$p},"\n";
			} else {
				print STDERR "Unexpected state, missing coverage at $pp of A2M. Using missing.\n";
				print OUT $gene,$miss_suffix,$pp,"\tX\n";
			}

			# Add insertions		
			if ( defined($insertions{$p}) ) {
				my @insertedBases = split('', $insertions{$p}[0] );
				my $total = $insertions{$p}[2];
				foreach my $insert ( @insertedBases ) {
					print OUT $gene,"\tNA\t",$total,"\t",lc($insert),"\tNA\tNA\tI\t$pp\n";
				}		
			}
		}
	} else {
		print OUT $header,"\n";
		my $cursor = 1;
		foreach my $line ( @coverages ) {
			my @fields = split("\t",$line);

			if ( $fields[$iCon] eq '.' ) {
				next;
			}
			
			my $p = $fields[$iPos]-1;
			if ( $fields[$iDepth] < $minConsensusSupport || $fields[$iQuality] < $minConsensusQuality ) {
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
				my @insertedBases = split('', $insertions{$p}[0] );
				my $total = $insertions{$p}[2];
				foreach my $insert ( @insertedBases ) {
					print OUT $fields[0],"\t",($cursor++),"\t",$total,"\t",$insert,"\tNA\tNA\n";
				}		
			}
		}
	}
	close(OUT);	
}

# amended consensus is written to a file
open(OUT,'>',$prefix.'/'.$outHdr . $ext) or die("$0 ERROR: cannot open $outHdr.fa for writing.\n");
print OUT '>',$outHdr,"\n";
for my $p (0..$#seq) {
	if ( defined($a2mReference) || ($seq[$p] ne '-' && $seq[$p] ne '.') ) {
		print OUT $seq[$p];
	}

	if ( defined($insertions{$p}) ) {
		print OUT lc($insertions{$p}[0]);
	}
}
print OUT "\n";
close(OUT);
