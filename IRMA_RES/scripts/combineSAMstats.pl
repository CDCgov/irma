#!/usr/bin/env perl
# combineSAMstats.pl
# Sam Shepard - 9.2014

use Storable;
use POSIX;
use Getopt::Long;
GetOptions(
		'name|N=s' => \$name,
		'insertion-threshold|I=f' => \$insertionThreshold,
		'deletion-threshold|D=f' => \$deletionThreshold,
		'store-stats|S=s' => \$storeStats,
		'alternative-threshold|A=f' => \$alternativeThreshold,
		'alternative-count|C=i' => \$alternativeCount
	);

if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\t$0 [options] <REF> <STAT1> <...>\n";
	$message .= "\t\t-I|--insertion-threshold <#>\tInsertion frequency where consensus is altered. Default = 0.15 or 15%.\n";
	$message .= "\t\t-I|--deletion-threshold <#>\tDeletion frequency where consensus is altered. Default = 0.75 or 75%.\n";
	$message .= "\t\t-A|--alternative-threshold <#>\tFrequency where alternative reference allele is changed. Default is off.\n";
	$message .= "\t\t-C|--alternative-count <#>\tVariant count where alternative reference allele is changed. Default is off.\n";
	$message .= "\t\t-N|--name <STR>\t\t\tName of consensus sequence.\n";
	die($message."\n");
}

open(REF,'<',$ARGV[0]) or die("Cannot open REF $ARGV[0] for reading.\n");
$/ = ">";
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$REF_SEQ = join('',@lines);
	if ( length($REF_SEQ) < 1 ) {
		next;
	}
	$N = length($REF_SEQ);
	last;
}
close(REF);
if ( !defined($N) ) { die("ERROR: no reference found in $ARGV[0].\n"); }

if ( !defined($insertionThreshold) ) {
	$insertionThreshold = 0.15;
}

if ( !defined($deletionThreshold) ) {
	$deletionThreshold = 0.75;
}

if ( !defined($alternativeThreshold) ) {
	$alternativeThreshold = 2;
}

if ( !defined($alternativeCount) ) {
	$alternativeCount = LONG_MAX;
}

@bigTable = ();
@statRef = ();
%insTable = ();
for($i=1;$i<scalar(@ARGV);$i++) {
	@statRef = @{retrieve($ARGV[$i])};
	for $p ( 0 .. ($N-1) ) {
		foreach $base ( keys(%{$statRef[0][$p]}) ) {
			$bigTable[$p]{$base} += $statRef[0][$p]{$base};
		}
		if ( %{$statRef[1]{$p}} ) {
			foreach $insert ( keys(%{$statRef[1]{$p}}) ) {
				$insTable{$p}{$insert} += $statRef[1]{$p}{$insert};
			}
		}
	}
}

@cons = @totals = ();
for $p ( 0 .. ($N-1) ) {
	$total = 0;
	$con = '';
	my $max;
	foreach $allele ( keys(%{$bigTable[$p]}) ) {
		$count = $bigTable[$p]{$allele};
		$total += $count;
		if ( !defined($max) || $count > $max ) {
			$max = $count;
			$con = $allele;
		}
	}
	$cons[$p] = $con;
	$totals[$p] = $total;
}

$header = $header2 = '';
$consensus = $alternative = '';
if ( $name ) {
	$header = '>'.$name."\n";
	$header2 = $header;
} else {
	$header = ">consensus\n";
	$header2 = ">alternative\n";
}
for $p ( 0 .. ($N-1) ) {
	if ( $cons[$p] ne '-' ) {
		$consensus .= $cons[$p];
		# alternative non-gap
		@alleles = keys(%{$bigTable[$p]});
		if ( scalar(@alleles) > 1 ) {
			@sortedAlleles = sort { $bigTable[$p]{$b} <=> $bigTable[$p]{$a} } @alleles;
			if ( $sortedAlleles[1] eq '-' ) {
				$alternative .= $cons[$p];
			} else {
				$altCount = $bigTable[$p]{$sortedAlleles[1]};
				$altFreq = $altCount / $totals[$p];
				
				if ( $altFreq < $alternativeThreshold || $altCount < $alternativeCount ) {
					$alternative .= $cons[$p];
				} else {
					$alternative .= $sortedAlleles[1];
				}
			}
		} else {
			$alternative .= $cons[$p];
		}
	} else {
		$freq = $bigTable[$p]{$cons[$p]} / $totals[$p];
		@alleles = keys(%{$bigTable[$p]});
		if ( $freq < $deletionThreshold && scalar(@alleles) > 1 ) {
			@sortedAlleles = sort { $bigTable[$p]{$b} <=> $bigTable[$p]{$a} } @alleles;
			$consensus .= $sortedAlleles[1];
			# alternative non-gap
			if ( scalar(@alleles) > 2 ) {
				$altCount = $bigTable[$p]{$sortedAlleles[2]};
				$altFreq = $altCount / $totals[$p];
				
				if ( $altFreq < $alternativeThreshold || $altCount < $alternativeCount ) {
					$alternative .= $sortedAlleles[1];
				} else {
					$alternative .= $sortedAlleles[2];
				}
			} else {
				$alternative .= $sortedAlleles[1];
			}
		}
	}

	if ( defined($insTable{$p}) ) {
		@sortedIns = sort { $insTable{$p}{$b} <=> $insTable{$p}{$a} } keys(%{$insTable{$p}});
		if ( $p < ($N-1) ) {
			$avgTotal = int(($totals[$p] + $totals[$p+1])/2);
		} else {
			$avgTotal = $totals[$p];
		}

		$freq = $insTable{$p}{$sortedIns[0]} / $avgTotal;
		if ( $freq >= $insertionThreshold ) {
			$consensus .= lc($sortedIns[0]);
			$alternative .= lc($sortedIns[0]);
		}
	}
}

# print the consensus sequences
print $header,$consensus,"\n";
if ( $alternative ne '' && $alternative ne $consensus ) {
	print $header2,$alternative,"\n";
}
