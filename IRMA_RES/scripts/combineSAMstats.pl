#!/usr/bin/env perl
# combineSAMstats.pl
# Sam Shepard - 2020.10

use Storable;
use POSIX;
use Getopt::Long;
Getopt::Long::Configure('no_ignore_case');
GetOptions(
		'name|N=s' => \$name,
		'insertion-threshold|I=f' => \$insertionThreshold,
		'deletion-threshold|D=f' => \$deletionThreshold,
		'insertion-depth-threshold|i=i' => \$insertionDepthThreshold,
		'deletion-depth-threshold|d=i' => \$deletionDepthThreshold,
		'alternative-threshold|A=f' => \$alternativeThreshold,
		'alternative-count|C=i' => \$alternativeCount,
		'store-stats|S=s' => \$storeStats,
		'mark-deletions|M' => \$markDeletions,
        'min-dropout-edge-support|E=i' => \$minDropoutEdgeSupport
	);

if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\t$0 [options] <REF> <STAT1> <...>\n";
	$message .= "\t\t-N|--name <STR>\t\t\t\tName of consensus sequence.\n";
	$message .= "\t\t-I|--insertion-threshold <#>\t\tInsertion frequency where consensus is altered. Default = 0.15 or 15%.\n";
	$message .= "\t\t-D|--deletion-threshold <#>\t\tDeletion frequency where consensus is altered. Default = 0.75 or 75%.\n";
	$message .= "\t\t-i|--insertion-depth-threshold <#>\tInsertion coverage depth where consensus is edit (given frequency). Default = 1.\n";
	$message .= "\t\t-d|--deletion-depth-threshold <#>\tDeletion coverage depth where consensus is edited (given frequency). Default = 1.\n";
	$message .= "\t\t-A|--alternative-threshold <#>\t\tFrequency where alternative reference allele is changed. Default is off.\n";
	$message .= "\t\t-C|--alternative-count <#>\t\tVariant count where alternative reference allele is changed. Default is off.\n";
	$message .= "\t\t-S|--store-stats <FILE>\t\t\tSave aggregate stats to a .sto file.\n";
	$message .= "\t\t-M|--mark-deletions\t\t\tOutput '-' for deletions in the consensus instead ommitting the deleted states.\n";
	$message .= "\t\t-E|--min-dropout-edge-support\t\tMinimum site coverage depth in order mask dropouts (91+) with flanking regions (6bp) with N.\n";
	die($message."\n");
}

open(REF,'<',$ARGV[0]) or die("$0 ERROR: cannot open REF $ARGV[0] for reading.\n");
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
if ( !defined($N) ) { die("$0 ERROR: no reference found in $ARGV[0].\n"); }

# Insert if >= thresholds
# Given as >= T AND >= C
if ( !defined($insertionThreshold) ) 		{ $insertionThreshold = 0.15; }
if ( !defined($insertionDepthThreshold) ) 	{ $insertionDepthThreshold = 1; }

# Delete if >= thresholds
# Given as < T OR < C implies >= T AND >= C
if ( !defined($deletionThreshold) ) 		{ $deletionThreshold = 0.75; }
if ( !defined($deletionDepthThreshold) ) 	{ $deletionDepthThreshold = 1; }

# Alternative if >= T AND C
# Given as < T OR < C implies >= T AND >= C
if ( !defined($alternativeThreshold) ) 		{ $alternativeThreshold = 2; }
if ( !defined($alternativeCount) ) 		{ $alternativeCount = LONG_MAX; }


if ( !defined($minDropoutEdgeSupport) || $minDropoutEdgeSupport < 0 ) { 
    $minDropoutEdgeSupport = 0 
} else {
    $minDropoutEdgeSupport = int($minDropoutEdgeSupport);
}


# Flags
$markDeletions = defined($markDeletions) ? 1 : 0; 
$storeStats = defined($storeStats) ? 1 : 0;



# Aggregate data
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
for my $p ( 0 .. ($N-1) ) {
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

if ( $minDropoutEdgeSupport > 0 ) {
    my $plurality_sequence = join('',( map { $_ eq '' ? '.' : $_ } @cons ));
    while( $plurality_sequence =~ m/([ATGCNatcgn]{6}[.]{91,}[ATCGNatcgn]{6})/g ) {
        for my $p ( $-[1] .. ($+[1]-1) ) {
            if ( $totals[$p] < $minDropoutEdgeSupport ) {
                $cons[$p] = 'N';
            }
        }
    }
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
for my $p ( 0 .. ($N-1) ) {
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
	# Plurality consensus is '-'
	} else {
		$freq = $bigTable[$p]{$cons[$p]} / $totals[$p];
		@alleles = keys(%{$bigTable[$p]});
		# Ignore deletion if below threshold and there exists another allele
		if ( ($bigTable[$p]{$cons[$p]} < $deletionDepthThreshold || $freq < $deletionThreshold) && scalar(@alleles) > 1 ) {
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
		# >= Thresholds for deletion OR just the deletion allele is found
		# Skip unless we are to mark the deletion
		} elsif ( $markDeletions ) {
			# Consensus is a deletion. 
			# Note: the alternative consensus cannot differ in length, so the alternative allele is not evaluated.
			$consensus	.= $cons[$p];
			$alternative 	.= $cons[$p];
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
		if ( $freq >= $insertionThreshold && $insTable{$p}{$sortedIns[0]} >= $insertionDepthThreshold ) {
			$consensus 	.= lc($sortedIns[0]);
			$alternative 	.= lc($sortedIns[0]);
		}
	}
}

# print the consensus sequences
print $header,$consensus,"\n";
if ( $alternative ne '' && $alternative ne $consensus ) {
	print $header2,$alternative,"\n";
}
