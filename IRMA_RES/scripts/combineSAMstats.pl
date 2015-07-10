#!/usr/bin/env perl
# combineSAMstats.pl
# Sam Shepard - 9.2014

use Storable;
use Getopt::Long;
GetOptions(
		'name|N=s' => \$name,
		'insertion-threshold|I=f' => \$insertionThreshold,
		'deletion-threshold|D=f' => \$deletionThreshold,
		'store-stats|S=s' => \$storeStats
	);

if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\t$0 [options] <REF> <STAT1> <...>\n";
	$message .= "\t\t-I|--insertion-threshold <#>\tInsertion frequency where consensus is altered. Default = .15 or 15%.";
	$message .= "\t\t-I|--deletion-threshold <#>\tDeletion frequency where consensus is altered. Default = .75 or 75%.";
	$message .= "\t\t-N|--name <STR>\t\tName of consensus sequence.\n";
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
}
close(REF);
if ( !defined($N) ) { die("ERROR: no reference found in $ARGV[0].\n"); }

if ( !defined($insertionThreshold) ) {
	$insertionThreshold = 0.15;
}

if ( !defined($deletionThreshold) ) {
	$deletionThreshold = 0.75;
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
		if ( defined(%{$statRef[1]{$p}}) ) {
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

if ( $name ) {
	print '>',$name,"\n";
} else {
	print ">consensus\n";
}
for $p ( 0 .. ($N-1) ) {
	if ( $cons[$p] ne '-' ) {
		print $cons[$p];
	} else {
		$freq = $bigTable[$p]{$cons[$p]} / $totals[$p];
		@alleles = keys(%{$bigTable[$p]});
		if ( $freq < $deletionThreshold && scalar(@alleles) > 1 ) {
			@sortedAlleles = sort { $bigTable[$p]{$b} <=> $bigTable[$p]{$a} } @alleles;
			print $sortedAlleles[1];
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
			print lc($sortedIns[0]);
		}
	}
}
print "\n";
