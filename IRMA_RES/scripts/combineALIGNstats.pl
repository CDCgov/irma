#!/usr/bin/env perl
# combineA2Mstats.pl
# Sam Shepard - 8.2014

use Storable;
use Getopt::Long;
GetOptions(	'name|N=s' => \$name, 'min-pad-count|M=i' => \$minPadCount, 
		'delete-by-ambiguity|A' => \$deleteByAmbig,
		 'skip-elongation|S' => \$skipExtension
	);

if ( scalar(@ARGV) < 1 ) {
	die("Usage:\t$0 <STAT> <...>\n");
}

if (!defined($minPadCount) ) {
	$minPadCount = 10;
}

if ( defined($skipExtension) ) {
	$notSkipExtension = 0;
} else {
	$notSkipExtension = 1;
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

# PROCESS fasta data
$/ = ">";
%count3 = %count5 = ();
@count = @statRef = ();
for($i=0;$i<scalar(@ARGV);$i++) {
	@statRef = @{retrieve($ARGV[$i])};
	for $p ( 0 .. (scalar(@{$statRef[0]})-1) ) {
		foreach $base ( keys(%{$statRef[0][$p]}) ) {
			$count[$p]{$base} += $statRef[0][$p]{$base};
		}
	}

	if ( $notSkipExtension ) {
		# LEADER
		foreach $p ( keys(%{$statRef[1]}) ) {
			while( ($base, $leaderCount) = each(%{$statRef[1]{$p}}) ) {
				$count5{$p}{$base} += $leaderCount;
			}
		} 

		# TRAILER
		foreach $p ( keys(%{$statRef[2]}) ) {
			while( ($base, $trailerCount) = each(%{$statRef[2]{$p}}) ) {
				$count3{$p}{$base} += $trailerCount;
			}

		} 
	}
}

if ( $name ) {
	print '>',$name,"\n";
} else {
	print ">consensus\n";
}
$seq ='';

$Ncount = scalar(@count);
$max5 = 0; $maxB = '';
while( ($base, $matchCount) = each(%{$count[0]}) ) {
	if ( $base ne '-' ) {
		if ($matchCount > $max5) {
			$max5 = $matchCount;
			$maxB = $base;
		}
	}
}
if ( $deleteByAmbig && $maxB eq '' ) {
	$seq .= 'N';
} else {
	$seq .= $maxB;
}

for($j = 1; $j < $Ncount - 1; $j++ ) {
	$max = 0; $maxB = ''; $total = 0;
	while( ($base, $matchCount) = each(%{$count[$j]}) ) {
		if ( $base ne '-' ) {
			if ($matchCount > $max) {
				$max = $matchCount;
				$maxB = $base;
			}
		}
	}
	if ( $deleteByAmbig && $maxB eq '' ) {
		$seq .= 'N';
	} else {
		$seq .= $maxB;
	}
}

$max3 = 0; $maxB = '';
while( ($base, $matchCount) = each(%{$count[$Ncount-1]}) ) {
	if ( $base ne '-' ) {
		if ($matchCount > $max3) {
			$max3 = $matchCount;
			$maxB = $base;
		}
	}
}
if ( $deleteByAmbig && $maxB eq '' ) {
	$seq .= 'N';
} else {
	$seq .= $maxB;
}


if ( $notSkipExtension ) {
	# LEADER
	$leader =''; $threshold = max($minPadCount,int($max5/4)+1);
	foreach $p ( 1 .. scalar(keys(%count5)) ) {
		$max = 0; $maxB = ''; $total = 0;
		while( ($base, $leaderCount) = each(%{$count5{"-$p"}}) ) {
			if ($leaderCount > $max) {
				$max = $leaderCount;
				$maxB = $base;
			}
			$total += $leaderCount;
		}
		
		if ( $max < $threshold ) {
			last;
		} else {
			$leader .= $maxB;
		}
	}
	$leader = reverse($leader);

	# TRAILER
	$trailer = ''; $threshold = max($minPadCount,int($max3/4)+1);
	foreach $p ( 0 .. scalar(keys(%count3)) - 1 ) {
		$max = 0; $maxB = ''; $total = 0;
		while( ($base, $trailerCount) = each(%{$count3{$p}}) ) {
			if ($trailerCount > $max) {
				$max = $trailerCount;
				$maxB = $base;
			}
			$total += $trailerCount;
		}
		
		if ( $max < $threshold ) {
			last;
		} else {
			$trailer .= $maxB;
		}
	}
	print $leader,$seq,$trailer,"\n";
} else {
	print $seq,"\n";
}
