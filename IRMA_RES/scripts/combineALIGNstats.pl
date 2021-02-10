#!/usr/bin/env perl
# combineA2Mstats.pl
# Sam Shepard - 8.2014

use Storable;
use POSIX;
use Getopt::Long;
GetOptions(	'name|N=s' => \$name,
		'min-pad-count|M=i' => \$minPadCount, 
		'delete-by-ambiguity|A' => \$deleteByAmbig,
		'skip-elongation|S' => \$skipExtension,
		'keep-deleted|K=s' => \$referenceSeq,
		'debug-mode|D' => \$debug,
		'count-alt|C=i' => \$altCount,
		'count-freq|F=f' => \$altFreq,
		'denominator|O=i' => \$denom
	);

if ( scalar(@ARGV) < 1 ) {
	die("Usage:\t$0 <STAT> <...>\n");
}

if ( !defined($denom) ) {
	$denom = 2;
} else {
	if ( $denom < 2 ) {
		$denom = 2;
	} else {
		$denom = int($denom);
	}
}

if ( defined($referenceSeq) ) {
	$/ = ">"; $keepDeleted = 1;
	open(REF,'<',$referenceSeq) or die("$0 ERROR: cannot open $referenceSeq for reading.\n");
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
	if ( !defined($REF_LEN) ) {
		print STDERR "WARNING (combineALIGNstats): reference $referenceSeq has no length, turning off keep-deleted.\n";
		$keepDeleted = 0;
	}
	@REF_SITES = split('',lc($REF_SEQ));
} else {
	$keepDeleted = 0;
}

if (!defined($minPadCount) ) {
	$minPadCount = 10;
}

if ( defined($skipExtension) ) {
	$notSkipExtension = 0;
} else {
	$notSkipExtension = 1;
}

if ( !defined($altFreq) ) {
	$altFreq = 2;
}

if ( !defined($altCount) ) {
	$altCount = LONG_MAX;
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

if ( $debug ) {
	foreach $p ( sort { $a <=> $b } keys(%count5) ) {
		print STDERR sprintf("%5d::",$p);
		foreach $base ( sort { $count5{$p}{$b} <=> $count5{$p}{$a} || $a cmp $b } keys(%{$count5{$p}}) ) {
			print STDERR "\t",$base,":",$count5{$p}{$base};
		}
		print STDERR "\n";
	} 

	print STDERR "5'\n";
	for $p ( 0 .. (scalar(@count)-1) ) {
		print STDERR sprintf("%5d::",$p);
		foreach $base ( sort { $count[$p]{$b} <=> $count[$p]{$a} || $a cmp $b } keys(%{$count[$p]}) ) {
			print STDERR "\t",$base,":",$count[$p]{$base};
		}
		print STDERR "\n";
       }

	print STDERR "3'\n";
	foreach $p ( sort { $a <=> $b } keys(%count3) ) {
		print STDERR sprintf("+%4d::",$p);
		foreach $base ( sort { $count3{$p}{$b} <=> $count3{$p}{$a} || $a cmp $b } keys(%{$count3{$p}}) ) {
			print STDERR "\t",$base,":",$count3{$p}{$base};
		}
		print STDERR "\n";
	} 

}

if ( $name ) {
	$header =  '>'.$name."\n";
	$header2 = '>'.$name.'{alt}'."\n";
} else {
	$header = ">consensus\n";
	$header2 = ">alternative\n";
}
$seqAlt = $seq ='';
$Ncount = scalar(@count);

if ( $keepDeleted && $Ncount != $REF_LEN ) {
	print STDERR "WARNING (combineALIGNstats): $Ncount != $REF_LEN, bad reference ($REF_NAME), turning off keep-deleted.\n";
	$keepDeleted = 0;
}

# middle
$max5 = $max3 = 0;
for($j = 0; $j < $Ncount; $j++ ) {
	$alt = $max = 0; $altB = $maxB = ''; $total = 0;
	while( ($base, $matchCount) = each(%{$count[$j]}) ) {
		if ( $base ne '-' ) {
			$total += $matchCount;
			if ($matchCount > $max) {
				$alt = $max;
				$altB = $maxB;

				$max = $matchCount;
				$maxB = $base;
			} elsif ( $matchCount > $alt ) {
				$alt = $matchCount;
				$altB = $base;
			}
		}
	}

	if ( $maxB ne '' ) {
		$seq .= $maxB;
		if ( $altB eq '' || $alt < $altCount || ($alt/$total) < $altFreq ) {
			$seqAlt .= $maxB;
		} else {
			$seqAlt .= $altB;
		}		
	} elsif ( $deleteByAmbig ) {
		$seq .= 'N';
		$seqAlt .= 'N';
	} elsif ( $keepDeleted ) {
		$seq .= $REF_SITES[$j];
		$seqAlt .= $REF_SITES[$j];	
	} # else deleted

	if ( $j == 0 ) {
		$max5 = $max;
	}
	
	if ( $j == ($Ncount-1) ) {
		$max3 = $max;
	}
}


if ( $notSkipExtension ) {
	# LEADER
	$leader =''; $threshold = max($minPadCount,int($max5/$denom)+1);
	foreach $p ( 1 .. scalar(keys(%count5)) ) {
		$max = 0; $maxB = ''; $total = 0;
		while( ($base, $leaderCount) = each(%{$count5{"-$p"}}) ) {
			if ($leaderCount > $max) {
				$max = $leaderCount;
				$maxB = $base;
			}
			$total += $leaderCount;
		}
		
		if ( $max < $threshold || $maxB !~ /[ACTGactg]/ ) {
			last;
		} else {
			$leader .= $maxB;
			$threshold = max($minPadCount,int($max/$denom)+1);
		}
	}
	$leader = reverse($leader);

	# TRAILER
	$trailer = ''; $threshold = max($minPadCount,int($max3/$denom)+1);
	foreach $p ( 0 .. scalar(keys(%count3)) - 1 ) {
		$max = 0; $maxB = ''; $total = 0;
		while( ($base, $trailerCount) = each(%{$count3{$p}}) ) {
			if ($trailerCount > $max) {
				$max = $trailerCount;
				$maxB = $base;
			}
			$total += $trailerCount;
		}
		
		if ( $max < $threshold || $maxB !~ /[ACTGactg]/ ) {
			last;
		} else {
			$trailer .= $maxB;
			$threshold = max($minPadCount,int($max/$denom)+1);
		}
	}
	print $header,$leader,$seq,$trailer,"\n";
	if ( $seqAlt ne '' && $seq ne $seqAlt ) {
		print $header2,$leader,$seqAlt,$trailer,"\n";
	}
} else {
	print $header,$seq,"\n";
	if ( $seqAlt ne '' && $seq ne $seqAlt ) {
		print $header2,$seqAlt,"\n";
	}
}
