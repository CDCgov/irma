#!/usr/bin/env perl
# Sam Shepard - fastQ_converter - 1.2014
use warnings;
use Getopt::Long;
Getopt::Long::Configure('no_ignore_case');
GetOptions(	'read-quality|T=i'=> \$qualityThreshold,
		'use-median|M' => \$useMedian,
		'fastQ-ouput|Q' => \$fastQformat,
		'base-quality|B=i' => \$baseQuality,
		'min-length|L=i' => \$minLength,
		'complement-add|C' => \$complementAndAdd,
		'ordinal-headers|O' => \$ordinal,
		'file-id|F=s' => \$fileID,
		'save-quality|S=s' => \$saveFile,
		'save-stats|A=s' => \$saveStats,	 
		'skip-remaining|K' => \$skipRemaining,
		'log-file|G=s' => \$logFile,
		'log-id|g=s' => \$logID,
		'keep-header|H' => \$keepHeader,
		'mask-adapter|m=s' => \$maskAdapter,
		'clip-adapter|c=s' => \$clipAdapter,
		'fuzzy-adapter|Z' => \$fuzzyAdapter,
		'uracil-to-thymine|U' => \$uracilToThymine,
		'enforce-clipped-length|E' => \$clippedMinLength
	);


if ( -t STDIN && scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <fastQ> [options]\n";
	$message .= "\t\t-T|--read-quality <threshold>\t\tSpecify the average read quality.\n";
	$message .= "\t\t-M|--use-median\t\t\t\tInterprets the threshold (-T) as the median, not average.\n";
	$message .= "\t\t-Q|--fastQ-output\t\t\tOutputs fastQ instead of fastA format.\n";
	$message .= "\t\t-L|--min-length <threshold>\t\tMinimum length of sequence data, default = 0.\n";
	$message .= "\t\t-C|--complement-add\t\t\tTake the reverse complement and add to data.\n";
	$message .= "\t\t-O|--ordinal-headers\t\t\tReplace header with strict ordinals.\n";
	$message .= "\t\t-F|--file-id <STR>\t\t\tFile id for ordinals.\n";
	$message .= "\t\t-S|--save-quality <STR>\t\t\tSave quality file for back-mapping.\n";
	$message .= "\t\t-A|--save-stats <STR>\t\t\tSave quality length file for analysis.\n";
	$message .= "\t\t-K|--skip-remaining\t\t\tDo not output data FASTA/FASTQ data (assumes -A).\n";
	$message .= "\t\t-H|--keep-header\t\t\tKeep header as usual.\n";
	$message .= "\t\t-c|--clip-adapter <STR>\t\t\tClip adapter.\n";
	$message .= "\t\t-m|--mask-adapter <STR>\t\t\tMask adapter.\n";
	$message .= "\t\t-Z|--fuzzy-adapter\t\t\tAllow one mismatch.\n";
	$message .= "\t\t-U|--uracil-to-thymine\t\t\tCovert uracil to thymine.\n";
	$message .= "\t\t-E|--enforce-clipped-length\t\tThe minimum length threshold (-L) is enforced when adapter clipped (-c).\n";
	die($message."\n");
}

$uracilToThymine = defined($uracilToThymine) ? 1 : 0;
$clippedMinLength = defined($clippedMinLength) ? 1 : 0;

if ( defined($clipAdapter) ) {
	$forwardAdapter = $clipAdapter;
	$reverseAdapter = reverse($forwardAdapter);
	$reverseAdapter =~ tr/ATGC/TACG/;
	$clipAdapter = 1;
} else {
	$clipAdapter = 0;
}

if ( defined($maskAdapter) ) {
	$forwardAdapter = $maskAdapter;
	$reverseAdapter = reverse($forwardAdapter);
	$reverseAdapter =~ tr/ATGC/TACG/;
	$adapterMask = 'N' x length($forwardAdapter);
	$maskAdapter = 1;
} else {
	$maskAdapter = 0;
}

if ( defined($keepHeader) ) {
	$keepHeader = 1;
	$notKeepHeader = 0;
} else {
	$keepHeader = 0;
	$notKeepHeader = 1;
}

if ( !defined($useMedian) ) {
	$useMedian = 0;
} else {
	$useMedian = 1;
}

if ( !defined($fileID) ) {
	$fileID = '';
} else {
	if ( $complementAndAdd ) {
		$fileID = '|'.$fileID.'|';
	} else {
		$fileID = $fileID.'|';
	}
}

@fuzzyAdaptersFWD = (); @fuzzyAdaptersREV = ();
if ( defined($fuzzyAdapter) && ($maskAdapter || $clipAdapter) ) {
	$L = length($forwardAdapter);
	if ( $L > 0 ) {
		$N = '[ATCGN]'; $L--;
		for $i (0 .. $L) {
			$tmp = $forwardAdapter;
			substr($tmp,$i,1,$N);
			push(@fuzzyAdaptersFWD, $tmp);
			$tmp = $reverseAdapter;
			substr($tmp,$i,1,$N);
			push(@fuzzyAdaptersREV, $tmp);
		}
	}
	$fuzzyAdapter = 1;
} else { $fuzzyAdapter = 0; }

if ( $saveFile ) {
	open(QUA,'>',$saveFile) or die("Cannot open $saveFile.\n");
	$saveQuality = 1;
} else {
	$saveQuality = 0;
}

if ( !defined($qualityThreshold) ) {
	$qualityThreshold = 0;
}

if ( defined($baseQuality) ) {
	$maskBases = 1;
	$baseQuality = int($baseQuality);
	die("Not implemented baseQuality.\n");
} else {
	$maskBases = 0;
}

if ( !defined($minLength) ) {
	$minLength = 0;
} elsif ( $minLength < 0 ) {
	die("ERROR: minimum length must be a non-negative integer.\n");
} else {
	$minLength = int($minLength);
}


$/ = "\n";
$act = $id = 0;

if ( $complementAndAdd ) {
	$pp = 'P';
	$nn = 'N';
	$dnp = '_';
} else {
	$pp = '';
	$nn = '';
	$dnp = '';
}

if ( defined($saveStats) ) {
	open(STAT,'>',$saveStats) or die("Cannot open $saveStats for writing.\n"); 
}

while($hdr=<>) {
	chomp($hdr);
	$seq=<>; chomp($seq);
	$junk=<>; chomp($junk);
	$quality=<>; chomp($quality);

	if ( length($seq) < $minLength ) {
		next;
	}

	if ( $uracilToThymine ) {
		$seq =~ tr/uU/tT/;
	}

	if ( $clipAdapter ) {
		if ( $seq  =~ /$reverseAdapter/i ) {
			$seq = substr($seq,0,$-[0]);
			$quality = substr($quality,0,$-[0]);
		} elsif ( $seq =~ /$forwardAdapter/i ) {
			$seq = substr($seq,$+[0]);	
			$quality = substr($quality,$+[0]);
		} elsif ( $fuzzyAdapter ) {
			my $findFuzzyFWD = 1;
			foreach $tmp ( @fuzzyAdaptersREV ) {
				if ( $seq =~ /$tmp/i ) {
					$seq = substr($seq,0,$-[0]);
					$quality = substr($quality,0,$-[0]);
					$findFuzzyFWD = 0;
					last;
				}
			}

			if ( $findFuzzyFWD ) {
				foreach $tmp ( @fuzzyAdaptersFWD ) {
					if ( $seq =~ /$tmp/i ) {
						$seq = substr($seq,$+[0]);	
						$quality = substr($quality,$+[0]);
						last;
					}
				}
			}
		}

		if ( $clippedMinLength &&  length($seq) < $minLength ) {
			next;
		}
	} elsif ( $maskAdapter ) {
		if ( $seq  =~ /$reverseAdapter/i ) {
			$seq =~ s/$reverseAdapter/$adapterMask/i;
		} elsif ( $seq =~ /$forwardAdapter/i ) {
			$seq =~ s/$forwardAdapter/$adapterMask/i;
		} elsif ( $fuzzyAdapter ) {
			my $findFuzzyFWD = 1;
			foreach $tmp ( @fuzzyAdaptersREV ) {
				if ( $seq =~ /$tmp/i ) {
					$seq =~ s/$tmp/$adapterMask/i;
					$findFuzzyFWD = 0;
					last;
				}
			}

			if ( $findFuzzyFWD ) {
				foreach $tmp ( @fuzzyAdaptersFWD ) {
					if ( $seq =~ /$tmp/i ) {
						$seq =~ s/$tmp/$adapterMask/i;
						last;
					}
				}
			}
		}
	}

	@a = unpack("c* i*",$quality);
	$q = 0; $n = scalar(@a);
	$id++;

	if ( $useMedian ) {
		@sorted = sort(@a);
		if ( $n % 2 == 0 ) {
			$q = ($sorted[$n/2] + $sorted[($n/2)-1]) / 2 - 33;
		} else {
			$q = $sorted[($n-1)/2] - 33;
		}
	
	# use average
	} else {
		foreach $x ( @a ) {
			$q += $x;
		}
		$q = ($q-$n*33)/$n;
	}

	if ( defined($saveStats) ) {
		print STAT $id,"\t",$q,"\t",$n,"\n";
		if ( defined($skipRemaining) ) {
			next;
		} 
	}

	if ( $q >= $qualityThreshold ) {
		$act++;
		$hdr = substr($hdr,1);

		if ( $notKeepHeader ) { $hdr =~ tr/ /_/; }

		if ( $ordinal ) {
			if ( $fastQformat ) {
				print $pp,$fileID,$id,"\n",$seq,"\n",$junk,"\n",$quality,"\n";
			} else {
				print '>',$pp,$fileID,$id,"\n",$seq,"\n";
				if ( $saveQuality ) {
					print QUA $pp,$fileID,$id,"\t",$hdr,"\t",$quality,"\n";
				}
			}
		} else {
			if ( $fastQformat ) {
				print '@',$hdr,$dnp,$pp,"\n",$seq,"\n",$junk,"\n",$quality,"\n";
			} else {
				print '>',$hdr,$dnp,$pp,'|',sprintf("%.1f",$q),'|',$n,"\n",$seq,"\n";
				if ( $saveQuality ) {
					print QUA $hdr,$dnp,$pp,"\t",$quality,"\n";
				}
			}
		}

		# Take the reverse complement and add it
		# courtesy http://reverse-complement.com/ambiguity.html
		if ( $complementAndAdd ) {
			$seq = reverse( $seq );
			$seq =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
			$quality=reverse($quality);
			if ( $ordinal ) {
				if ( $fastQformat ) {
					print $nn,$fileID,$id,"\n",$seq,"\n",$junk,"\n",$quality,"\n";
				} else {
					print '>',$nn,$fileID,$id,"\n",$seq,"\n";
					if ( $saveQuality ) {
						print QUA $nn,$fileID,$id,"\t",$hdr,"\t",$quality,"\n";
					}
				}
			} else {
				if ( $fastQformat ) {
					print $hdr,$dnp,$nn,"\n",$seq,"\n",$junk,"\n",$quality,"\n";
				} else {
					print '>',$hdr,$dnp,$nn,'|',sprintf("%.1f",$q),'|',$n,"\n",$seq,"\n";
					if ( $saveQuality ) {
						print QUA $hdr,$dnp,$nn,"\t",$quality,"\n";
					}
				}
			}
		}
	}
}


if ( defined($saveStats) ) {
	close(STAT);
}

if ( $logFile ) {
	if ( defined($logID) ) {
		$logID = ':'.$logID;
	} else {
		$logID = '';
	}
	open(OUT,'>>',$logFile) or die("Cannot open $logFile for appending.\n");
	print OUT $logFile,"$logID\t",($id),"\t",($act),"\t",$qualityThreshold,"\t",$minLength,"\t",$useMedian,"\n";
	close(OUT);
}

if ( $saveQuality ) {
	close(QUA);
}
