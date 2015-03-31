#!/usr/bin/env perl
# Sam Shepard - fastQ_converter - 1.2014

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
		'log-file|G=s' => \$logFile,
		'log-id|g=s' => \$logID
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
	die($message."\n");
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

while($hdr=<>) {
	chomp($hdr);
	$seq=<>; chomp($seq);
	$junk=<>; chomp($junk);
	$quality=<>; chomp($quality);
	@a = unpack("c* i*",$quality);

	$q = 0; $n = scalar(@a);

	$id++;
	if ( length($seq) < $minLength ) {
		next;
	}

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

	if ( $q >= $qualityThreshold ) {
		$act++;
		$hdr = substr($hdr,1);
		$hdr =~ s/ /_/;

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
