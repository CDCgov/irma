#!/usr/bin/env perl
# Sam Shepard - xflate - 1.2014

use Getopt::Long;
GetOptions(	'inflate|I' => \$inflate,
		'cluster-all|C' => \$clusterAll,
		'file-label|L=s' => \$label,
		'reverse-complement-inflate|R' => \$rci,
		'fastQ|Q' => \$fastQ
	);

if ( scalar(@ARGV) < 2 ) {
	$message = "Usage:\n\tperl $0 <xflate.txt> [options] <fasta1.fa> ...\n";
	$message .= "\t\t-I|--inflate\t\tInflate sequence files.\n";
	$message .= "\t\t-C|--cluster-all\tUse a cluster for all sequences.\n";
	$message .= "\t\t-L|--file-label <STR>\tLabels the clusters with an identifier.\n";
	$message .= "\t\t-R|--rev-comp-inf\tReverse complement inflation.\n";
	$message .= "\t\t-Q|--fastQ\t\tExpects fastQ input (deflate) and fastQ output (inflate).\n";
	die($message."\n");
}

if ( !defined($label) ) {
	$label = '';
} else {
	$label = '|'.$label;
}

$tableFile = shift(@ARGV);
# INFLATE
if ( $inflate ) {
	open(TABLE,'<',$tableFile) or die("Cannot open $tableFile for reading (trying to inflate without a file?).\n");
	$/ = "\n"; %headersByCluster = ();

	if ( $fastQ ) {
		while($line=<TABLE>) {
			chomp($line);
			@fields = split("\t",$line);
#			@{$headersByCluster{$fields[0]}} = ();
			for($i=1;$i<scalar(@fields);$i+=2) {
#				push(@{$headersByCluster{$fields[0]}},$fields[$i]);
				$headersByCluster{$fields[0]}[int($i/2)] = $fields[$i];
				$qualityByHeader{$fields[$i]} = $fields[$i+1];
			}
		}
	} else {
		while($line=<TABLE>) {
			chomp($line);
			@fields = split("\t",$line);
			@{$headersByCluster{$fields[0]}} = @fields[1..$#fields];
		}
	}
	close(TABLE);

	foreach $file ( @ARGV ) {
		$/ = '>';
		open( IN, '<', $file ) or die("Cannot open $file.\n");
		while( $record = <IN> ) {
			chomp($record);
			@lines = split(/\r\n|\n|\r/, $record);
			$header = shift(@lines);
			$sequence = lc(join('',@lines));

			$length = length($sequence);
			if ( $length == 0 ) {
				next;
			}
			
			if ( $header =~ /^(C\d+%\d+%[^{]*)/ ) {
				$cluster = $1;
				if ( $rci && $header =~ /{c}/ ) {
					$sequence = reverse($sequence);
					$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
				}
				
				if ( $fastQ ) {
					if ( defined($headersByCluster{$cluster}) ) {
						foreach $hdr ( @{$headersByCluster{$cluster}} ) {
							print STDOUT $hdr,"\n",$sequence,"\n+\n",$qualityByHeader{$hdr},"\n";
						}
					} else {
						$quality = '?' x length($sequence);
						print STDOUT '@',$header,"\n",$sequence,"\n+\n",$quality,"\n";
					}
				} else {
					if ( defined($headersByCluster{$cluster}) ) {
						foreach $hdr ( @{$headersByCluster{$cluster}} ) {
							print STDOUT '>',$hdr,"\n",$sequence,"\n";
						}
					} else {
						print STDOUT '>',$header,"\n",$sequence,"\n";
					}
				}
			} else {
				if ( $fastQ ) {
					$quality = '?' x length($sequence);
					print STDOUT '@',$header,"\n",$sequence,"\n+\n",$quality,"\n";
				} else {
					print STDOUT '>',$header,"\n",$sequence,"\n";
				}
			}
		}
	}
	
# DEFLATE
} else {
	if ( -e $tableFile ) {
		print STDERR "WARNING, using ${tableFile}2 since file exists.\n";
		$tableFile .= '2';
	}
	open(TABLE,'>',$tableFile) or die("Cannot open $tableFile for reading (trying to inflate without a file?).\n");
	%clustersBySequence = ();
	foreach $file ( @ARGV ) {
		open( IN, '<', $file ) or die("Cannot open $file.\n");
		if ( $fastQ ) {
			$/ = "\n";
			while($header=<IN>) {
				chomp($header);
				$sequence=<IN>; chomp($sequence);
				$junk=<IN>; chomp($junk);
				$quality=<IN>; chomp($quality);

				if ( !defined($clustersBySequence{$sequence}) ) {
					$clustersBySequence{$sequence}[0] = $header;
				} else {
					push(@{$clustersBySequence{$sequence}},$header);
				}

				$qualityByHeader{$header} = $quality;
			}
		} else {
			$/ = '>';
			while( $record = <IN> ) {
				chomp($record);
				@lines = split(/\r\n|\n|\r/, $record);
				$header = shift(@lines);
				$sequence = lc(join('',@lines));

				$length = length($sequence);
				if ( $length == 0 ) {
					next;
				}
				
				if ( !defined($clustersBySequence{$sequence}) ) {
					$clustersBySequence{$sequence}[0] = $header;
				} else {
					push(@{$clustersBySequence{$sequence}},$header);
				}
			}
		}
		close(IN);
	}

	$i = 0;
	foreach $seq ( keys(%clustersBySequence) ) {
		$clusterSize=scalar(@{$clustersBySequence{$seq}});
		if ( $clusterAll || $clusterSize > 1 || $clustersBySequence{$seq}[0] =~ /^C\d+?%/ ) {
			print TABLE "C${i}%",$clusterSize,'%',$label;
			if ( $clustersBySequence{$seq}[0] =~ /\{(.+?)\}/ ) {
				$annot = $1;	
				print STDOUT '>C',$i,'%',$clusterSize,"%$label\{$annot\}\n",$seq,"\n";
			} else {
				print STDOUT '>C',$i,'%',$clusterSize,"%$label\n",$seq,"\n";
			}

			if ( $fastQ ) {
				foreach $hdr ( @{$clustersBySequence{$seq}} ) {
					print TABLE "\t",$hdr,"\t",$qualityByHeader{$hdr};
				}
			} else {
				foreach $hdr ( @{$clustersBySequence{$seq}} ) {
					print TABLE "\t",$hdr;
				}
			}
			print TABLE "\n";
			$i++;
		} else {
			print STDOUT '>',$clustersBySequence{$seq}[0],"\n",$seq,"\n";
		}
	}
	close(TABLE);
}
