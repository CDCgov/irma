#!/usr/bin/env perl
# Sam Shepard - 3.18.2014

use File::Basename;
use Getopt::Long;
use Storable;
GetOptions(	'separate-ha-na-og|T' => \$triplet,
		'classify|C' => \$classify,
		'groups|G=s' => \$classificationGroups, 
		'include-chimera|I' => \$includeChimera,
		'align-to-ref|A' => \$alignSequences,
		'skip-elongation|S' => \$skipExtension,
		'prefix|P=s' => \$prefix
		);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <blat.txt> <blat.fasta>\n";
	$message .= "\t\t-T|--separate-ha-na-og\t\tSeparates the data into three groups for influenza.\n";
	$message .= "\t\t-C|--clasify\t\t\tUse BLAT scores to sort/classify the sequences (best match).\n";
	$message .= "\t\t-G|--groups <STR>\t\tGroup by pattern in string for classifications.\n";
	$message .= "\t\t-I|--include-chimera\t\tInclude chimeric data in match data instead of skipping it.\n";
	$message .= "\t\t-A|--align-to-ref\t\tAlign data to ref using BLAT matches.\n";
	$message .= "\t\t-S|--skip-elongation\t\tSkip the elongation of the reference to find 5prime and 3prime regions.\n";
	$message .= "\t\t-P|--prefix <STR>\t\tPrefix to store gene-wise stats.\n";
	die($message."\n");
}


if ( defined($skipExtension) ) {
	$elongateReference = 0;
} else {
	$elongateReference = 1;
}

sub alignedBLAT($$$) {
	my @v = split("\t",$_[0]);
	my $s = $_[1];
	my $h = $_[2];
	my $qSize = $v[10];
	my $tSize =  $v[14];
	my $bCount = $v[17];
	my @bLengths = split(',',$v[18]);
	my @qStarts = split(',',$v[19]);
	my @tStarts = split(',',$v[20]);
	my $left = 0;
	my $right = 0;
	my $seq = '';
	my $leader = '';
	my $trailer = '';
	
	if ( scalar(@v) < 22 ) {
		if ( $bCount == 1 ) {
			$left = $tStarts[0];
			$right = $tSize-($tStarts[0]+$bLengths[0]);

			if ( $qStarts[0] > 0 ) {
				$leader = substr($s,0,$qStarts[0]);	
			}
			$seq .=('-'x$left).uc(substr($s,$qStarts[0],$bLengths[0])).('-'x$right);
			$right = $qStarts[0]+$bLengths[0];
			if ( $right < $qSize  ) {
				$trailer = substr($s,$right,$qSize-$right);
			}
		} else {
			if ( $qStarts[0] > 0 ) {
				$leader = substr($s,0,$qStarts[0]);	
			}

			$left = $tStarts[0];
			$seq .=('-'x$left).uc(substr($s,$qStarts[0],$bLengths[0]));
			for(my $i = 1; $i < $bCount; $i++) {
				$left = $tStarts[$i] - $tStarts[$i-1] - $bLengths[$i-1];
				$seq .= ('-'x$left).uc(substr($s,$qStarts[$i],$bLengths[$i]));
				
			}

			$right = $tSize-($tStarts[$bCount-1]+$bLengths[$bCount-1]);
			$seq .= '-'x$right;
			
			$right = $qStarts[$bCount-1]+$bLengths[$bCount-1];
		
			if ( $right < $qSize  ) {
				$trailer = substr($s,$right,$qSize-$right);
			}
		}
	# PSLX mode
	} else {
		my @qSeqs = split(',',$v[21]);
		if ( $bCount == 1 ) {
			$left = $tStarts[0];
			$right = $tSize-($tStarts[0]+$bLengths[0]);
			if ( $qStarts[0] > 0 ) {
				$leader = substr($s,0,$qStarts[0]);	
			}
			$seq .=('-'x$left).uc($qSeqs[0]).('-'x$right);
			$right = $qStarts[0]+$bLengths[0];
			if ( $right < $qSize  ) {
				$trailer = substr($s,$right,$qSize-$right);
			}
		} else {
			if ( $qStarts[0] > 0 ) {
				$leader = substr($s,0,$qStarts[0]);	
			}

			$left = $tStarts[0];
			$seq .=('-'x$left).uc($qSeqs[0]);
			for(my $i = 1; $i < $bCount; $i++) {
				$left = $tStarts[$i] - $tStarts[$i-1] - $bLengths[$i-1];
				$seq .= ('-'x$left).uc($qSeqs[$i]);
			}

			$right = $tSize-($tStarts[$bCount-1]+$bLengths[$bCount-1]);
			$seq .= '-'x$right;
			$right = $qStarts[$bCount-1]+$bLengths[$bCount-1];
			if ( $right < $qSize  ) {
				$trailer = substr($s,$right,$qSize-$right);
			}
		}
	}

#	if ( $h =~ /MP/ ) {print STDERR '>',$v[9],"\n",$seq,"\n"; }
	return ($leader,$seq,$trailer);
}

sub recordStats($$@) {
	my $hashRef = $_[0];
	my $gene = $_[1];
	my $leader = $_[2];
	my $sequence = $_[3];
	my $trailer = $_[4];
	my $x = '';
	my $gap = '';
	my $length = 0;
	my $leaderLen = 0;
	my $trailerLen = 0;

	$length = length($sequence);

	# repair as necessary
	$trailerLen = length($trailer);
	$leaderLen = length($leader);
	if ( $leaderLen && $sequence =~ /^([-]{1,7})[ACTGN]/ ) {
		$gap = length($1);
		if ( $leaderLen >= $gap ) {
			substr($sequence,0,$gap) = uc(substr($leader,-$gap));
			$leader = substr($leader,0,$leaderLen-$gap);
		}
	}

	if ( $trailerLen > 0 && $sequence =~ /[ACTGN]([-]{1,7})$/ ) {
		$gap = length($1);
		if ( $trailerLen >= $gap ) {
			substr($sequence,$-[1],$gap) = uc(substr($trailer,0,$gap));
			$trailer = substr($trailer,-($trailerLen-$gap));
		}
	}

	for $x ( 0 .. ($length-1) ) {
		 $hashRef->{$gene}[0][$x]{substr($sequence,$x,1)}++;
	}

	if ( $elongateReference ) {
		if ( $sequence =~ /^[ACTGN]/ ) {
			for $x ( - length($leader) .. -1 ) {
				$hashRef->{$gene}[1]{$x}{substr($leader,$x,1)}++;
			}
		}

		if ( $sequence =~ /[ATCGN]$/ ) {
			for $x ( 0..length($trailer)-1 ) {
				$hashRef->{$gene}[2]{$x}{substr($trailer,$x,1)}++;
			}
		}
	}
}

$/ = "\n";
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
for(1..5) { $junk=<IN>; }
$prevID = 'ID'; $prevStr = 'strand'; 
%stats = ();
while($line=<IN>) {
	chomp($line);
	@v = split("\t",$line);
	($match,$mismatch,$strand,$query,$target) = ($v[0],$v[1],$v[8],$v[9],$v[13]);
	$score = $match - $mismatch;
	if ( $score > $maxGene{$query}[1] ) {
		$maxGene{$query}[1] = $score;
	       	$maxGene{$query}[0] = $target;
		$maxGene{$query}[2] = $line;
		if ( $triplet ) {
			if ( $target =~ /_HA/ ) {
				$class = 'HA';
			} elsif ( $target =~ /_NA/ ) {
				$class = 'NA';
			} else {
				$class = 'OG';
			}

			$C{$query} = $class;
		}
	}
	$total{$query}++;
	$stats{$query}{$target}{$strand}++;
}
close(IN);

foreach $q ( keys(%stats) ) {
	#print $q;
	if ( $total{$q} > 1 ) {
		@queryGenes = keys(%{$stats{$q}});
		if ( scalar(@queryGenes) == 1 ) {
			@geneStrands = keys(%{$stats{$q}{$queryGenes[0]}});
			if ( scalar(@geneStrands) == 1 ) {
				#print "\t",'multihit';
				$Q{$q} = 'm';
				$S{$q} = $geneStrands[0];
			} else {
				$Q{$q} = 'c';
				#print "\t",'chimeric';
			}
		} else {
			$gene = $maxGene{$q}[0];
			@geneStrands = keys(%{$stats{$q}{$gene}});
			if( scalar(@geneStrands) == 1 ) {
				if ( $stats{$q}{$gene}{$geneStrands[0]} == 1 ) {
					#print "\t",'ambiguous regular';
					$Q{$q} = 'ar';
				} else {
					#print "\t",'ambiguous multihit';
					$Q{$q} = 'am';
				}
				$S{$q} = $geneStrands[0];
			} else {
				#print "\t",'ambiguous chimeric';
				$Q{$q} = 'ac';
			}
		}
	} else {
		#print "\t",'regular';
		$Q{$q} = 'r';
		$gene = $maxGene{$q}[0];
		@geneStrands = keys(%{$stats{$q}{$gene}});
		$S{$q} = $geneStrands[0];
	}
	#print "\n";
}

$/ = ">";
$name = basename($ARGV[0],'.blat');
$path = dirname($ARGV[0]);

open(CHIM,'>',$path.'/'.$name.'.chim') or die("Cannot open $path/$name.chim for writing.\n");
open(MATCH,'>',$path.'/'.$name.'.match') or die("Cannot open $path/$name.match for writing.\n");
open(NOMATCH,'>',$path.'/'.$name.'.nomatch') or die("Cannot open $path/$name.nomatch for writing.\n");

if ( $classify ) { 
	open(CLASS,'>',"$path/$name.class") or die("Cannot write $path/$name.class.\n");
} elsif ( $triplet ) {
	open(HA,'>',$path.'/'.$name.'.match.HA') or die("Cannot open $path/$name.match.HA for writing.\n");
	open(OG,'>',$path.'/'.$name.'.match.OG') or die("Cannot open $path/$name.match.OG for writing.\n");
	open(NA,'>',$path.'/'.$name.'.match.NA') or die("Cannot open $path/$name.match.NA for writing.\n");
}

if ( $alignSequences ) { 
	%alignStats = ();
}

if ( defined($includeChimera) ) {
	$sortOutChim = 0;
} else {
	$sortOutChim = 1;
}

open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	}

	if ( defined($Q{$header})) {
		if ( $sortOutChim && $Q{$header} =~ /c/o ) {
			print CHIM '>',$header,"\n",$sequence,"\n";
		} else {
			if ( $S{$header} eq '+' ) {
				print MATCH '>',$header,"\n",$sequence,"\n";
				if ( defined($alignSequences) ) {
					recordStats(\%alignStats,$maxGene{$header}[0], alignedBLAT($maxGene{$header}[2],$sequence,$maxGene{$header}[0]) );
				}

				if ( $classify ) {
					print CLASS $header,"\t",$maxGene{$header}[0],"\t",$maxGene{$header}[1],"\n";
				}
				if ( $triplet ) {
					if ( $C{$header} eq 'HA' ) {
						print HA '>',$header,"\n",$sequence,"\n";
					} elsif ( $C{$header} eq 'NA' ) {
						print NA '>',$header,"\n",$sequence,"\n";
					} else {
						print OG '>',$header,"\n",$sequence,"\n";
					}
				}
			} else {
				$sequence = reverse( $sequence );
				$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;

				print MATCH '>',$header,"{c}\n",$sequence,"\n";
				if ( defined($alignSequences) ) {
					recordStats(\%alignStats, $maxGene{$header}[0],alignedBLAT($maxGene{$header}[2],$sequence,$maxGene{$header}[0]));
				}
 
				if ( $classify ) {
					print CLASS $header,"{c}","\t",$maxGene{$header}[0],"\t",$maxGene{$header}[1],"\n";
				}
				if ( $triplet ) {
					if ( $C{$header} eq 'HA' ) {
						print HA '>',$header,"{c}\n",$sequence,"\n";
					} elsif ( $C{$header} eq 'NA' ) {
						print NA '>',$header,"{c}\n",$sequence,"\n";
					} else {
						print OG '>',$header,"{c}\n",$sequence,"\n";
					}
				}
			}
		}
	} else {
		print NOMATCH '>',$header,"\n",$sequence,"\n";
	}
}
close(IN);
close(NOMATCH);
close(MATCH);
close(CHIM);
if ( $classify ) { 
	close(CLASS);
} elsif ( $triplet ) {
	close(HA); close(NA); close(OG);
}

if ( !defined($prefix) ) {
	$prefix = '';
} else {
	$prefix = $prefix . '-';
}

if ( $name =~ /_(\d{4,})$/ ) {
	$suffix = '_'.$1;
} else {
	$suffix = '';
}

if ( defined($alignSequences) ) {
	foreach $gene ( keys(%alignStats) ) {
		$filename = $path.'/'.$prefix.$gene.$suffix.'.sto';
		@count = ();
		@count = @{$alignStats{$gene}};
		store(\@count, $filename);
	}
} 

