#!/usr/bin/env perl
# Sam Shepard - 2017-04-12

use File::Basename;
use Getopt::Long;
use Storable;
GetOptions(	'separate-ha-na-og|T' => \$triplet,
		'classify|C' => \$classify,
		'groups|G=s' => \$classificationGroups, 
		'include-chimera|I' => \$includeChimera,
		'align-to-ref|A' => \$alignSequences,
		'skip-elongation|S' => \$skipExtension,
		'prefix|P=s' => \$prefix,
		'min-match-length|L=i' => \$minMatchLength
		);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <blat.txt> <blat.fasta>\n";
	$message .= "\t\t-T|--separate-ha-na-og\t\tSeparates the data into three groups for influenza.\n";
	$message .= "\t\t-C|--classify\t\t\tUse BLAT scores to sort/classify the sequences (best match).\n";
	$message .= "\t\t-G|--groups <STR>\t\tGroup by pattern in string for classifications.\n";
	$message .= "\t\t-I|--include-chimera\t\tInclude chimeric data in match data instead of skipping it.\n";
	$message .= "\t\t-A|--align-to-ref\t\tAlign data to ref using BLAT matches.\n";
	$message .= "\t\t-S|--skip-elongation\t\tSkip the elongation of the reference to find 5prime and 3prime regions.\n";
	$message .= "\t\t-P|--prefix <STR>\t\tPrefix to store gene-wise stats.\n";
	$message .= "\t\t-L|--min-match-length <INT>\t\tMinimum match length after accounting for Chimera.\n";
	die($message."\n");
}

# Post chimera match length, gets returned to UNMATCHED reads
if ( defined($minMatchLength) && int($minMatchLength) > 0 ) {
	$minMatchLength = int($minMatchLength);
	$filterMatchLength = 1;
} else {
	$minMatchLength = 0;
	$filterMatchLength = 0;
}

$elongateReference = defined($skipExtension) ? 0 : 1;

if ( $triplet || defined($classificationGroups) ) {
	$GroupFN = sub { if ( $_[0] =~ /([HN]A)/ ) { return $1; } else { return 'OG'; } };
	if ( defined($classificationGroups) ) {
		$triplet = 1;
		@list = split(':',$classificationGroups);
		if ( scalar(@list) == 2 ) {
			$other = $list[$#list];
		} elsif ( scalar(@list) == 1 ) {
			$other = 'X';
		} else {
			die("Expected only one semi-colon per gene_group list. You entered: $classificationGroups\n");
		}
		@geneGroups = split(',',$list[0]);

		$GroupFN = sub {
			my $fnGene = $_[0];
			my $fnPat = '';
			my @fnPats = @geneGroups;
			my $fnOtherwise = $other;
			foreach $fnPat ( @fnPats ) {
				if ( $fnGene =~ /$fnPat/ ) { return $fnPat; }
			}
			
			return $fnOtherwise;
		};
	}
}


# Reverse complement
sub rc($) {
	my $sequence = reverse( $_[0] );
	$sequence =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
	return $sequence;
}

sub getLength($) {
	my @v = split("\t",$_[0]);
	my @bLengths = split(',',$v[18]);
	my $sum = 0;
	foreach my $l ( @bLengths ) {
		$sum += $l;
	}

	return $sum;
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
	if ( $leaderLen && $sequence =~ /^([-]{1,9})[ACTGN]/ ) {
		$gap = length($1);
		if ( $leaderLen >= $gap ) {
			substr($sequence,0,$gap) = uc(substr($leader,-$gap));
			$leader = substr($leader,0,$leaderLen-$gap);
		}
	}

	if ( $trailerLen > 0 && $sequence =~ /[ACTGN]([-]{1,9})$/ ) {
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
%stats = %maxScore = %supergroups = ();
while($line=<IN>) {
	chomp($line);
	@v = split("\t",$line);
	($match,$mismatch,$strand,$query,$target) = ($v[0],$v[1],$v[8],$v[9],$v[13]);
	$target =~ s/\{S\d+\}$//;
	$score = $match - $mismatch;
	if ( $score >= $maxScore{$query} && $score > 0 ) {
	       	$maxGene{$query}{$target}[0] = $strand;
		$maxGene{$query}{$target}[1] = $line;
		$maxGene{$query}{$target}[2] = $score;
		$maxScore{$query} = $score;
		if ( $triplet ) {
			$supergroup = $GroupFN->($target); 
			$maxGene{$query}{$target}[3] = $supergroup;
			$supergroups{$supergroup} = 1;
		 }
	}
	$total{$query}++;
	$stats{$query}{$target}{$strand}++;
}
close(IN);

foreach $q ( keys(%stats) ) {
	@targets = keys(%{$maxGene{$q}});
	foreach $t ( @targets ) {
		if ( $maxGene{$q}{$t}[2] < $maxScore{$q} ) {
			delete($maxGene{$q}{$t});
		}
	}

	if ( $total{$q} > 1 ) {
		@queryGenes = keys(%{$stats{$q}});
		if ( scalar(@queryGenes) == 1 ) {
			@geneStrands = keys(%{$stats{$q}{$queryGenes[0]}});
			if ( scalar(@geneStrands) == 1 ) {
				$Q{$q} = 'm'; # multihit
			} else {
				$Q{$q} = 'c'; # chimera
			}
		} else {
			# there be ties
			@maxGenes = keys(%{$maxGene{$q}});
			if ( scalar(@maxGenes) > 1 ) {
				%tmStrands = ();
				foreach $gene ( @maxGenes ) {
					@geneStrands = keys(%{$stats{$q}{$gene}});
					if( scalar(@geneStrands) == 1 ) {
						if ( $Q{$q} ne 'ac' ) {
							$tmStrands{$geneStrands[0]} = 1;
							$Q{$q} = 'tm'; # tied multiples
						}
					} else {
						# ambiguous chimeric is enough for one of the genes to doom it
						$Q{$q} = 'ac';  # ambiguous chimera
					}
				}

				if ( scalar(keys(%tmStrands)) > 1 ) {
					$Q{$q} = 'xc';	# Reversible complements
				}

			} else {
				$gene = $maxGenes[0];
				@geneStrands = keys(%{$stats{$q}{$gene}});
				if( scalar(@geneStrands) == 1 ) {
					if ( $stats{$q}{$gene}{$geneStrands[0]} == 1 ) {
						$Q{$q} = 'ar'; # ambiguous regular
					} else {
						$Q{$q} = 'am'; # ambiguous multihit
					}
				} else {
					$Q{$q} = 'ac'; # ambiguous chimera
				}
			}
		}
	} else {
		$Q{$q} = 'r'; # regular
	}
}

$/ = ">";
$name = basename($ARGV[0],'.blat');
$path = dirname($ARGV[0]);


open(CHIM,'>',$path.'/'.$name.'.chim') or die("Cannot open $path/$name.chim for writing.\n");
open(MATCH,'>',$path.'/'.$name.'.match') or die("Cannot open $path/$name.match for writing.\n");
open(NOMATCH,'>',$path.'/'.$name.'.nomatch') or die("Cannot open $path/$name.nomatch for writing.\n");

%supergroupHandles = ();
if ( $classify ) { 
	open(CLASS,'>',"$path/$name.class") or die("Cannot write $path/$name.class.\n");
} elsif ( $triplet ) {
	foreach $sg ( keys(%supergroups) ) {
		open($supergroupHandle{$sg},'>',"$path/$name.match.$sg") or die("Cannot open $path/$name.match.$sg for writing.\n");
	}
}

$sortOutChim = defined($includeChimera) ? 0 : 1;
if ( $alignSequences ) { 
	%alignStats = ();
}

open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$q = $header = shift(@lines);
	$sequence = lc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	}

	%written = ();
	if ( defined($Q{$header})) {
		# chimera
		if ( $sortOutChim && index($Q{$header},'c') != -1 ) {
			print CHIM '>',$header,"\n",$sequence,"\n";
		} else {
			foreach $gene ( keys(%{$maxGene{$q}}) ) {
				($strand,$matchLine,$score,$sg) = (@{$maxGene{$q}{$gene}});
				$seq2 = $strand eq '+' ? $sequence : rc($sequence);
				$tag = $strand eq '+' ? '' : '{c}';

				if ( $filterMatchLength && getLength($matchLine) < $minMatchLength ) {
					print NOMATCH '>',$header,"\n",$sequence,"\n";
					next;	
				}

				if ( defined($alignSequences) ) { recordStats(\%alignStats,$gene, alignedBLAT($matchLine,$seq2,$gene) ); }
				if ( $classify ) { print CLASS $header,$tag,"\t",$gene,"\t",$score,"\n"; }
				# write MATCH
				if ( !defined($written{$strand}) ) { 
					print MATCH '>',$header,$tag,"\n",$seq2,"\n";
					$written{$strand} = 1;
				}

				# write SUPERGROUP
				if ( $triplet && !defined($written{$strand.$sg}) ) {
					$handle = $supergroupHandle{$sg};
					print $handle '>',$header,$tag,"\n",$seq2,"\n";
					$written{$strand.$sg} = 1;
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
	foreach $sg ( keys(%supergroups) ) {
		close($supergroupHandle{$sg});
	}
}

$prefix = !defined($prefix) ? '' : $prefix.'-';
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
