#!/usr/bin/env perl
# parseSORTresults.pl
# Sam Shepard 9.2014

use Getopt::Long;
GetOptions(	'pattern-list|P=s' => \$patternList, 'ignore-annotations|G' => \$ignoreAnnotations,
		'min-read-count|C=i' => \$minimumRcount, 'min-read-patterns|D=i' => \$minimumRPcount,
		'ban-list|B=s' => \$banList );

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <SORT_results.tab> <match.FASTA> <PREFIX> [options]";
	die($message."\n");
}

if ( !defined($minimumRcount)  || $minimumRcount < 1 )  { $minimumRcount  = 1; }
if ( !defined($minimumRPcount) || $minimumRPcount < 1 ) { $minimumRPcount = 1; }

my %counts = (); 
my %IDs = ();
my %rCounts = ();
$/ = "\n";
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($line=<IN>) {
	chomp($line);
	if ( length($line) == 0 ) { next; }
	($ID,$target,$score) = split("\t",$line);
	if ( $ignoreAnnotations && $target =~ /^([^{]+)\{[^}]*}$/ ) { $target = $1; }

	$counts{$target}++;
	$IDs{$ID}{$target} = $score;
	if ( $ID =~ /C\d+%(\d+)%/ ) { $rCounts{$target} += $1; }
}
close(IN);


foreach $target ( %counts ) {
	if ( !defined( $rCounts{$target} ) ) {
		$rCounts{$target} = $counts{$target};
	}
}

# choose primary or secondary between groups
open(OUT,'>',$ARGV[2].'.txt') or die("Cannot open $ARGV[2].txt\n");
@genes = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
foreach $gene ( @genes ) {
	print OUT $gene,"\t",$counts{$gene},"\t",$rCounts{$gene},"\n";
	# By default, anything with the valid minimum is okay
	if  ( $counts{$gene} >= $minimumRPcount && $rCounts{$gene} >= $minimumRcount ) {
		$valid{$gene} = 1;
	} else {
		$valid{$gene} = 0;
	}
}
close(OUT);

# We remove genes that are explicitly banned
if ( defined($banList) && length($banList) > 0 ) {
	@banPatterns = split(',',$banList);
	foreach $banPat ( @banPatterns ) {
		foreach $gene ( @genes ) {
			if ( $gene =~ /$banPat/ ) {
				$valid{$gene} = 0;
			}
		}
	}
}

# If we have a pattern list we must divide into groups
if ( defined($patternList) && length($patternList) > 0) {
	$patternList =~ tr/ //d;
	$patternList =~ tr/;/,/;
	$patternList =~ tr/\//|/;
	@patterns = split(',',$patternList);

	# If ALL genes are in a group, get the best (2nd best & on are invalid)
	if ( scalar(@patterns) == 1 && $patterns[0] eq '__ALL__' ) {
		for($i=1;$i<scalar(@genes);$i++) {
			$valid{$genes[$i]} = 0;
		}
		@patterns = ();
	} else {
		# Divide into group maximums
		%genesByPat = ();
		foreach $pat ( @patterns ) {
			foreach $gene ( @genes ) {
				# Can enforce FCFS uniqueness for each pattern
				if ( $gene =~ /$pat/ ) {
					$genesByPat{$pat}{$gene} = $counts{$gene};
				}
			}
			@geneList = sort { $genesByPat{$pat}{$b} <=> $genesByPat{$pat}{$a} } keys(%{$genesByPat{$pat}});
			for($i=1;$i<scalar(@geneList);$i++) {
				$valid{$geneList[$i]} = 0;
			}
		}
	}
}

# Set up handles for primary and secondary data	
my %handles = ();
foreach $gene ( @genes ) {
	# primary
	if ( $valid{$gene} > 0 ) {
		$file = $ARGV[2].'-'.$gene.'.fa';
	# secondary
	} else {
		$file = $ARGV[2].'-'.$gene.'.fa.2';
	}
	open($handles{$gene},'>',$file) or die("Cannot open $file for writing.\n");
}

if ( scalar(@genes) == 0 ) {
	die("$0 ERROR: no classification output found! Aborting.\n");
}

open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n"); $/ = '>';
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 ) { next; }

	%secondary = (); $primaryWritten = 0;
	foreach $gene ( keys(%{$IDs{$header}}) ) {
		if ( $valid{$gene} > 0 ) {
			$primaryWritten = 1;
			$handle = $handles{$gene};
			print $handle '>',$header,"\n",$sequence,"\n";
		} else {
			$secondary{$gene} = $IDs{$header}{$gene};
		}
	}

	if ( $primaryWritten == 0 ) {
		@genes = sort { $secondary{$b} <=> $secondary{$a} } keys(%secondary);
		$handle = $handles{$genes[0]};
		print $handle '>',$header,"\n",$sequence,"\n";
	}
}
close(IN);

foreach $gene ( keys(%handles) ) {
	$fh = $handles{$gene};
	close($fh);
}
