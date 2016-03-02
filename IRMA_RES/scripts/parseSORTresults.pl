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

if ( !defined($minimumRcount) || $minimumRcount < 1 ) {
	$minimumRcount = 1;
}

if ( !defined($minimumRPcount) || $minimumRPcount < 1 ) {
	$minimumRPcount = 1;
}

%counts = ();
$/ = "\n";
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($line=<IN>) {
	chomp($line);
	($ID,$target,$sig_threshold) = split("\t",$line);
	if ( $ignoreAnnotations && $target =~ /^([^{]+){[^}]*}$/ ) {
		$target = $1;
	}

	$counts{$target}++;
	$IDs{$ID} = $target;
	if ( $ID =~ /C\d+%(\d+)%/ ) {
		$rCounts{$target} += $1;
	}
}
close(IN);


foreach $target ( %counts ) {
	if ( !defined( $rCounts{$target} ) ) {
		$rCounts{$target} = $counts{$target};
	}
}

# choose primary or secondary between groups
open(OUT,'>',$ARGV[2].'.txt') or die("Cannot open $ARGV[0].txt\n");
@genes = sort { $counts{$b} <=> $counts{$a} } keys(%counts);
foreach $gene ( @genes ) {
	print OUT $gene,"\t",$counts{$gene},"\t",$rCounts{$gene},"\n";
	if  ( $counts{$gene} >= $minimumRPcount && $rCounts{$gene} >= $minimumRcount ) {
		$valid{$gene} = 1;
	} else {
		$valid{$gene} = 0;
	}
}
close(OUT);

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

if ( defined($patternList) && length($patternList) > 0) {
	@patterns = split(',',$patternList);
	if ( scalar(@patterns) == 1 && $patterns[0] eq '__ALL__' ) {
		for($i=1;$i<scalar(@genes);$i++) {
			$valid{$genes[$i]} = 0;
		}
		@patterns = ();
	}

	%genesByPat = ();
	foreach $pat ( @patterns ) {
		foreach $gene ( @genes ) {
			# can enforce FCFS uniqueness for each pat
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
	
%handles = ();
foreach $gene ( @genes ) {
	if ( $valid{$gene} > 0 ) {
		$file = $ARGV[2].'-'.$gene.'.fa';
	} else {
		$file = $ARGV[2].'-'.$gene.'.fa.2';
	}
	open($handles{$gene},'>',$file) or die("Cannot open $file for writing.\n");
}

open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n"); $/ = '>';
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));

	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	}

	$gene = $IDs{$header};
	$handle = $handles{$gene};
	print $handle '>',$header,"\n",$sequence,"\n";
}
close(IN);

foreach $gene ( keys(%handles) ) {
	$fh = $handles{$gene};
	close($fh);
}
