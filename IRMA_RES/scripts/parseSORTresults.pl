#!/usr/bin/env perl
# parseSORTresults.pl
# Sam Shepard 9.2014

use Getopt::Long;
GetOptions(	'pattern-list|P=s' => \$patternList, 'min-read-count|C=i' => \$minimumRcount, 'min-read-patterns|D=i' => \$minimumRPcount );

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <SORT_results.tab> <match.FASTA> <PREFIX> [significance]";
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
	($ID,$annot,$sig_threshold) = split("\t",$line);
	$counts{$annot}++;
	$IDs{$ID} = $annot;
	if ( $ID =~ /C\d+%(\d+)%/ ) {
		$rCounts{$annot} += $1;
	}
}
close(IN);

if ( defined($patternList) ) {
	@patterns = split(',',$patternList);
} else {
	@patterns = ('HA','NA','PA','PB1','PB2','NS','MP','NP');
}

open(OUT,'>',$ARGV[2].'.txt') or die("Cannot open $ARGV[0].txt\n");
@sorted = sort { $counts{$a} <=> $counts{$b} } keys(%counts);
%best = %bCount = ();
foreach $gene ( @sorted ) {
	print OUT $gene,"\t",$counts{$gene},"\t",$rCounts{$gene},"\n";
	foreach $pat ( @patterns ) {
		if ( $gene =~ /$pat/ && $counts{$gene} > $bCount{$pat} ) {
			$best{$pat} = $gene;
			$bCount{$pat} = $counts{$gene};
		}
	}
}
close(OUT);

foreach $pat ( @patterns ) {
	if ( $bCount{$pat} > 0 ) {
		$valid{$best{$pat}} = 1;
	}
}

%handles = ();
foreach $gene ( @sorted ) {
	if ( $valid{$gene} > 0 && $counts{$gene} >= $minimumRPcount && $rCounts{$gene} >= $minimumRcount ) {
		$file = $ARGV[2].'-'.$gene.'.fa';
	} else {
		$file = $ARGV[2].'-'.$gene.'.fa.2';
	}
	open($handles{$gene},'>',$file) or die("Cannot open $file for writing.\n");
}

open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
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

	$gene = $IDs{$header};
	$handle = $handles{$gene};
	print $handle '>',$header,"\n",$sequence,"\n";
}
close(IN);
