#!/usr/bin/env perl
# Sam Shepard -- Winnow SAM 
# 8.2015

use Getopt::Long;
GetOptions(
		'use-matches|M' => \$useMatches, 'in-place|I' => \$inPlace, 'interleaved-pairs|P' => \$interleavedPairs,
		'score-field|F=i' => \$scoreField
	);

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <sam>\n";
	die($message."\n");
}

if ( !defined($scoreField) || $scoreField < 0 ) {
	$scoreField = 0;
} else {
	if ( $scoreField >= 12 ) {
		$scoreField -= 12;
	}
}


# FUNCTIONS #
sub countMatch($) {
	my $cig = $_[0];
	my $count = 0;
	while( $cig =~ /(\d+)([MIDNSHP])/g ) {
		if ( $2 eq 'M' ) {
			$count += $1;
		}
	}
	return $count;
}

#############
$pair = 0;
$/ = "\n"; $previous = $header = '';
open(SAM,'<',$ARGV[0]) or die("$0 ERROR: cannot open $ARGV[0] for reading.\n");
while($line = <SAM>) {
	if ( substr($line,0,1) eq '@' ) {
		if ( $line ne $previous ) {
			$header .= $line;
		}
		$previous = $line;
		next;
	}

	
	chomp($line);
	my ($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,@extra_fields) = split("\t",$line);
	my $AS = $extra_fields[$scoreField];	

	if ( $interleavedPairs ) {
		$qname = $qname . '_' . ($pair%2);
		$pair++;
	}

	if ( $cigar eq '*' ) {
		next;
	}

	if ( $useMatches ) {
		$score = countMatch($cigar);
	} elsif ( $AS =~ /AS:\w:(\d+)/ ) {
		$score = $1;
	} else {
		print STDERR "Warning, using simple matches as back-up: $qname\n";
		$score = countMatch($cigar);
	}

	if ( !defined($scoreByQuery{$qname}) || $scoreByQuery{$qname} < $score ) {
		$scoreByQuery{$qname} = $score;
		$recordByQuery{$qname} = $line;
	}
}
close(SAM);

if ( $inPlace ) {
	open(SAM,'>',$ARGV[0]) or die("$0 ERROR: cannot open $ARGV[0] for writing.\n");
	print SAM $header;
	foreach $query ( keys(%recordByQuery) ) {
		print SAM $recordByQuery{$query},"\n";
	}
	close(SAM);
} else {
	print $header;
	foreach $query ( keys(%recordByQuery) ) {
		print $recordByQuery{$query},"\n";
	}
}
