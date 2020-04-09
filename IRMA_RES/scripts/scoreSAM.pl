#!/usr/bin/env perl
# Sam Shepard -- scoreSAM
# 8.2015

use Getopt::Long;
my $score_field;
GetOptions( 'score-field|F=i' => \$score_field );

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <sam> [--score-field|-F <int>]\n";
	die($message."\n");
}

$/ = "\n";
my $score1 = 0;

if ( !defined($score_field) || $score_field < 0 ) {
	$score_field = 11;
} else {
	$score_field = int($score_field) - 1;
}

open(SAM,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($line = <SAM>) {
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	chomp($line);
	
	#($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,$AS) = split("\t",$line);
	my @sam_fields = split("\t",$line);
	if ( $score_field > $#sam_fields ) {
		die("$0: ERROR in $ARGV[0], score_field $score_field > $#sam_fields\n");
	}

	my @sub_fields = split(':',$sam_fields[$score_field]);
	$score1 += $sub_fields[2];
}
close(SAM);
print $score1;
