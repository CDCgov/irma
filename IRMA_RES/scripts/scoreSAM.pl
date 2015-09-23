#!/usr/bin/env perl
# Sam Shepard -- scoreSAM
# 8.2015

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <sam>\n";
	die($message."\n");
}

$/ = "\n";
$score1 = 0;
open(SAM,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($line = <SAM>) {
	if ( substr($line,0,1) eq '@' ) {
		next;
	}
	chomp($line);
	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual,$AS) = split("\t",$line);

	@f = split(':',$AS);
	$score1 += $f[2];
}
close(SAM);
print $score1;
