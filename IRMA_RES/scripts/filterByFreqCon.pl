#!/usr/bin/env perl
# Sam Shepard

if ( scalar(@ARGV) != 2 ) {
	$message = "perl $0 <all> <vars>";
	die("$message\n");
}

open(IN,'<',$ARGV[0]) or die("cannot open $ARGV[0] for reading.\n");
$allHdr = <IN>;
@allRecs = <IN>;
close(IN);

open(IN,'<',$ARGV[1]) or die("cannot open $ARGV[1] for reading.\n");
$varHdr = <IN>;
@varRecs = <IN>;
close(IN);

$max = 0;
foreach $allRec ( @allRecs ) {
	chomp($allRec);
	@fields = split("\t",$allRec);
	($frequency,$confidence,$type) = ($fields[5],$fields[7],$fields[10]);
	if ( $confidence eq 'NA' ) { next; }
	if ( $confidence == 0 && $frequency > $max ) {
		$max = $frequency;
	}
}

print $varHdr;
foreach $varRec ( @varRecs ) {
	chomp($varRec);
	@fields = split("\t",$varRec);
	$frequency = $fields[8];
	if ( $frequency > $max ) {
		print $varRec,"\n";	
	}
}

