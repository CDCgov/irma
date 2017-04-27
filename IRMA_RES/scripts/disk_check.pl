#!/usr/bin/env perl
# Sam Shepard

if ( scalar(@ARGV) < 2 ) {
	die("Usage:\n\tperl $0 <PPATH> <file1.fastq.gz> <file2.fastq.gz> <... fastq.gz>\n");
}

$ppath = shift(@ARGV);
if ( -e $ppath ) {
	open(DF, "df -k $ppath|") or die("Cannot run df.\n");
	@lines = <DF>; chomp(@lines); close(DF);
	$available_k = ( split(/\s+/,$lines[$#lines] ) )[3] ;
	$available_m = $available_k / 1024;
	$total_size = 0;
	foreach $f ( @ARGV ) {
		if ( -e $f ) {
			$total_size += ( stat($f) )[7];
		}
	}
	$file_m = $total_size / 1024 / 1024;
	$estimated_m = $file_m * 5 + 10;

	if ( $available_m < $estimated_m ) {
		print STDERR "WARNING! Need around $estimated_m to complete the project but only have $available_m available on disk.\n";
	}
	print $total_size,"\t",$available_k,"\n"; exit 0;
} else {
	print "-1"; exit 1;
}
