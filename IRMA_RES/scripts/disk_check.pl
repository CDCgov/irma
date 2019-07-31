#!/usr/bin/env perl
# Sam Shepard

if ( scalar(@ARGV) < 2 ) {
	die("Usage:\n\tperl $0 <PPATH> <file1.fastq.gz> <file2.fastq.gz> <... fastq.gz>\n");
}

$ppath = shift(@ARGV);
if ( -e $ppath ) {
	open(DF, "df -k \"$ppath\"|") or die("Cannot run df.\n");
	@lines = <DF>; chomp(@lines); close(DF);
	$available_k = ( split(/\s+/,$lines[$#lines] ) )[3] ;
	# K to M
	$available = $available_k / 1024;		
	$total_size = 0;
	foreach $f ( @ARGV ) {
		if ( $f =~ /^\s*$/ ) { next; }
		if ( -e $f ) {
			$total_size += ( stat($f) )[7];
		} else {
			die("$0: bad file input $f\n");
		}
	}
	# convert bytes to Mebibytes
	$file_m = $total_size / 1024 / 1024;
	# based on empirical data
	$estimated_m = $file_m * 15 + 5;

	# if the available resource are not enough for what is estimated, complain
	if ( $available <= $estimated_m ) {
		print sprintf("needed ~%.2fM to execute, but only %.2fM available on disk\n",$estimated_m,$available);
		exit 3;
	} else {
		$units = 'M';
		if ( $available > 2048 ) { 
			$available /= 1024; $units = 'G';
			if ( $available > 2048 ) {
				$available /= 1024; $units = 'T';
				if ( $available > 2048 ) {
					$available /= 1024; $units = 'P';
				}
			}
		}
		
		print sprintf("found %.1f%s free space, only needed ~%.1fM\n",$available,$units,$estimated_m);
		exit 0;
	}
} else {
	die("$0: bad path input $ppath\n");
}
