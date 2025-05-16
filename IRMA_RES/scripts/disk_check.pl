#!/usr/bin/env perl
# disk_check.pl
#
# Samuel S. Shepard
#
# Description: check to see if we have enough disk space to process a run

use English qw(-no_match_vars);

if ( scalar(@ARGV) < 2 ) {
    die("Usage:\n\tperl $PROGRAM_NAME <PPATH> <file1.fastq.gz> <file2.fastq.gz> <... fastq.gz>\n");
}

$ppath = shift(@ARGV);
if ( -e $ppath ) {
    open( my $DF, "df -k \"$ppath\"|" ) or die("Cannot run df.\n");
    @lines = <$DF>;
    chomp(@lines);
    close($DF) or die("Could not close: $OS_ERROR\n");
    $available_k = ( split( /\s+/smx, $lines[$#lines] ) )[3];

    # K to M
    $available  = $available_k / 1024;
    $total_size = 0;
    foreach my $f (@ARGV) {
        if ( $f =~ /^\s*$/smx ) { next; }
        if ( -e $f ) {
            $total_size += ( stat($f) )[7];
        } else {
            die("$PROGRAM_NAME: bad file input $f\n");
        }
    }

    # convert bytes to Mebibytes
    $file_m = $total_size / 1024 / 1024;

    # based on empirical data
    $estimated_m = $file_m * 15 + 5;

    # if the available resource are not enough for what is estimated, complain
    if ( $available <= $estimated_m ) {
        print sprintf( "needed ~%.2fM to execute, but only %.2fM available on disk. The project path was: '$ppath'\n",
                       $estimated_m, $available );
        exit 3;
    } else {
        $units = 'M';
        if ( $available > 2048 ) {
            $available /= 1024;
            $units = 'G';
            if ( $available > 2048 ) {
                $available /= 1024;
                $units = 'T';
                if ( $available > 2048 ) {
                    $available /= 1024;
                    $units = 'P';
                }
            }
        }

        print sprintf( "found %.1f%s free space, only needed ~%.1fM\n", $available, $units, $estimated_m );
        exit 0;
    }
} else {
    die("$PROGRAM_NAME: bad path input $ppath\n");
}
