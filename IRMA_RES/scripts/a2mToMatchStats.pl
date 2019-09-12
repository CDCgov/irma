#!/usr/bin/env perl
# a2mToMatchStats.pl
# Sam Shepard - 8.2014


use Getopt::Long;
GetOptions(	'skip-elongation|S' => \$skipExtension );
use Storable;
if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\t$0 <A2M> <OUT> [options]\n";
	$message .= "\t\t-S|--skip-elongation\tSkip the reference elongation algorithm.\n";
	die($message."\n");
}

if ( defined($skipExtension) ) {
	$elongateReference = 0;
} else {
	$elongateReference = 1;
}

# PROCESS fasta data
$/ = ">";
$i = 0; 
@count = (); 
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].");
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$sequence =~ tr/.//d;
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	} else {
		$leader = $trailer = "";
		if ( $sequence =~ /^([actgn]+)[ACTGN-]/ ) {
			$leader = $1;
		}

		if ($sequence =~ /[ACTGN-]([actgn]+)$/ ) {
			$trailer = $1;
		}
		
		$sequence =~ tr/[a-z]//d;
		$length = length($sequence);

		# repair as necessary
		$trailerLen = length($trailer);
		$leaderLen = length($leader);
		if ( $leaderLen && $sequence =~ /^([-]{1,9})[ACTGN]/ ) {
			$gap = length($1);
			if ( $leaderLen >= $gap ) {
				substr($sequence,0,$gap) = uc(substr($leader,-$gap));
				$leader = substr($leader,0,$leaderLen-$gap);
			}
		}

		if ( $trailerLen > 0 && $sequence =~ /[ACTGN]([-]{1,9})$/ ) {
			$gap = length($1);
			if ( $trailerLen >= $gap ) {
				substr($sequence,$-[1],$gap) = uc(substr($trailer,0,$gap));
				$trailer = substr($trailer,-($trailerLen-$gap));
			}
		}

		for $p ( 0 .. ($length-1) ) {
			 $count[0][$p]{substr($sequence,$p,1)}++;
		}

		if ( $elongateReference ) {
			if ( $sequence =~ /^[ACTGN]/ ) {
				for $x ( - length($leader) .. -1 ) {
					$count[1]{$x}{substr($leader,$x,1)}++;
				}
			}

			if ( $sequence =~ /[ATCGN]$/ ) {
				for $x ( 0..length($trailer)-1 ) {
					$count[2]{$x}{substr($trailer,$x,1)}++;
				}
			}
		}
		$i++;
	}
}
close(IN);
store(\@count, $ARGV[1]);
