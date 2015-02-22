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
	$notSkipExtension = 0;
} else {
	$notSkipExtension = 1;
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
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	} else {
		if ( $notSkipExtension ) {
			if ( $sequence =~ /^([actgn.]*)([-]{0,2}[ACTGN])/ ) {
				$leader = $1;
				$leader =~ tr/.//d;
				for $x ( - length($leader) .. -1 ) {
					$count[1]{$x}{substr($leader,$x,1)}++;
				}
			}

			if ( $sequence =~ /[ATCGN][-]{0,2}([actgn.]*)$/ ) {
				$trailer = $1;
				$trailer =~ tr/.//d;
				for $x ( 0..length($trailer)-1 ) {
					$count[2]{$x}{substr($trailer,$x,1)}++;
				}
			}
		} else {
			if ( $sequence =~ /^([actg]{1,2}[actgn.]*)([-]{1,2})[ACTGN]/ ) {
				$leader = $1; $gap = length($2);
				$leader =~ tr/.//d;
				if ( length($leader) >= $gap ) {
					substr($sequence,$+[2]-$gap,$gap) = uc(substr($leader,-$gap));
				}
			}

			if ( $sequence =~ /[ATCGN]([-]{1,2})([actg]{1,2}[actgn.]*)$/ ) {
				$trailer = $2; $gap = length($1);
				$trailer =~ tr/.//d;
				if ( length($trailer) >= $gap ) {
					substr($sequence,$-[1],$gap) = uc(substr($trailer,0,$gap));
				}
			}
		}

		$length -= ($sequence =~ tr/[a-z.]//d)+1;
		for $p ( 0 .. $length) {
			 $count[0][$p]{substr($sequence,$p,1)}++;
		}
		$i++;
	}
}
close(IN);
store(\@count, $ARGV[1]);
