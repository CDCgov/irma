#!/usr/bin/env perl
# makePluralityConsensus
# Sam Shepard - vfn4@cdc.gov - 2016-11-04
# This script is in the public domain.

# DICLAIMER & LIMITATION OF LIABILITY.
# The materials embodied in this software are "as-is" and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose. In no event shall the Centers for Disease Control and Prevention (CDC) or the United States (U.S.) Government be liable to you or anyone else for any direct, special, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not CDC or the U.S. Government has been advised of the possibility of such loss, however caused and on any theory of liability, arising out of or in connection with the possession, use or performance of this software. In no event shall any other party who modifies and/or conveys the program as permitted according to GPL license [www.gnu.org/licenses/], make CDC or the U.S. government liable for damages, including any general, special, incidental or consequential damages arising out of the use or inability to use the program, including but not limited to loss of data or data being rendered inaccurate or losses sustained by third parties or a failure of the program to operate with any other programs. Any views, prepared by individuals as part of their official duties as United States government employees or as contractors of the United States government and expressed herein, do not necessarily represent the views of the United States government. Such individuals' participation in any part of the associated work is not meant to serve as an official endorsement of the software. The CDC and the U.S. government shall not be held liable for damages resulting from any statements arising from use of or promotion of the software that may conflict with any official position of the United States government.

use Getopt::Long;
GetOptions(	'name|N=s' => \$name, 'allow-deletions|D' => \$allowDeletions );

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <input.fasta> [-N|--name <STRING>] [-D|--allow-deletions]\n";
	die($message."\n");
}

# PROCESS fasta data
open(IN, '<', $ARGV[0]) or die("Cannot open $ARGV[0].\n");
$/ = ">"; %count = (); $L = 0;
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = uc(join('',@lines));
	$length = length($sequence);

	if ( $length == 0 ) {
		next;
	}

	# unify gap characters
	$sequence =~ tr/:~./-/;
	@col = split('', $sequence);
	foreach $p ( 0..($length-1) ) {
		$count[$p]{$col[$p]}++;
	}
}
close(IN);
$L = $length;

if ( $name ) {
	print '>',$name,"\n";
} else {
	print ">consensus\n";
}

foreach $p ( 0..($length-1) ) {
	@nts = sort { $count[$p]{$b} <=> $count[$p]{$a} } keys( %{$count[$p]} );
	foreach $nt ( 0..$#nts ) {
		if ( $nts[$nt] ne '-' || $allowDeletions ) {
			print $nts[$nt];
			last;
		}
	}
}
print "\n";
