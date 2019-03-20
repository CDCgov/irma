#!/usr/bin/env perl
# makePatchworkConsensus
# Sam Shepard - vfn4@cdc.gov - 2019-03
# This script is in the public domain.

# DICLAIMER & LIMITATION OF LIABILITY.
# The materials embodied in this software are "as-is" and without warranty of any kind, express, implied or otherwise, including without limitation, any warranty of fitness for a particular purpose. In no event shall the Centers for Disease Control and Prevention (CDC) or the United States (U.S.) Government be liable to you or anyone else for any direct, special, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including without limitation, loss of profit, loss of use, savings or revenue, or the claims of third parties, whether or not CDC or the U.S. Government has been advised of the possibility of such loss, however caused and on any theory of liability, arising out of or in connection with the possession, use or performance of this software. In no event shall any other party who modifies and/or conveys the program as permitted according to GPL license [www.gnu.org/licenses/], make CDC or the U.S. government liable for damages, including any general, special, incidental or consequential damages arising out of the use or inability to use the program, including but not limited to loss of data or data being rendered inaccurate or losses sustained by third parties or a failure of the program to operate with any other programs. Any views, prepared by individuals as part of their official duties as United States government employees or as contractors of the United States government and expressed herein, do not necessarily represent the views of the United States government. Such individuals' participation in any part of the associated work is not meant to serve as an official endorsement of the software. The CDC and the U.S. government shall not be held liable for damages resulting from any statements arising from use of or promotion of the software that may conflict with any official position of the United States government.

use strict;
use warnings;
use Getopt::Long;

my ($name,$message);
GetOptions(	'name|N=s' => \$name );

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <a2m> <keyword_for_pairing>\n";
	die($message."\n");
}

# PROCESS a2m data
my @count = (); 
my @default = ();
my @lines = ();
my @col = ();
my @nts = ();
my ($header,$sequence,$length,$nt,$p);
my $L = 0;
my $keyword = $ARGV[1];

$/ = ">"; 
open(IN, '<', $ARGV[0]) or die("Cannot open $ARGV[0].\n");
while(my $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = join('',@lines);
	$length = length($sequence);
	if ( $length == 0 ) { next; }

	@col = split('', $sequence);
	if ( $header =~ /\Q$keyword\E/ ) {
		foreach $p ( 0..($length-1) ) {
			$default[$p]{$col[$p]}++;
		}

	} else {
		foreach $p ( 0..($length-1) ) {
			$count[$p]{$col[$p]}++;
		}
	}
}
close(IN);
$L = $length;

if ( defined($name) && $name ne '' ) { 
	print STDOUT '>',$name,"\n";
} else { 
	print STDOUT ">consensus\n";
}

$sequence = '';
foreach my $p ( 0..($length-1) ) {
	@nts = sort { $count[$p]{$b} <=> $count[$p]{$a} } keys( %{$count[$p]} );
	$nt = defined($nts[0]) ? $nts[0] : '.';
	if ( $nt =~ /[.-]/ ) {
		@nts = sort { $default[$p]{$b} <=> $default[$p]{$a} } keys( %{$default[$p]} );
		$nt = defined($nts[0]) ? $nts[0] : '.';
	}
	$sequence .= uc($nt);
}
$sequence =~ tr/.-//d;
print STDOUT $sequence,"\n";
