#!/usr/bin/env perl
# Sam Shepard with Kristine Lacek
# 2020.02

use strict;
use warnings;

my ($expected_length, $trim_ends,$message,$remove_inserted_N,$remove_deletions);

use Getopt::Long;
GetOptions(
		'expected-length|L=i' => \$expected_length,
		'trim-ends|T' => \$trim_ends,
		'remove-deletions|D' => \$remove_deletions,
		'remove-inserted-N|I' => \$remove_inserted_N
	);

if ( scalar(@ARGV) != 3 && scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <a2m> <sam>\n";
	$message .= "\t-L|--expected-length <INT>\tExpected length of coordinate space for reference.\n";
	$message .= "\t-T|--trim-ends\t\t\tTrim padded sequence ends of '-' and 'N'.\n";
	$message .= "\t-D|--remove-deletions\t\tRemoves all deletions (-) from the padded sequence.\n";
	$message .= "\t-I|--remove-inserted-N\t\tRemove inserted Ns (n) from the padded sequence.\n";
	die($message."\n");
}


## Read A2M in FASTA format. 
$/ = '>';
my ($ref_header,$ref_sequence) = ('','');
open(A2M,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while(my $record = <A2M>) {
	chomp($record);
	my @lines = split(/\r\n|\n|\r/, $record);
	$ref_header = shift(@lines);
	$ref_sequence = join('',@lines);
	my $length = length($ref_sequence);
	if ( $length < 1 ) {
		next;
	}
}
close(A2M);

my $coordinate_reference = $ref_sequence;
$coordinate_reference =~ tr/atcgn.//d;
my $coordinate_length = length($coordinate_reference);
my $coordinate_mask = $coordinate_reference;
$coordinate_mask =~ tr/-/*/;

if ( defined($expected_length) && int($expected_length) != $coordinate_length ) {
	die("Reference length mismatch: $expected_length != $coordinate_length\n");
}

sub getLength($) {
	my $cigar = $_[0];
	my $length = 0;
	# Consider if 'N' states should be included or not for use-case
	while($cigar =~ /(\d+)([MIDNSHP])/g ) {
		if ( $2 eq 'M' || $2 eq 'D' ) {
			$length += $1;
		}
	}
	return $length;
}


$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while(my $line = <SAM>) {
	if ( substr($line,0,1) eq '@' ) {  next; }
	chomp($line);
	my ($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	if ( $cigar eq '' || $pos eq '' ) { next; }

	my $rStart	= $pos - 1;			# refrence position is 0-based
	my $rLength	= getLength($cigar);		# Length of 1 is a single position

	substr($coordinate_mask, $rStart, $rLength, substr($coordinate_reference, $rStart, $rLength) );	
}

# Restore insertion states and ambiguate amplicon drop out
my $padded_reference = $coordinate_mask;
while ( $ref_sequence =~ /([atcgn]+)/g ) {
	my $insert = $1;
	substr($padded_reference, $-[0], 0) = $1
}
while ( $padded_reference =~ /(?<=[*N])([ACTG])?[-]+|[-]+([ACTG])?(?=[*N])/g ) {
        my $replacement_string = '*' x length($0);
        substr($padded_reference, $-[0], length($0)) = $replacement_string
}
$padded_reference =~ tr/*/N/;


if ( defined($remove_inserted_N) ) {
	$padded_reference =~ tr/n//d;
}

if ( defined($remove_deletions) ) {
	$padded_reference =~ tr/.-//d;
}

if ( defined($trim_ends) ) {
	$padded_reference =~ s/^[N.]+//;
	$padded_reference =~ s/[N.]+$//;
}

print STDOUT '>',$ref_header,"\n",$padded_reference,"\n";
