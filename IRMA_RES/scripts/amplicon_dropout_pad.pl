#!/usr/bin/env perl
# Sam Shepard with Kristine Lacek
# 2020.03

use strict;
use warnings;

my ($expected_length, $trim_ends,$message,$remove_inserted_N,$remove_deletions,$pad_adjacent,$coverage_in,$coverage_out);

use Getopt::Long;
GetOptions(
		'expected-length|L=i' => \$expected_length,
		'trim-ends|T' => \$trim_ends,
		'remove-deletions|D' => \$remove_deletions,
		'remove-inserted-N|I' => \$remove_inserted_N,
		'pad-adjacent-deletions|A' => \$pad_adjacent,
		'a2m-coverage-in|C=s' => \$coverage_in,
		'pad-coverage-out|O=s' => \$coverage_out
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <a2m> <sam> [options]\n";
	$message .= "\t-L|--expected-length <INT>\tExpected length of coordinate space for reference.\n";
	$message .= "\t-T|--trim-ends\t\t\tTrim padded sequence ends of '-' and 'N'.\n";
	$message .= "\t-D|--remove-deletions\t\tRemoves all deletions (-) from the padded sequence.\n";
	$message .= "\t-I|--remove-inserted-N\t\tRemove inserted Ns (n) from the padded sequence.\n";
	$message .= "\t-A|--pad-adjacent-deletions\tPads deletions adjacent to dropouts.\n";
	$message .= "\t-C|--a2m-coverage-in <file>\tA2M coverage table..\n";
	$message .= "\t-O|--pad-coverage-out <file>\tOutput padded coverage table..\n";
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

if ( defined($pad_adjacent) ) {
	# Safely mask patterns found adjacent
	sub mask($$$) {
		my $all = defined($_[0]) ? $_[0] : '';
		my $left = defined($_[1]) ? $_[1] : '';
		my $right = defined($_[2]) ? $_[2] : '';
		return ($left . '*' x (length($all)-length($left)-length($right)) . $right);
	}

	$padded_reference =~ s/(?<=[*N])([ACTG])?[-]+|[-]+([ACTG])?(?=[*N])/mask($&,$1,$2)/ige;
	$padded_reference =~ s/^([*]+)/'-' x length($1)/e;
	$padded_reference =~ s/([*]+)$/'-' x length($1)/e;
}

if ( defined($coverage_in) ) {
	open(COV,'<',$coverage_in) or die("Cannot open $coverage_in for reading.\n");
	my $header = <COV>;
	my @coverage = <COV>; 
	chomp(@coverage);
	close(COV);

	if ( !defined($coverage_out) ) {
		$coverage_out = $coverage_in;
	}
	open(COV,'>',$coverage_out) or die("Cannot open $coverage_out for writing.\n");
	print COV $header;
	my ($iPos, $iCon) = (1, 3);
	foreach my $idx ( 0 .. $#coverage ) {
		my $line = $coverage[$idx];
		my @fields = split("\t",$line);

		# Will not line up with insertions
		# my $cPos = $fields[$iPos] - 1;
		if ( substr($padded_reference,$idx,1) eq '*' ) {
			if ( $fields[$iCon] eq '-' ) {
				print COV $fields[0],"\t",$fields[1],"\t0\tN\t0\t0\t0\t0\t$fields[8]\tP\n";
			} else {
				print STDERR "Expected '-' in table at $fields[$iPos], but found $fields[$iCon] instead, printing existing line.\n";
				print COV $line,"\n";
			}
		} else {
			print COV $line,"\n";
		}
	}	
	close(COV);
}

$padded_reference =~ tr/*/N/;

if ( defined($remove_inserted_N) ) {
	$padded_reference =~ tr/n//d;
}

if ( defined($remove_deletions) ) {
	$padded_reference =~ tr/.-//d;
}

if ( defined($trim_ends) ) {
	$padded_reference =~ s/^[N.]+(?=[^N.])//;
	$padded_reference =~ s/(?<=[^N.])[N.]+$//;
}

print STDOUT '>',$ref_header,"\n",$padded_reference,"\n";
