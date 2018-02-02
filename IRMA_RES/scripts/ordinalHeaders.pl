#!/usr/bin/env perl
# Sam Shepard - 8.2015

use Getopt::Long;
GetOptions(	'keep-annotation|A'=> \$keepAnnot, 
		'annot-with-header|H' => \$headerAnnot,
		'name-annotation-format|N=s' => \$nameAnnotation,
		'only-name-format|O=s' => \$onlyName,
		'pipe-format|P' => \$pipeFormat
	);
open( IN, '<', $ARGV[0] ) or die("$0 ERROR: Cannot open $ARGV[0].\n");

if ( scalar(@ARGV) != 1 ) {
	die("Usage:\n\tperl $0 <input.fasta> [-A] [-H] [-N <STR>] [-O <STR>]\n");
}

# PROCESS fasta data
$/ = ">"; $i = 0;
while($record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$oldHeader = shift(@lines);
	$ordinal = 'S'. ($i++);
	$sequence = join('',@lines);

	if ( $keepAnnot ) {
		if ( $oldHeader =~ /{(.+)}/ ) {
			$ordinal .= "{$1}";
		}
	} elsif ( $headerAnnot ) {
		$ordinal .="{$oldHeader}";
	} elsif ( $nameAnnotation ) {
		$ordinal = $nameAnnotation.'{'.$ordinal.'}';
	} elsif ( $onlyName ) {
		$ordinal = $onlyName;
	} elsif ( $pipeFormat ) {
		$ordinal .= '|'.$oldHeader;
	}

	if ( length($sequence) == 0 ) {
		next;
	} else {
		print '>',$ordinal,"\n";
		print $sequence,"\n";
	}
}
close(IN);
####################
