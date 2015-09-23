#!/usr/bin/env perl

use Getopt::Long;
GetOptions(	'no-header|H' => \$skipHeader,
		'table-field|F=i' => \$field,
		'table-delim|D=s' => \$delim,
		'ref-name|N=s' => \$name,
		'prefix|P=s' => \$prefix,
		'inserts-to-ref|I' => \$insertsToRef
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <pairwiseAlignment> <table> [options]\n";
	$message .= "\t\t--table-field|-F <INT>\tField number in table (1-base) containing positional information. DEFAULT = 1\n";
	$message .= "\t\t--table-delim|-D <CHA>\tField delimiter for table. DEFAULT = <TAB>\n";
	$message .= "\t\t--ref-name|-N <STR>\tName of the reference sequence. DEFAULT = first record\n";
	$message .= "\t\t--prefix|-P <STR>\tPrefix to table.\n";
	$message .= "\t\t--inserts-to-ref|-I\tInserts relative to reference.\n";

	die($message."\n");
}

if ( !defined($delim) ) {
	$delim = "\t";
}

if ( !defined($field) ) {
	$field = 0;
} elsif ( $field > 0 ) {
	$field = int($field) - 1;
} else {
	$field = 0;
}

if ( !defined($name) ) {
	$firstRecord = 1;
} else {
	$firstRecord = 0;
}

if ( defined($prefix) ) {
	$prefix = $prefix . "\t";
} else {
	$prefix = '';
}

open(IN,'<',$ARGV[0]) or die("ERROR: Cannot open $ARGV[0] for reading.\n");
$/ = ">"; $K = 0;
%matchByRef = ();
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = uc(join('',@lines));
	$length = length($sequence);
	if ( $length == 0 ) {
		next;	
	} else {
		$K++;
	}

	if ( ($firstRecord && $K == 1) || (!$firstRecord && $header =~ /$name/) ) {
		$refSeq = $sequence;
		$refHeader = $header;
		$refLength = $length;
	} else {
		$altSeq = $sequence;
		$altHeader = $header;
		$altLength = $length;
	}
}
close(IN);

if ( $refLength != $altLength ) {
	die("ERROR $0: Unequal lengths ($refLength <> $altLength).\n");
}

@ref = split('',$refSeq);
@alt = split('',$altSeq);
$aCoord = $rCoord = 0;
for($i=0;$i<$refLength;$i++) {
	$r = $ref[$i];
	$a = $alt[$i];

	if ( $r ne '-' ) {
		$insertIndex = 0;
		$rCoord++;
	}

	if ( $a ne '-' ) {
		$aCoord++;
	}

	# MATCH
	if ( $a ne '-' && $r ne '-' ) {
		$refByAlt{$aCoord} = $rCoord;
	# INSERTION relative to reference
	} elsif ( $a ne '-' ) {
		$insertIndex++;
		$refByAlt{$aCoord} = sprintf('%d.%04d',$rCoord,$insertIndex);
	}
}

open(IN,'<',$ARGV[1]) or die("ERROR: Cannot open $ARGV[1] for reading.\n");
$/ = "\n";
$header = <IN>;
if ( !$skipHeader) {
	print $header;
}
@data = <IN>; 
chomp(@data);
close(IN);

if ( $insertsToRef ) {
	$removeInserts = 0;
} else {
	$removeInserts = 1;
}

$aCoord = $rCoord = 0;
foreach $line (@data) {
	@fields = split($delim,$line);
	$aCoord = $fields[$field];

	if ( !defined($refByAlt{$aCoord}) ) {
		print STDERR $aCoord,"\n";
	}
	$rCoord = $refByAlt{$aCoord};
	if ( $removeInserts && $rCoord =~ /\./) {
		next;		
	}

	$fields[$field] = $rCoord;
	print $prefix,join($delim,@fields),"\n";
}
