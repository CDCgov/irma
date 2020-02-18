#!/usr/bin/env perl
# alignPositionByPSA
# Samuel Shepard - 2018.10
# Version 1.0

use Getopt::Long;
GetOptions(	'no-header|H' => \$skipHeader,
		'table-field|F=i' => \$field,
		'table-delim|D=s' => \$delim,
		'ref-name|N=s' => \$name,
		'prefix|P=s' => \$prefix,
		'inserts-to-ref|I' => \$insertsToRef,
		'show-original-position|O' => \$showOriginal,
		'skip-comments|S' => \$skipComments,
		'reprint-comments|R' => \$reprintComments
	);

if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <pairwiseAlignment.fasta> <table> [options]\n";
	$message .= "\t\t-H|--no-header\t\t\tDo not print input table header line. (first line)\n";
	$message .= "\t\t-F|--table-field <INT>\t\tField number in table (1-base) containing positional information. DEFAULT = 1\n";
	$message .= "\t\t-D|--table-delim <CHA>\t\tField delimiter for table. DEFAULT = <TAB>\n";
	$message .= "\t\t-N|--ref-name <STR>\t\tName of the reference sequence. DEFAULT = first record\n";
	$message .= "\t\t-P|--prefix <STR>\t\tPrefix to table.\n";
	$message .= "\t\t-I|--inserts-to-ref\t\tInserts relative to reference.\n";
	$message .= "\t\t-O|--show-original-position\tShow original position column.\n";
	$message .= "\t\t-S|--skip-comments\t\tSkip lines beginning with #.\n";
	$message .= "\t\t-R|--reprint-comments\t\tSimply reprint lines beginning with #.\n";
	die($message."\n");
}


if ( defined($skipComments) || defined($reprintComments) ) {
	$checkComments = 1;
} else {
	$checkComments = 0;
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
if ( $checkComments && substr($header,0,1) eq '#' ) {
	if ( $reprintComments ) {
		print $header;
	}
} else {
	if ( !$skipHeader) {
		if ( defined($showOriginal) ) {
			chomp($header);
			@H = split($delim,$header);
			$H[$field] = $H[$field].'_original'.$delim.$H[$field].'_revised';
			print $prefix,join($delim,@H),"\n";
		} else {
			print $header;
		}
	}
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
	if ( $checkComments && substr($line,0,1) eq '#' ) {
		if ( $reprintComments ) {
			print $line,"\n";
		}
		next;
	}	
	@fields = split($delim,$line);
	$aCoord = $fields[$field];

	if ( !defined($refByAlt{$aCoord}) ) {
		print STDERR $aCoord," not found\n";
		next;
	} else {
		$rCoord = $refByAlt{$aCoord};
	}

	if ( $removeInserts && $rCoord =~ /\./) { next; }

	$fields[$field] = defined($showOriginal) ? $aCoord.$delim.$rCoord : $rCoord;
	print $prefix,join($delim,@fields),"\n";
}
