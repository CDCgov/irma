#!/usr/bin/env perl
# Sam Shepard - 5.2015
# Accumulator for delimited data.

use Getopt::Long;
GetOptions(
		'ignore-header|H' => \$header,
		'comment-char|C=s' => \$cc,
		'field-number|F=i' => \$f,
		'field-delim|D=s' => \$d,
		'average|A' => \$takeAverage,
		'sum-count|S' => \$sumAndCount,
		'take-count|T' => \$takeCount,
		'newline|N' => \$newline,
		'minimum|M' => \$takeMinimum,
		'maximum|X' => \$takeMaximum,
		'everything|E' => \$takeEverything,
		'help|?' => \$printHelp
	);

if ( $printHelp ) {
		$message = "Usage:\nsumField.pl <file or STDIN> [options]\n";
		$message .= "\t\t-H|--ignore-header|H\t\tSkips header line for calculations.\n";
		$message .= "\t\t-C|--comment-char <CHAR>\tSkip lines beginning with the comment field.\n";
		$message .= "\t\t-F|--field-number <NUM>\t\tField number (1..N) to calculate on.\n";
		$message .= "\t\t-D|--field-delim <CHAR>\t\tDelimiter for fields, default <tab>.\n";
		$message .= "\t\t-A|--average\t\t\tTake the average of the column, default sum.\n";
		$message .= "\t\t-S|--sum-count\t\t\tCalculate the sum and the count.\n";
		$message .= "\t\t-C|--take-count\t\t\tTake the count of the column.\n";
		$message .= "\t\t-N|--newline\t\t\tPrint trailing newline.\n";
		$message .= "\t\t-M|--minimum\t\t\tTake the minimum of the column.\n";
		$message .= "\t\t-X|--maximum\t\t\tTake the maximum of the column.\n";
		$message .= "\t\t-E|--everything\t\t\tTake the count, sum, min, max, and average.\n";
		die($message);
}

if ( !defined($newline) ) {
	$newline = '';
} else {
	$newline = "\n";
}

if ( !defined($d) ) {
	$d = "\t";
} else {
	if ( $d eq '|' ) {
		$d = "/\|/";
	}
}

if ( !defined($f) || int($f) <= 0 ) {
	$f = 0;
} else {
	$f = int($f)-1;
}

$/ = "\n";
if ( $header ) { $junk=<>; }

$N = $total = 0;
my $min;
my $max;

if ( $sumAndCount ) {
	$takeAverage = 1;
}


if ( $takeEverything ) {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		}
		@fields = split($d,$line);
		if ( $f > $#fields ) {
			next;
		}

		if ( !defined($max) || $fields[$f] > $max ) {
			$max = $fields[$f];
		}

		if ( !defined($min) || $fields[$f] < $min ) {
			$min = $fields[$f];
		}

		$total += $fields[$f];
		$N++;
	}
	print $N,"\t",$total,"\t",$min,"\t",$max,"\t",sprintf("%.2f",$total/$N),$newline;
} elsif ( $takeMaximum ) {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		}
		@fields = split($d,$line);
		if ( $f > $#fields ) {
			next;
		}

		if ( !defined($max) || $fields[$f] > $max ) {
			$max = $fields[$f];
		}
	}
	print $max,$newline;
} elsif ( $takeMinimum ) {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		}
		@fields = split($d,$line);
		if ( $f > $#fields ) {
			next;
		}

		if ( !defined($min) || $fields[$f] < $min ) {
			$min = $fields[$f];
		}
	}
	print $min,$newline;
} elsif ( $takeAverage ) {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		}
		@fields = split($d,$line);
		if ( $f > $#fields ) {
			next;
		}

		$total += $fields[$f];
		$N++;
	}

	if ( $sumAndCount ) {
		print sprintf("%d\t%d",$total,$N),$newline;
	} else {
		print sprintf("%.2f",$total/$N),$newline;
	}
} elsif ( $takeCount ) {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		} else {
			$N++;
		}
	}
	print $N,$newline;
} else {
	while($line=<>) {
		chomp($line);
		if ( defined($cc) && $line =~ /^\Q$cc\E/ ) {
			next;
		}
		@fields = split($d,$line);
		if ( $f > $#fields ) {
			next;
		}

		$total += $fields[$f];
	}
	print $total,$newline;
}
