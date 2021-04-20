#!/usr/bin/env perl
# Samuel Shepard - 2021-03
# Perform left join

use strict;
use warnings;

my ($delim,$fieldSet,$message,$flexibleHeader,$inPlace,$impalaNull);
use Getopt::Long;
GetOptions(
		'delim|D=s' => \$delim,
		'field|F=s' => \$fieldSet,
		'flexible-header|H' => \$flexibleHeader,
		'in-place|I' => \$inPlace,
		'impala-null|N' => \$impalaNull
	);


if ( scalar(@ARGV) < 1 ) {
	$message = "\nUsage:\n\tperl $0 <main_table> [<join1> <join2> ...]\n";
	$message .= "\t\t-D|--delim <CHAR>\tDelimiter for the key, the column delimiter many only be tab.\n";
	$message .= "\t\t-F|--field <STR>\tComma-delimited set of fields to use for group. Default: column 1.\n";
	$message .= "\t\t-H|--flexible-header\tAllows header group fields to not match.\n";
	$message .= "\t\t-I|--in-place\t\tOverwrite main table using joined table.\n";
	$message .= "\t\t-N|--impala-null\tUse Apache Impala style nulls (\\N) rather than R style nulls (NA).\n";
	die($message."\n");
}

sub complementArray($$) {
	my ($A1,$A2) = ($_[0],$_[1]);
	my %H2 = map { $_ => 1 } @{$A2};
	my $key = '';

	my @indices = grep { !defined($H2{$_}) } (0..$#{$A1});
	return ( @{$A1}[@indices] );
}

my @fields = (); 
my $numberSelected = 0; 
my $maxSelected = 1;
if ( defined($fieldSet) ) {
	@fields = split(',', $fieldSet);
	$numberSelected = scalar(@fields);
	foreach my $x (@fields ) {
		if ( $x > $maxSelected ) { $maxSelected = $x; }
		if ( $x == 0 ) {
			die("$0 ERROR: field must be specified.\n");
		} elsif( $x < 0 ) {
			die("$0 ERROR: field must be a positive number.\n");
		}
	}
	foreach my $x ( 0 .. ($numberSelected-1) ) { 
		$fields[$x]--; 
	}
} else {
	$fields[0] = 0;
}

if ( !defined($delim) ) {
	$delim = '|';
} elsif( $delim eq '' ) {
	die("$0 ERROR: No delimiter argument detected.\n");
} elsif( length($delim) > 1 ) {
	die("$0 ERROR: single character delimiter expected instead of '$delim'.\n");
}

$flexibleHeader = defined($flexibleHeader) ? 1 : 0;
$inPlace = defined($inPlace) ? 1 : 0;
my $NULL = defined($impalaNull) ? '\N' : 'NA';

my $fileLimit = scalar(@ARGV) - 1;
my @data = (); 
my @lengthRemaining = ();
my @header = ();
foreach my $i ( 1 .. $fileLimit ) {
	my $first = 1;
	open(my $IN,'<',$ARGV[$i]) or die("Cannot open file $ARGV[$i].\n");
	while(my $line = <$IN> ) {
		chomp($line);
		my @values = split("\t",$line);
		my $id = '';
		my $numberFound = scalar(@values);
		if ( $numberSelected > 0 ) {
			if ( $maxSelected > $numberSelected ) {
				die("$0 ERROR: non-existant field specified. Wanted $numberSelected (max: $maxSelected) but found $numberFound\n");
			}
			$id = join($delim, (@values[@fields]) );
		} else {
			$id = $values[ $fields[0] ];
		}

		if ( $id ne '' ) {
			my @remainingColumns = map { $_ eq '' ? $NULL : $_ } complementArray(\@values,\@fields);
			my $N = scalar(@remainingColumns);
			if ( !defined($lengthRemaining[$i-1]) || $N > $lengthRemaining[$i-1] ) {
				$lengthRemaining[$i-1] = $N;
			}
			$data[$i-1]{$id} = [@remainingColumns];
			if ( $first ) {
				$first = 0;
				$header[$i-1] = [@remainingColumns];
			}
		}
	}
	close($IN);
}

sub leftJoin($$) {
	my $line = $_[0]; chomp($line);
	my $first = $_[1];
	my @values = split("\t",$line);
	my $numberFound = scalar(@values);
	my $id = '';
	if ( $numberSelected > 0 ) {
		if ( $maxSelected > $numberSelected ) {
			die("$0 ERROR: non-existant field specified. Wanted $numberSelected (max: $maxSelected) but found $numberFound\n");
		}
		$id = join($delim, (@values[@fields]) );
	} else {
		$id = $values[ $fields[0] ];
	}

	if ( $id ne '' ) {
		foreach my $i ( 1 .. $fileLimit ) {
			if ( defined($data[$i-1]{$id}) ) {
				my $N = scalar(@{$data[$i-1]{$id}});
				if ( $N < $lengthRemaining[$i-1] ) {
					if ( $N > 0 ) {
						$line .= "\t".join("\t",@{$data[$i-1]{$id}});
					}
					foreach( 1 .. ($lengthRemaining[$i-1] - $N) ) {
						$line .= "\t$NULL";
					}
				} else {
					$line .= "\t".join("\t",@{$data[$i-1]{$id}});
				}
			} elsif ( $first && $flexibleHeader ) {
				$line .= "\t".join("\t",@{$header[$i-1]});
			} else {
				
				$line .= "\t$NULL" for 1 .. $lengthRemaining[$i-1];
			}
		}
	}	

	return $line;
}


my $first = 1;
if ( $inPlace ) {
	open(my $IN,'<',$ARGV[0]) or die("Cannot open main table for reading: $ARGV[0].\n");
	my @lines = <$IN>;
	close($IN);

	open(my $OUT,'>',$ARGV[0]) or die("Cannot open main table for writing: $ARGV[0].\n");
	if ( scalar(@lines) > 0 ) {
		print $OUT leftJoin($lines[0],1),"\n";
	}

	if ( scalar(@lines) > 1 ) {
		foreach my $i ( 1 .. $#lines ) {
			print $OUT leftJoin($lines[$i],0),"\n";
		}
	}
	close($OUT);
} else {
	open(my $IN,'<',$ARGV[0]) or die("Cannot open main table for reading: $ARGV[0].\n");
	print STDOUT leftJoin(<$IN>,1),"\n";
	while(my $line = <$IN>) {
		print STDOUT leftJoin($line,0),"\n";
	}
	close($IN);
}
