#!/usr/bin/env perl
# Samuel Shepard - 3.14.2019
# Perform left join

use strict;
use warnings;

my ($delim,$fieldSet,$message);
use Getopt::Long;
GetOptions(
		'delim|D=s' => \$delim,
		'field|F=s' => \$fieldSet
	);

if ( scalar(@ARGV) < 1 ) {
	$message = "Usage:\n\tperl $0 <table> [<left1> <left2> ...]\n";
	$message .= "\t\t--delim|-D <CHAR>\tDelimiter for the key, the column delimiter many only be tab.\n";
	$message .= "\t\t--field|-F <STR>\tComma-delimited set of fields to use for group. Default: column 1.\n";
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


my $fileLimit = scalar(@ARGV) - 1;
my @data = (); 
my @lengthRemaining = ();
foreach my $i ( 1 .. $fileLimit) {
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
			my @remainingColumns = map { $_ eq '' ? '\N' : $_ } complementArray(\@values,\@fields);
			my $N = scalar(@remainingColumns);
			if ( !defined($lengthRemaining[$i-1]) || $N > $lengthRemaining[$i-1] ) {
				$lengthRemaining[$i-1] = $N;
			}
			$data[$i-1]{$id} = [@remainingColumns];
		}
	}
	close($IN);
}

open(my $IN,'<',$ARGV[0]) or die("Cannot open main table: $ARGV[0].\n");
while(my $line = <$IN>) {
	chomp($line);
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
						$line .= "\t\\N";
					}
				} else {
					$line .= "\t".join("\t",@{$data[$i-1]{$id}});
				}
			} else {
				$line .= "\t\\N" for 1 .. $lengthRemaining[$i-1];
			}
		}
	}	
	print STDOUT $line,"\n";
}
close($IN);
