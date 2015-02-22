#!/usr/bin/env perl
# Sam Shepard - 2.2014

use Getopt::Long;
GetOptions(	'name|N=s' => \$name
	);

if ( scalar(@ARGV) < 1 ) {
	$message = "Usage:\n\tperl $0 <input.fasta>\n";
	die($message."\n");
}

# PROCESS fasta data
$/ = ">"; $i = 0; %count = (); $L = 0;
while($record = <> ) {
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
	for($j = 0; $j < $length; $j++) {
		if ( defined($count[$j]{$col[$j]}) ) {
			 $count[$j]{$col[$j]}++;
		} else {
			 $count[$j]{$col[$j]} = 1;
		}
	}
	$i++;
}
$L = $length;

if ( $name ) {
	print '>',$name,"\n";
} else {
	print ">consensus\n";
}
for($j = 0; $j < $L; $j++ ) {
	@nts = sort { $count[$j]{$b} <=> $count[$j]{$a} } keys( %{$count[$j]} );
	for($i = 0; $i < scalar(@nts); $i++) {
		if ( $nts[$i] ne '-' ) {
			print $nts[$i];
			last;
		}
	}
}
print "\n";
