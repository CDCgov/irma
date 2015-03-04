#!/usr/bin/env perl



$line=<>; print $line;
$line=<>; print $line;

while($line=<>) {
	if ( $line !~ /^@/ ) {
		print $line;
	}
}
