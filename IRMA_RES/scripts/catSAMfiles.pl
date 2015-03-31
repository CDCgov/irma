#!/usr/bin/env perl
# Sam Shepard - 3.16.2015

while($line=<>) {
	print $line;
	if ( $line !~ /^@/ ) { last; }	
}

while($line=<>) {
	if ( $line !~ /^@/ ) {
		print $line;
	}
}
