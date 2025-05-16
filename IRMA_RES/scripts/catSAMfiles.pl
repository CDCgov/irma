#!/usr/bin/env perl
# catSAMfiles
#
# Samuel S. Shepard - 2015
#
# Description: concatenate SAM records without the header

while ( $line = <> ) {
    print $line;
    if ( $line !~ /^@/ ) { last; }
}

while ( $line = <> ) {
    if ( $line !~ /^@/ ) {
        print $line;
    }
}
