#!/usr/bin/env perl

use English qw(-no_match_vars);
use warnings;
use strict;

open my $logfile, '<', 'CHANGELOG.md' or die "Cannot open 'CHANGELOG.md': $OS_ERROR\n";

local $RS = "## [";
<$logfile>;    # skip
my $changes = <$logfile>;
chomp($changes);

my $version = q{};
if ( $changes =~ /^([\d.]+)\]/smx ) {
    $version = $1;
}

print STDOUT '## Release Notes for v[', $changes, "\n";

local $RS = "\n";
while ( my $line = <$logfile> ) {
    if ( $line =~ /^\[$version\]/smx ) {
        print STDOUT $line, "\n";
        last;
    }
}

