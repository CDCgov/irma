#!/usr/bin/env perl
## no critic (RequireExtendedFormatting, RequireDotMatchAnything)
# In general:
# - Extended formatting has comments and what we have is readable enough.
# - Our use of `.` matching should not consume newlines.

use English qw(-no_match_vars);
use Carp    qw(croak);
use strict;
use warnings;

local $RS = undef;
open my $fh, '<', 'CHANGELOG.md' or croak "Can't open CHANGELOG.md: $OS_ERROR";
my $changelog = <$fh>;
close $fh or croak "Can't close CHANGELOG.md: $OS_ERROR";

local $RS = "\n";
my $script_version = q{};
open my $HEAD, '<', 'IRMA' or die "Cannot open IRMA: $OS_ERROR!\n";
while ( my $line = <$HEAD> ) {
    if ( $line =~ /^VERSION="([^"]+)"$/m ) {
        $script_version = $1;
        last;
    }
}
close $HEAD or croak "Can't close IRMA: $OS_ERROR!\n";

if ( $changelog =~ /^## \[(\S+?)\]: (\S+?)$/m ) {
    my ( $version, $date ) = ( $1, $2 );

    if ( $date ne 'TBD' ) {
        if ( $date !~ /^20\d{2}-[01]\d-[0-3]\d$/m ) {
            die "Version $version has invalid date format: $date (expected yyyy-mm-dd or TBD)\n";
        }

        if ( $version ne $script_version ) {
            die "IRMA header ($script_version) mismatches changelog ($version / $date)!\n";
        }

        # Here we do need to consume newlines
        if ( $changelog !~ /<!-- Versions -->.*?^\[$version\]:/ms ) {
            die "Version $version was not linked in CHANGELOG!\n";
        }
    }
}
