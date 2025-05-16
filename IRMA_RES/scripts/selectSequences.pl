#!/usr/bin/env perl
# selectSequences - Version 1.1
#
# Samuel S. Shepard - 2012
#
# Description: Subsets fasta data using inclusion or exclusion pattern searches (w.r.t. the header).

####################
# PROCESS parameters
use Getopt::Long;
GetOptions(
            'case-sensitive|C' => \$case,
            'include|I:s'      => \$include,
            'append|A:s'       => \$append,
            'exclude|E:s'      => \$exclude,
            'min-length|L=i'   => \$minLength,
            'show-excluded|S'  => \$showExcluded,
            'to-upper|U'       => \$toUpper
);

if ( defined($minLength) ) {
    $byLength = 1;
} else {
    $byLength = 0;
}

if ( -t STDIN && !scalar(@ARGV) ) {
    $message = "Usage:\n\tperl $0 [FASTA ...] [OPTIONS]\n";
    $message .= "\t\t-C|--case-sensitive\tMatches ignore case.\n";
    $message .= "\t\t-I|--include <STRING>\tInclude headers with given string.\n";
    $message .= "\t\t-E|--exclude <STRING>\tExclude headers with given string.\n";
    $message .= "\t\t-A|--append <STRING>\tAppend STRING to output headers.\n";
    $message .= "\t\t-L|--min-length <INT#>\tMinimum sequence length permitted.\n";
    $message .= "\t\t-S|--show-excluded\tPrint excluded sequences to STDERR.\n";
    $message .= "\t\t-U|--to-upper\t\tTo uppercase.\n";
    $message .= "\t<STDIN> read if no FASTA given.\n";

    die($message);
}

# PROCESS fasta data
$/ = ">";
while ( $record = <> ) {
    chomp($record);
    @lines    = split( /\r\n|\n|\r/, $record );
    $header   = $h = shift(@lines);
    $sequence = join( '', @lines );

    if ( $append ne '' ) {
        $header .= $append;
    }

    $length = length($sequence);
    if ( defined($toUpper) ) {
        $sequence = uc($sequence);
    }

    if ( $length == 0 ) {
        next;
    } elsif ( $byLength && $length < $minLength ) {
        next;

        # both
    } elsif ( $include ne '' && $exclude ne '' ) {
        if ( $h =~ /\Q$include\E/ && $h !~ /\Q$exclude\E/ ) {
            print '>', $header, "\n", $sequence, "\n";
        }

        # include
    } elsif ( $include ne '' && $exclude eq '' ) {
        if ( $h =~ /\Q$include\E/ ) {
            print '>', $header, "\n", $sequence, "\n";
        }

        # exclude
    } elsif ( $include eq '' && $exclude ne '' ) {
        if ( $h !~ /\Q$exclude\E/ ) {
            print '>', $header, "\n", $sequence, "\n";
        } elsif ($showExcluded) {
            print STDERR "Excluded ($exclude):\t", $header, "\n";
        }
    } else {
        print '>', $header, "\n", $sequence, "\n";
    }
}
####################
