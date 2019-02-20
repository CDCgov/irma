#!/usr/bin/env perl
# Sam Shepard
# Replace complicated code from IRMA main script.

use strict;
use warnings;

use File::Basename;
use Getopt::Long;
my $paired;
GetOptions( 'paired-end|P' => \$paired );


sub sum_field($$) {
	my @lines = @{$_[0]};
	my $field = $_[1];
	my $sum = 0;
	my @f = ();
	foreach my $line ( @lines ) {
		@f = split("\t",$line);
		$sum += $f[$field];
	} 
	return $sum;
}

sub sum($$$) {
	my @lines = @{$_[0]};
	my $field = $_[1];
	my $delim = $_[2];
	my $sum = 0;
	my @f = ();
	foreach my $line ( @lines ) {
		@f = split($delim,$line);
		$sum += $f[$field];
	} 
	return $sum;
}

my %data = ();
my @tmp;
my ($pats,$reads,$type,$key);
foreach my $file ( @ARGV ) {
	my $filename = basename($file);
	open(IN,'<',$file) or die("Cannot open $file for reading.");
	if ( $filename =~ /^QC_log/ ) {
		$/ = "\n"; my @lines = <IN>; chomp(@lines);
		
		if ( defined($paired) ) {
			@tmp = grep( /LEFT\.log/,@lines );
			$data{'0-R1'} = [sum_field(\@tmp,1),'NA','NA'];
			@tmp = grep( /RIGHT\.log/,@lines );
			$data{'0-R2'} = [sum_field(\@tmp,1),'NA','NA'];
		}
		my $initial = sum_field(\@lines,1);	
		$data{'1-initial'} = [$initial,'NA','NA'];
		$data{'2-failQC'} = [$initial - sum_field(\@lines,2),'NA','NA'];
	} elsif ( $filename =~ /\.(nomatch|chim|fa2|fa)$/ ) {
		$type = $1; $/ = "\n"; my @lines = <IN>; chomp(@lines);
		@lines = grep(/^>/,@lines);
		$reads = scalar(@lines);
		$pats = sum(\@lines,1,'%');

		if ( $filename =~ /R0\.fa/ ) {
			$data{'2-passQC'} = [$reads,$pats,'NA'];
		} else {
			if ( $type eq 'chim' ) {
				$key = '3-chimeric';
			} elsif ( $type eq 'nomatch' ) {
				$key = '3-nomatch';
			} elsif ( $filename =~ /UNRECOGNIZABLE/ ) {
				$key = '3-unrecognizable';
			}

			if ( !defined($data{$key}) ) {
				$data{$key} = [$reads,$pats,'NA'];
			} else {
				$data{$key}[0] += $reads;
				$data{$key}[1] += $pats;
			}
		}
	}
	close(IN);
}

print STDOUT "Record\tReads\tPatterns\tPairsAndWidows\n";
foreach my $d ( sort(keys(%data)) ) {
	print STDERR $d,"\t",join("\t",@{$data{$d}}),"\n";
}
