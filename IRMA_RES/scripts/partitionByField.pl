#!/usr/bin/env perl
# partitionByField -  Version 1.0
# Partitions sequences by fields in the fasta header given a delimiter.
#
# Copyright (C) 2012, Centers for Disease Control & Prevention
# Author: Samuel S. Shepard (vfn4@cdc.gov)
#
# GPL version 3
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

####################
# PROCESS parameters
use Getopt::Long;
GetOptions(
		'max-groups|M=i'=> \$maxGroups,
		'delimiter|D:s'=> \$delim,
		'field|F:i' => \$field,
		'strip-file-braces|S' => \$stripBraces,
		'suffix|X=s' => \$suffix,
		'use-header|H' => \$byHeader,
		'ignore-annotation|A' => \$ignoreAnnotation
		);

if ( scalar(@ARGV) < 2 ) {
	$message = "\n$0 <input.fasta> <output_prefix> [OPTIONS]\n";
	$message .= "\t-S|--strip-file-braces\t\t\tStrip braces from annotations on filenames.\n";
	$message .= "\t-X|--suffix STRING\t\t\tExtension to add to the file. Default = source file.\n";
	$message .= "\t-M|--max-groups POSITIVE_NUMBER\t\tMaximum number of groups allowed..\n";
	$message .= "\t-D|--delimiter CHARACTER\t\tSingle character delimiter within FASTA header. Default is '|'.\n";
	$message .= "\t-F|--field POSITIVE_NUMBER\t\tHeader field to use for pairing records. Default is 1.\n";
	$message .= "\t-A|--ignore-annotation\t\t\tIgnore annotation when getting field values. Think {annot}.\n";
	$message .= "\t-H|--use-header\t\t\t\tUse header as the partition field. (Dangerous)\n";
	die($message."\n");
}

if ( defined($maxGroups) && $maxGroups <= 0 ) {
	die("$0 ERROR: The number of groups must be a POSITIVE integer.\n");
} elsif ( !defined($maxGroups) ) {
		$maxGroups = 100;
}

if ( !defined($field) ) {
	$field = 0;
} elsif($field == 0 ) {
	die("$0 ERROR: field must be specified.\n");
} elsif($field < 0) {
	die("$0 ERROR: field must be a positive number.\n");
} else {
	$field -= 1;
}

if ( !defined($delim) ) {
	$delim = '|';
} elsif( $delim eq '' ) {
	die("$0 ERROR: No delimiter argument detected.\n");
} elsif( length($delim) > 1 ) {
	die("$0 ERROR: single character delimiter expected instead of '$delim'.\n");
}

if ( !defined($suffix) ) {
	@parts = split(/\./, $ARGV[0]);
	if ( scalar(@parts) == 1 ) {
		$suffix = '';
	} else {
		$suffix = '.'.$parts[$#parts];
	}
} else {
	$suffix =~ tr/.//d;
	$suffix = '.'.$suffix;
}
$prefix = $ARGV[1];

# PROCESS fasta data
$separator = '_'; if ( defined($byHeader) ) { $separator = ''; }
open( IN, '<', $ARGV[0] ) or die("$0 ERROR: Could not open $ARGV[0].\n");
@headers = (); @sequences = (); %fileHandles = (); %identifiers = ();
$/ = ">"; $i = 0; $numberGroups = 0;
while( $record = <IN> ) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$header = shift(@lines);
	$sequence = lc(join('',@lines));
	
	@fields = split( /\Q$delim\E/, $header );
	$length = length($sequence);

	if ( $length == 0 ) {
		next;	
	} elsif ( $header !~ /\Q$delim\E/ && $delim ne '' && !defined($byHeader) ) {
		die("$0 ERROR: delimiter '$delim' not found in header:\n\t$header\n");
	} elsif ( $#fields < $field && 0 ) {
		$x = $#fields + 1; $field++;
		die("$0 ERROR: Only $x fields found in header but field $field selected:\n\t$header\n");
	}	
	$i++;

	if ( $byHeader ) {
		$key = $header;
	} else {
		$key = $fields[$field];
	}

	if ( $key ne '' ) {
		$key =~ tr/\//-/;
		$key =~ tr/+//d;
		$key =~ /^\s*(.*?)\s*$/;
		$key = $1;
		$key =~ tr/ /_/;
	} else {
		$key = 'BLANK';
	}

	if ( $ignoreAnnotation ) {
		$key =~ s/{.+?}//;
	}

	if ( !exists($identifiers{$key}) ) {
		$numberGroups++;
		if ( $numberGroups > $maxGroups ) {
			die("$0 ERROR: Number of groups ($numberGroups) exceeds the maximum permitted ($maxGroups). Use -M option to adjust this value.\n");
		}

		$group = $key;
		if ( $stripBraces ) {
			$group =~ tr/{}//d;
		}

		open($fileHandles{$key},'>', $prefix.$separator.$group.$suffix) or die("$0 ERROR: Cannot open $prefix.'_'.$group.$suffix\n");
		$handle = $fileHandles{$key};
		print $handle '>',$header,"\n",$sequence,"\n";
		$identifiers{$key} = 1;
	} else {
		$identifiers{$key}++;
		$handle = $fileHandles{$key};
		print $handle '>',$header,"\n",$sequence,"\n";
	}
}
close(IN);

foreach $handle ( keys( %fileHandles ) ) {
	close($handle);
}
####################
