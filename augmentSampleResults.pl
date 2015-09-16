#!/usr/bin/env perl
# Sam Shepard

use File::Basename;
use Getopt::Long;
GetOptions(	'word|W=s'=> \$word, 'suffix|S=s' => \$suffix, 'ID-only|I' => \$idOnly, 'infix|X=s' => \$infix,
		'dir-field|F=i' =>  \$field, 'ignore-dir|G' => \$ignoreIDfield, "no-header|H" => \$noHeader
 );

if ( scalar(@ARGV) != 1 ) {
	$message = "Usage:\n\tperl $0 <key.txt> [options]\n";
	$message .= "\t\t-H|--no-header\t\tDo not output header.\n";
	$message .= "\t\t-W|--word <STR>]\tContains <word> in table filename.\n";
	$message .= "\t\t-S|--suffix <STR>\tTable filename ends in suffix.\n";
	$message .= "\t\t-I|--ID-only\t\tKey just contains the ID.\n";
	$message .= "\t\t-X|--infix <STR>\tInfix is the intermediate path. Default is 'tables/'\n";
	$message .= "\t\t-F|--dir-field <INT>\tField containing dirname.\n"; 
	$message .= "\t\t-G|--ignore-dir\t\tDo not print out dirname into collated data.\n";
	die($message."\n");
}

if ( !defined($word) ) {
	$word='all';
}

if ( !defined($infix) ) {
	$infix = 'tables';
}

if ( defined($field) && $field > 0 ) {
	$sampleDirField = $field;
} else {
	$sampleDirField = 1;
}
$first = 1;
open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].\n");
$hdrs =<IN>; chomp($hdrs);
@tmp = split("\t",$hdrs);
if ( defined($idOnly) ) {
	$hdrs = $tmp[0];
}
if ( $ignoreIDfield ) {
	$line = ''; $first = 1;
	for($i = 0;$i<scalar(@tmp);$i++ ) {
		if ( $i != ($sampleDirField-1) ) {
			if ( $first ) {
				$first = 0;
				$line = $tmp[$i];
			} else {
				$line .= "\t".$tmp[$i];
			}
		}
	}
	$hdrs=$line;
}

$firstData = 1;
while($line=<IN>) {
	chomp($line);
	@fields = split("\t",$line);
	$id = $fields[$sampleDirField-1];
	@files=glob("$id/$infix/*$word*$suffix");
	if ( defined($idOnly) ) {
		$line = $fields[0];
	} 

	if ( $ignoreIDfield ) {
		$line = ''; $first = 1;
		for($i = 0;$i<scalar(@fields);$i++ ) {
			if ( $i != ($sampleDirField-1) ) {
				if ( $first ) {
					$first = 0;
					$line = $fields[$i];
				} else {
					$line .= "\t". $fields[$i];
				}
			}
		}
	}
	
	foreach $file ( @files ) {
		open(DAT,'<',$file) or die("Cannot open $file for reading.\n");
		$header = <DAT>; chomp($header);
		$numCols = scalar(split("\t",$header));
		$basename = basename($file);

		if ( $firstData ) {
			if ( ! $noHeader ) {
				print "$hdrs\t",$header,"\n";
			}
			$firstData = 0;
		}

		while($line2=<DAT>) {
			chomp($line2);
			print $line,"\t",$line2,"\n";
		}

		close(DAT);
	}
}
close(IN);
