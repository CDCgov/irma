#!/usr/bin/env perl
# Sam Shepard - 2014-7-11
# variant calling effects on translation



if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <reference> <alleles.txt> [OPTIONS]";

	die($message."\n");
}

%gc = (
	'TCA' => 'S', # Serine
	'TCC' => 'S', # Serine
	'TCG' => 'S', # Serine
	'TCT' => 'S', # Serine
	'TTC' => 'F', # Phenylalanine
	'TTT' => 'F', # Phenylalanine
	'TTA' => 'L', # Leucine
	'TTG' => 'L', # Leucine
	'TAC' => 'Y', # Tyrosine
	'TAT' => 'Y', # Tyrosine
	'TAA' => '*', # Stop
	'TAG' => '*', # Stop
	'TGC' => 'C', # Cysteine
	'TGT' => 'C', # Cysteine
	'TGA' => '*', # Stop
	'TGG' => 'W', # Tryptophan
	'CTA' => 'L', # Leucine
	'CTC' => 'L', # Leucine
	'CTG' => 'L', # Leucine
	'CTT' => 'L', # Leucine
	'CCA' => 'P', # Proline
	'CAT' => 'H', # Histidine
	'CAA' => 'Q', # Glutamine
	'CAG' => 'Q', # Glutamine
	'CGA' => 'R', # Arginine
	'CGC' => 'R', # Arginine
	'CGG' => 'R', # Arginine
	'CGT' => 'R', # Arginine
	'ATA' => 'I', # Isoleucine
	'ATC' => 'I', # Isoleucine
	'ATT' => 'I', # Isoleucine
	'ATG' => 'M', # Methionine
	'ACA' => 'T', # Threonine
	'ACC' => 'T', # Threonine
	'ACG' => 'T', # Threonine
	'ACT' => 'T', # Threonine
	'AAC' => 'N', # Asparagine
	'AAT' => 'N', # Asparagine
	'AAA' => 'K', # Lysine
	'AAG' => 'K', # Lysine
	'AGC' => 'S', # Serine
	'AGT' => 'S', # Serine
	'AGA' => 'R', # Arginine
	'AGG' => 'R', # Arginine
	'CCC' => 'P', # Proline
	'CCG' => 'P', # Proline
	'CCT' => 'P', # Proline
	'CAC' => 'H', # Histidine
	'GTA' => 'V', # Valine
	'GTC' => 'V', # Valine
	'GTG' => 'V', # Valine
	'GTT' => 'V', # Valine
	'GCA' => 'A', # Alanine
	'GCC' => 'A', # Alanine
	'GCG' => 'A', # Alanine
	'GCT' => 'A', # Alanine
	'GAC' => 'D', # Aspartic Acid
	'GAT' => 'D', # Aspartic Acid
	'GAA' => 'E', # Glutamic Acid
	'GAG' => 'E', # Glutamic Acid
	'GGA' => 'G', # Glycine
	'GGC' => 'G', # Glycine
	'GGG' => 'G', # Glycine
	'GGT' => 'G'  # Glycine
);

# physio-chemical factors
# Atchley et al. 2008
# "Solving the protein sequence metric problem."
# Proc Natl Acad Sci U S A. 2005 May 3;102(18):6395-400. Epub 2005 Apr 25.
%pcf = (
	'A' => [-0.59,-1.3,-0.73,1.57,-0.15],
	'C' => [-1.34,0.47,-0.86,-1.02,-0.26],
	'D' => [1.05,0.3,-3.66,-0.26,-3.24],
	'E' => [1.36,-1.45,1.48,0.11,-0.84],
	'F' => [-1.01,-0.59,1.89,-0.4,0.41],
	'G' => [-0.38,1.65,1.33,1.05,2.06],
	'H' => [0.34,-0.42,-1.67,-1.47,-0.08],
	'I' => [-1.24,-0.55,2.13,0.39,0.82],
	'K' => [1.83,-0.56,0.53,-0.28,1.65],
	'L' => [-1.02,-0.99,-1.51,1.27,-0.91],
	'M' => [-0.66,-1.52,2.22,-1.01,1.21],
	'N' => [0.95,0.83,1.3,-0.17,0.93],
	'P' => [0.19,2.08,-1.63,0.42,-1.39],
	'Q' => [0.93,-0.18,-3.01,-0.5,-1.85],
	'R' => [1.54,-0.06,1.5,0.44,2.9],
	'S' => [-0.23,1.4,-4.76,0.67,-2.65],
	'T' => [-0.03,0.33,2.21,0.91,1.31],
	'V' => [-1.34,-0.28,-0.54,1.24,-1.26],
	'W' => [-0.6,0.01,0.67,-2.13,-0.18],
	'Y' => [0.26,0.83,3.1,-0.84,1.51],
	'-' => [0,0,0,0,0]
	);

@aa =sort(keys(%pcf));
$N = scalar(@aa);
%distMat = ();
for($i = 0; $i < $N; $i++ ) {
	for($j = 0; $j < $N; $j++ ) {
		$dist = 0;
		for( $k = 0; $k < 5; $k++ ) {
			$dist += ($pcf{$aa[$i]}[$k] - $pcf{$aa[$j]}[$k]) **2;
		}
		$dist = sqrt($dist);
		$distMat{$aa[$i]}{$aa[$j]} = $dist;
	}
}
$maxDist = $distMat{'S'}{'Y'};
shift(@aa); # remove gaps
open(IN, '<', $ARGV[0] ) or die("Cannot open $ARGV[0] for reading.\n");
$/ = ">";
$record = <IN>;
$record = <IN>;
chomp($record);
close(IN);
@lines = split(/\r\n|\n|\r/, $record);
$header = shift(@lines);
$sequence = uc(join('',@lines));
$length = length($sequence);

$aaSeq[0] = $aaSeq[1] = $aaSeq[2] = '';
@seq = split('',$sequence);
for($i=0;$i < $length;$i+=3 ) {
	if ( ($i+0+2) < $length ) {
		$codon1 = substr($sequence,$i,3);
		$aaSeq[0] .= $gc{$codon1};
	}

	if ( ($i+1+2) < $length ) {
		$codon2 = substr($sequence,$i+1,3);
		$aaSeq[1] .= $gc{$codon2};
	}

	if ( ($i+2+2) < $length ) {
		$codon3 = substr($sequence,$i+2,3);
		$aaSeq[2] .= $gc{$codon3};
	}
}

$sc1 = ($aaSeq[0] =~ tr/*//);
$sc2 = ($aaSeq[1] =~ tr/*//);
$sc3 = ($aaSeq[2] =~ tr/*//);
if ( $sc1 <= $sc2 && $sc1 <= $sc3 ) {
	$aaRef = $aaSeq[0];
	$offset = 0;
	print STDERR "Frame 1 chosen.\n";
} elsif ( $sc2 <= $sc1 && $sc2 <= $sc3 ) {
	$aaRef = $aaSeq[1];
	$offset = 1;
	print STDERR "Frame 2 chosen.\n";
} else {
	$aaRef = $aaSeq[2];
	$offset = 2;
	print STDERR "Frame 3 chosen.\n";
}
$N = $length - $offset - ($length%3);
$UB = $N - 1;
@ref = split('',$aaRef);

print STDERR $UB,"\t",$offset,"\n";
open(IN,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
$/ = "\n";
$header = <IN>; chomp($header);
print $header,"\thasSubstitution\tMutation\tPCD\taaAllele\tSignificant\tNonsense\n";
while($record=<IN>) {
	chomp($record);
	@fields	= split("\t",$record);
	($position,$allele,$type) = ($fields[1],$fields[2],$fields[10]);
	($frequency,$pairedUB,$qualityUB) = ($fields[5],$fields[8],$fields[9]);
	$position--;
	$CN = int( ($position-$offset) / 3);
	$aaAllele = $ref[$CN];

	if ( $frequency > $pairedUB && $frequency > $qualityUB ) {
		$sig = 'Yes';
	} else {
		$sig = 'No';
	}

	if ( $type eq 'Majority' ) {
		$pcd = 0;
		$sub = 'No';
		$mutation = $ref[$CN];
		$nonsense = 'No';
	} elsif ( $position < $offset || $position > $UB ) {
		$pcd = 0;
		$sub = 'OOB';
		$mutation = '';
		$nonsense = 'No';
	} else {
		$CP = ($position-$offset)%3;
		$codon = substr($sequence,$position-$CP,3);
		substr($codon,$CP,1,$allele);
		$aa = $gc{$codon};
		if ( $aa eq $ref[$CN] ) {
			$pcd = 0;
			$sub = 'No';
			$mutation = '';
		} else {
			$pcd = $distMat{$aa}{$ref[$CN]};
			$sub = 'Yes';
			$mutation = $ref[$CN].($position+1).$aa;
		}
		$aaAllele = $aa;

		if ( $aa eq '*' ) {
			$nonsense = 'Yes';
		} else {
			$nonsense = 'No';
		}
	}
	
	print $record,"\t",$sub,"\t",$mutation,"\t",$pcd,"\t",$aaAllele,"\t",$sig,"\t",$nonsense,"\n";
}
close(IN);
