#!/usr/bin/env perl
# Sam Shepard -- phase.pl - 2.17.2016
# parallelized phasing algorithm

use Storable;
use constant NO_PHASING => 1;
use constant NO_PHASING_VEC => (1,1,1,1);
use constant FULLY_PHASED => 0;
use constant FULLY_PHASED_VEC => (0,0,0,0);

use Getopt::Long;
GetOptions(
		'array-size|S=i' => \$arraySize,
		'index|I=i' => \$index
	);

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <prefix> <var.sto> <pat.sto> [options]\n";
	$message .= "\t\t-S|--array-size\t\tArray size for array job (default 10).\n";
	die($message."\n");
}

$prefix = $ARGV[0]; $varFile = $ARGV[1]; $patFile = $ARGV[2];
if( !defined($arraySize) ) {
	$arraySize = 10;
} elsif ( $arraySize < 1  ) {
	die("Array size must be 1 or greater.\n");
}


if ( defined($index ) ) {
	if ( $index < 1 || $index > $arraySize ) {
		die("Index must be 1 to $arraySize.\n");
	} else {
		$index = int($index);
		$prefix = $prefix .'-'. sprintf("%04d",$index);
	}
} else {
	$index = 1;
	$arraySize = 1;
}


sub max($$) {
	if ( $_[0] > $_[1] ) {
		return $_[0];
	} else {
		return $_[1];
	}
}

sub min($$) {
	if ( $_[0] < $_[1] ) {
		return $_[0];
	} else {
		return $_[1];
	}
}

sub avg($$) {
	return (($_[0]+$_[1])/2);
}

%readPats = ();
%variants = %{retrieve($varFile)};
%readPats = %{retrieve($patFile)};

# vectorize
$N = 0; @alleles = @names = ();
@sites = sort { $a <=> $b } keys(%variants);
foreach $index ( 0 .. $#sites ) {
	$v = $sites[$index];
	foreach $b ( sort { $a cmp $b } keys(%{$variants{$v}}) ) {
		$alleles[$N] = [$index,$b,$variants{$v}{$b}];
		$names[$N] = ($v+1).$b;
		$N++;
	}
}

$O = ($N**2 - $N)/2;
if ( $arraySize > $O ) {
	print STDERR "WARNING: array size ($arraySize) greater than operations ($O)!\n";
	print STDERR "Setting array size to $O\n";
	$arraySize = $O;
	if ( $index > $arraySize ) {
		die("WARNING: index ($index) greater than adjusted array size($arraySize). Aborting.\n");
	}
}

open(EXP,'>',$prefix."-EXPENRD.sqm") or die("Cannot open ${prefix}_EXPENRD.sqm for writing.\n");
open(JAC,'>',$prefix."-JACCARD.sqm") or die("Cannot open ${prefix}_JACCARD.sqm for writing.\n");
open(MUT,'>',$prefix."-MUTUALD.sqm") or die("Cannot open ${prefix}_MUTUALD.sqm for writing.\n");
open(JOP,'>',$prefix."-NJOINTP.sqm") or die("Cannot open ${prefix}_NJOINTP.sqm for writing.\n");

$tN = $N;		  	# totalN left
$Qb = int($O/$arraySize); 	# base Quantity to do
$s = 1 + $Qb*($index - 1);	# starting operation 				
if ( $index == $arraySize ) {
	$Qb += $O % $arraySize;	# for last operation clean up remainder
} 

if( $index == 1 ) {
	printToMatrix(FULLY_PHASED_VEC,$names[0]."\t","\n");
}

# Derive initial coordinates using triangle series:
# 	http://www.mathsisfun.com/algebra/triangular-numbers.html
$row = int((1+sqrt(1+8*($s-1)))/2);	# starting row coord
$fr = ($row*($row-1))/2+1;		# first in row
$col = $s % $fr;			# starting col coord

if ( $col == 0 ) {
	printToMatrix('','','','',$names[$row],'');
}

while( $col < $row ) {
	printToMatrix(dist($row,$col),"\t",'');
	$Qb--; $col++;
	if ( $Qb == 0 ) {
		last;
	}
}
if ( $row == $col ) {
	printToMatrix(FULLY_PHASED_VEC,"\t","\n");
}
$row++;

while( $Qb > 0 ) {
	$col = 0;
	printToMatrix('','','','',$names[$row],'');
	while( $col<$row ) {
		printToMatrix(dist($row,$col),"\t",'');
		$Qb--; $col++;
		if ( $Qb == 0 ) {
			last;
		}
	}
	if ( $row == $col ) {
		printToMatrix(FULLY_PHASED_VEC,"\t","\n");
	}
	$row++;
}


sub printToMatrix($$$$$$) {
	my ($mut,$jac,$exp,$jop,$prefix,$suffix) = @_;

	print EXP $prefix,$exp,$suffix;
	print JAC $prefix,$jac,$suffix;
	print MUT $prefix,$mut,$suffix;
	print JOP $prefix,$jop,$suffix;
}

sub dist($$) {
	my $i = $_[0];
	my $j = $_[1];
	my ($s1,$b1,$Fb1) = @{$alleles[$i]};
	my ($s2,$b2,$Fb2) = @{$alleles[$j]};
	my ($total,$Eb1,$Eb2,$Fb1b2,$X) = (0,0,0,0,0);
	my ($p1,$p2,$pat) = ('','','');
	my ($mn1,$mn2,$mx1,$mx2,$mnA) = (0,0,0,0,0);
	my ($mutd,$jacc,$expd,$njop) = (0,0,0,0);

	# different site
	if ( $s1 != $s2 ) {
		if ( $Fb1 == 0 || $Fb2 == 0 ) {
			return (NO_PHASING_VEC);
		}

		foreach $pat ( keys(%readPats) ) {
			$p1 = substr($pat,$s1,1);
			$p2 = substr($pat,$s2,1);
			if ( $p1 eq '.' || $p2 eq '.' ) { next; }
			
			$X = $readPats{$pat};
			$total += $X;
			if ( $p1 eq $b1 && $p2 eq $b2 ) { 
				$Fb1b2 += $X;
				$Eb1 += $X;
				$Eb2 += $X;	
			} elsif ( $p1 eq $b1 ) {
				$Eb1 += $X;
			} elsif ( $p2 eq $b2 ) {
				$Eb2 += $X;
			}
		}

		if ( $total != 0 ) {
			$Fb1b2 /= $total;
			$Eb1 /= $total;
			$Eb2 /= $total;

			if ( $Fb1b2 == 0 || $Eb1 == 0 || $Eb2 == 0 ) {
				return (NO_PHASING_VEC);
			}
		} else {
			return (NO_PHASING_VEC);
		}

		$mn1 = min($Eb1,$Fb1); $mn2 = min($Eb2,$Fb2); $mnA = min($mn2,$mn1);
		$mx1 = max($Eb1,$Fb1); $mx2 = max($Eb2,$Fb2);
		$mutd = 1 - $Fb1b2**2/($mx1*$mx2);
		$jacc = 1 - $Fb1b2/($mx1+$mx2-$Fb1b2);
		if ( $total <= 20 || $mutd == 0 || $jacc == 0 ) { 
			$expd = 1 - (($Fb1b2*$mnA)/($mx1*$mx2));
		} else {
			$expd = $jacc;
		}
		$njop = 1 - 2*$Fb1b2;
		return ($mutd,$jacc,$expd,$njop);

	# same site, and, same allele
	} elsif ( $b1 eq $b2 ) {
		return (FULLY_PHASED_VEC);
	# same site, and, different allele
	} else {
		return (NO_PHASING_VEC);
	}
}

close(EXP); close(JAC);
close(MUT); close(JOP);
