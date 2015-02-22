#!/usr/bin/env perl

# Reference_Name // Position //Allele // Count // Total // Frequency // Average_Quality // ConfidenceNotMacErr


if ( scalar(@ARGV) != 2 ) {
	$message = "Usage:\n\tperl $0 <variants.txt> <pairing-stats.txt>";
	die($message."\n");
}

### NEGATIVE BINOMIAL ###
$b=1;
#$kappa = 1.644854;
#$kappa = 2.326348;
$kappa = 3.090232;
$kappa2 = $kappa ** 2;
$eta = $kappa2/3 + 1/6;
$gamma1 = $b*( $kappa2*(13/18) + 17/18);
$gamma2 = $kappa2*(1/18)+7/36;
#########################

open(IN2,'<',$ARGV[1]) or die("Cannot open $ARGV[1].\n");
while($line=<IN>) {
	chomp($line);
	($ref,$type,$value) = split("\t",$line);
	if ( $type eq 'ExpectedErrorRate' ) {
		$EEpaired = $value;
	}
}
close(IN2);

open(IN,'<',$ARGV[0]) or die("Cannot open $ARGV[0].\n");
$line =<IN>;
while($line=<IN>) {
	chomp($line);
	@fields = split("\t",$line);
	$F = $fields[5];
	$N = $fields[4];
	$Q = $fields[6];
	$X = $fields[3];

	$EE = 1/(10**($Q/10));
	
	$UB = UB($EE,$N);
	$UBpaired = UB($EEpaired,$N);

	print $fields[1],$fields[2],"\t",$fields[4];	
	if ( $F > $UB && $F > $UBpaired ) {
		print sprintf("\t%.8f in [0,%.8f or %.8f]\t%.8f+*\n",$EE,$UBpaired,$UB,$F);
	} elsif ( $F > $UB ) {
		print sprintf("\t%.8f in [0,%.8f or %.8f]\t%.8f*\n",$EE,$UBpaired,$UB,$F);
	} elsif ( $F > $UBpaired ) {
		print sprintf("\t%.8f in [0,%.8f or %.8f]\t%.8f+\n",$EE,$UBpaired,$UB,$F);
	} else {
		print sprintf("\t%.8f in [0,%.8f or %.8f]\t%.8f\n",$EE,$UBpaired,$UB,$F);
	}
}

sub UB($$) {
	my $p = $_[0];
	my $N = $_[1];
	my $u= $p/(1-$p);
	my $V= $u + $b*($u**2);
	my $u2= ($p*$N + $eta)/($N+ $b*(-2)*$eta);
	my $UB = $u2 + $kappa*sqrt($V+($gamma1*$V+$gamma2)/$N)*(1/sqrt($N));
	return $UB;
}
