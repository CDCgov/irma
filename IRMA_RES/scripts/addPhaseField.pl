#!/usr/bin/env perl
# Sam Shepard - 2017-06-06
# Add Phase Field - add to variant table the computed phase field

if ( scalar(@ARGV) != 2 ) {
	die("Usage:\n\t$0 <variants.txt> <phases.txt>\n");
}

open(VAR,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
@variants = <VAR>; chomp(@variants);
close(VAR);

if ( scalar(@variants) < 3 ) {
	open(VAR,'>',$ARGV[0]) or die("Cannot open $ARGV[0] for writing.\n");
	print VAR $variants[0],"\tPhase\n";
	if ( scalar(@variants) == 2 ) { print VAR $variants[1],"\t1\n"; }
	close(VAR);
} else {
	%phaseCount = %phaseByAllele = %rankByPhase = (); $i = 0;
	open(PHAZ,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
	while($line=<PHAZ>) {
		chomp($line);
		@t = split("\t",$line);
		($posAllele,$phase) = @t;
		$phaseByAllele{$posAllele} = $phase;
		$phaseCount{$phase}++;
		$i++;
	}
	close(PHAZ);
	if ( $i == 0 ) { die("No variant phases detected for file: $ARGV[0]\n"); }

	@phaseRank = sort { $phaseCount{$b} <=> $phaseCount{$a} || $a <=> $b } keys(%phaseCount);
	foreach $i ( 0 .. $#phaseRank ) { $rankByPhase{$phaseRank[$i]} = ($i+1); }	

	$headerLine = shift(@variants);
	@headers = split("\t",$headerLine);
	$pIndex = $aIndex = -1;
	foreach $i ( 0 .. $#headers ) {
		if ( $headers[$i] =~ /position/i ) { $pIndex = $i; }
		if ( $headers[$i] =~ /minority.allele/i ) { $aIndex = $i; }
	}
	if ( $aIndex == -1 || $pIndex == -1 ) { die("No position or minority allele field found in the variant table header.\n"); }

	open(VAR,'>',$ARGV[0]) or die("Cannot open $ARGV[0] for writing.\n");
	print VAR $headerLine,"\tPhase\n";
	foreach $varLine ( @variants ) {
		@t = split("\t",$varLine);
		
		if ( $#t < $#headers ) { die("Variant table line is shorter than the header.\n"); }
		$posAllele = $t[$pIndex] . $t[$aIndex];
		if ( defined($phaseByAllele{$posAllele}) ) {
			print VAR $varLine,"\t",$rankByPhase{$phaseByAllele{$posAllele}},"\n";
		} else {
			print VAR $varLine,"\tNULL\n";
		}
	}
	close(VAR);
}
