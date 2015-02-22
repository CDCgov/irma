#!/usr/bin/env perl
# Sam Shepard - 9.2014
# writePairingStats

use Storable;
#use Getopt::Long;
#GetOptions(	'no-gap-allele|G' => \$noGap );

if ( scalar(@ARGV) < 1 ) {
	$message = "Usage:\n\tperl $0 <stats1> <...>\n";
	die($message."\n");
}

@keys = ('dmv','fmv','tmv','obs');
%table = ();

foreach $file ( @ARGV ) {
	%data = %{retrieve($file)};

	foreach $rn ( keys(%data) ) {
		foreach $key ( @keys ) {
			$table{$rn}{$key} += $data{$rn}{$key};
		}
	}
}

foreach $rn ( keys(%table) ) {
	$obs = $table{$rn}{'obs'};
	$tmv = $table{$rn}{'tmv'};
	$fmv = $table{$rn}{'fmv'};
	$dmv = $table{$rn}{'dmv'};
	$insObs = $table{$rn}{'insObs'};
	$insErr = $table{$rn}{'insErr'};

	$TMJ = $obs - $fmv - $tmv;
	$hObs = $TMJ + $tmv;

	if ( $obs > 0 ) {
		print $rn,"\tObservations\t",$obs,"\n";
		print $rn,"\tExpectedErrorRate\t",($fmv/$obs),"\n";
		print $rn,"\tMinimumExpectedVariation\t",($tmv/$obs),"\n";
		print $rn,"\tMinimumDeletionErrorRate\t",($dmv/$obs),"\n";
	} else {
		print $rn,"\tObservations\t",$obs,"\n";
		print $rn,"\tExpectedErrorRate\t0\n";
		print $rn,"\tMinimumExpectedVariation\t0\n";
		print $rn,"\tMinimumDeletionErrorRate\t0\n";
	}

	if ( $insObs > 0 ) {
		if ( $insObs == $insErr ) {
			$insObs++;
		}
		$insErr /= $insObs;
		$IEbyRef{$rn} = $insErr;
		print $rn,"\tMinimumInsertionErrorRate\t$insErr\n";
	} else {
		print $rn,"\tMinimumInsertionErrorRate\t0\n";
	}

}
