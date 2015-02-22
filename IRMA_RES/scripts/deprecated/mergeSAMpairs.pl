#!/usr/bin/env perl
# Sam Shepard -- paired end merger
# 3.2014

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
	die($message."\n");
}

# FUNCTIONS #
sub condenseCigar($) {
	my $cig = $_[0];
	my $cigar = '';
	my $state = '';
	while( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+)/g ) {
		$state = $1;
		$cigar .= length($state);
		$cigar .= substr($state,0,1);
	}
	return $cigar;
}

sub lgg($) {
	return log($_[0])/log(10);
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
#############

$REisBase = qr/[ATCG]/;
$REgetMolID = qr/(.+?)[_ ]([12]):.+/;

%seqByRname = ();
$/ = ">";
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$rn = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	$seqByRname{$rn} = [split('',uc($seq))];
	$N = $lenByRname{$rn} = length($seq);
}
close(REF);

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
@sam = <SAM>; chomp(@sam);
close(SAM);

%pairs = %insByIndex = %tableByRname = ();
for($K=0;$K<scalar(@sam);$K++) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
	if ( $qname =~ $REgetMolID ) {
		$qMolID = $1;
		$qSide = $2;
	}

	if ( defined($seqByRname{$rn}) ) {
		$N = $lenByRname{$rn};
		@NTs = split('',uc($seq));
		@QCs = split('',$qual);
		@cigars = split('',$cigar);
		$rpos=$pos-1;
		$qpos=0;
		
		if ( $rpos > 0 ) {
			$aln = '.' x $rpos; 
			$qAln = ' ' x $rpos;
		} else {
			$aln = '';
			$qAln = '';
		}
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				for(1..$inc) {
					$qAln .= $QCs[$qpos];
					$aln .= $NTs[$qpos];
					$tableByRname{$rn}[$rpos]{$NTs[$qpos]}++;
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$qAln .= ' ' x $inc;
				$aln .= '-' x $inc;
				for(1..$inc) {
					$tableByRname{$rn}[$rpos]{'-'}++;
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insByIndex{$rn}{$K}{$rpos} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif( $op eq 'N' ) {
				$rpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}
		$aln .= '.' x (($N)-$rpos);
		$qAln .= ' ' x (($N)-$rpos);
		$pairs{$rn}{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
	}
}

%stats = ();
open(SAM,'>',$ARGV[2].'.sam') or die("Cannot open $ARGV[2].sam for writing.\n");
open(OUT,'>',$ARGV[2].'-merged.txt') or die("Cannot open $ARGV[2]-merged.txt\n");
print SAM $sam[0],"\n",$sam[1],"\n";
foreach $rn ( keys(%pairs) ) {
	$obs = $fmv = $hObs = $tmv = $TMJ = 0;
	@seqRef = @{$seqByRname{$rn}};
	foreach $mID ( keys(%{$pairs{$rn}}) ) {
		@mPairs = keys(%{$pairs{$rn}{$mID}});
		if ( scalar(@mPairs) == 2 ) {
			@a1 = @{$pairs{$rn}{$mID}{'1'}};			
			@a2 = @{$pairs{$rn}{$mID}{'2'}};
			($s1,$e1) = ($a1[3],$a1[4]);
			($s2,$e2) = ($a2[3],$a2[4]);
			
			$start = min($s1,$s2);
			$end = max($e1,$e2);
			$mSeq = '';
			$cigars = '';
			$qSeq = '';
			
			$K1 = $a1[2];
			$K2 = $a2[2];
			@bases1 = split('',$a1[0]);
			@bases2 = split('',$a2[0]);
			@quals1 = unpack("c* i*",$a1[1]);
			@quals2 = unpack("c* i*",$a2[1]);

			for($i=$start;$i<=$end;$i++) {
				$x = $bases1[$i];
				$y = $bases2[$i];
				$qx = $quals1[$i];
				$qy = $quals2[$i];
				$r = $seqRef[$i];

				if ( defined($insByIndex{$rn}{$K1}{$i}) && defined($insByIndex{$rn}{$K2}{$i}) ) {
					$ins1 = lc($insByIndex{$rn}{$K1}{$i}[0]);
					$ins2 = lc($insByIndex{$rn}{$K2}{$i}[0]);
					if ( $ins1 eq $ins2 ) {
						$mSeq .= $ins1;
						@qIns1 = split('',$insByIndex{$rn}{$K1}{$i}[1]);	
						@qIns2 = split('',$insByIndex{$rn}{$K2}{$i}[1]);
						for($qIndex=0;$qIndex<length($ins1);$qIndex++ ) {
							$qSeq .= chr(max(ord($qIns1[$qIndex]),ord($qIns2[$qIndex])));
						}
						$cigars .= 'I' x length($ins1);
					} elsif ( $ins2 =~ /$ins1/ ) {
						# 1 in 2
						$mSeq .= $ins1;
						$qSeq .= $insByIndex{$rn}{$K1}{$i}[1];
						$cigars .= 'I' x length($ins1);
					} elsif ( $ins1 =~ /$ins2/ ) {
						# 2 in 1
						$mSeq .= $ins2;
						$qSeq .= $insByIndex{$rn}{$K2}{$i}[1];
						$cigars .= 'I' x length($ins2);
					}
				}

				if ( $x ne '.' && $y ne '.' ) {
					$obs++;
					if ( $x eq $y ) {
						if ( $x eq '-' ) {
							$tmv++;
							$cigars .= 'D';
							$tableByRname{$rn}[$i]{'-'}--;
						} else {
							if ( $x ne $r ) { $tmv++; }
							$cigars .= 'M';
							$mSeq .= $x;	
							$qSeq .= chr(max($qx,$qy));
							$tableByRname{$rn}[$i]{$x}--;
						}
					} elsif( $x eq $r ) {
						$fmv++;
						$mSeq .= $x;
						$qSeq .= chr($qx);
						$cigars .= 'M';
						$tableByRname{$rn}[$i]{$y}--;
					} elsif( $y eq $r ) {
						$fmv++;
						$mSeq .= $y;
						$qSeq .= chr($qy);
						$cigars .= 'M';
						$tableByRname{$rn}[$i]{$x}--;
					} else {
						$fmv++;
						if ( $x =~ $REisBase && $y !~ $REisBase ) {
							$cigars .= 'M';
							$mSeq .= $x;
							$qSeq .= chr($qx);
							$tableByRname{$rn}[$i]{$y}--;
						} elsif ( $x !~ $REisBase && $y =~ $REisBase ) {
							$cigars .= 'M';
							$mSeq .= $y;
							$qSeq .= chr($qy);
							$tableByRname{$rn}[$i]{$x}--;
						} elsif ( $qx > ($qy+4) ) {
							$cigars .= 'M';
							$mSeq .= $x;
							$qSeq .= chr($qx);
							$tableByRname{$rn}[$i]{$y}--;
						} elsif ( $qy > ($qx+4) ) {
							$cigars .= 'M';
							$mSeq .= $y;
							$qSeq .= chr($qy);
							$tableByRname{$rn}[$i]{$x}--;
						} else {
							$cigars .= 'M';
							$mSeq .= 'N';
							$qSeq .= chr(int(avg($qx,$qy)));
							$tableByRname{$rn}[$i]{$x}--;
							$tableByRname{$rn}[$i]{$y}--;
							$tableByRname{$rn}[$i]{'N'}++;
						}
					}
				} elsif ( $x eq '.' && $y ne '.' ) {
					if ( $y eq '-' ) {
						$cigars .= 'D';
					} else {
						$cigars .= 'M';	
						$mSeq .= $y;
						$qSeq .= chr($qy);
					}
				} elsif( $x ne '.' && $y eq '.' ) {
					if ( $x eq '-' ) {
						$cigars .= 'D';
					} else {
						$cigars .= 'M';	
						$mSeq .= $x;
						$qSeq .= chr($qx);
					}
				} else {
					$cigars .= 'N';
				}
			}
			
			$qname = $a1[5];
			$mapq = avg($a1[6],$a2[6]);
			$qname =~ s/(.+?[_ ])[12](:.+)/${1}3${2}/;
			print SAM $qname,"\t",'0',"\t",$rn,"\t",($start+1),"\t",$mapq;
			print SAM "\t",condenseCigar($cigars),"\t*\t0\t0\t",$mSeq,"\t",$qSeq,"\n";
		} else {
			$K = $pairs{$rn}{$mID}{$mPairs[0]}[2];
			print SAM $sam[$K],"\n";
		}
	}

	$TMJ = $obs - $fmv - $tmv;
	$hObs = $TMJ + $tmv;
	print OUT $rn,"\t","TrueMajor","\t",$TMJ,"\n";
	print OUT $rn,"\t","FalseMinor","\t",$fmv,"\n";
	print OUT $rn,"\t","TrueMinor","\t",$tmv,"\n";
	print OUT $rn,"\t","HomogeneousObs","\t",$hObs,"\n";
	print OUT $rn,"\t","Observations","\t",$obs,"\n";
}
close(SAM);
close(OUT);

open(COVG,'>',$ARGV[2].'-coverage.txt') or die("Cannot open $ARGV[2]-coverage.txt for writing.\n");
open(CONS,'>',$ARGV[2].'.ref') or die("Cannot open $ARGV[2].ref for writing.\n");
print COVG "Position\tCoverage\n";
foreach $rn ( keys(%tableByRname) ) {
	$N = $lenByRname{$rn};
	print CONS '>',$rn,"\n";
	for($p=0;$p<$N;$p++) {
		$consensus = '.';
		$coverage = -1;
		$total = 0;
		foreach $base ( keys(%{$tableByRname{$rn}[$p]}) ) {
			if ( $tableByRname{$rn}[$p]{$base} > $coverage ) {
				$coverage = $tableByRname{$rn}[$p]{$base};
				$consensus = $base;
			}
			$total += $tableByRname{$rn}[$p]{$base};
		}
		
		print COVG ($p+1),"\t",$total,"\n";
		print CONS $consensus;
	}
	print CONS "\n";
}
close(CONS);
close(COVG);
