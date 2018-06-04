#!/usr/bin/env perl
# Sam Shepard -- mergeSAMpairs -- version 2
# 9.2014

use Storable;
use Getopt::Long;
GetOptions(
		'use-storable|S' => \$useStorable, 'bowtie-format|B' => \$bowtieFormat
	);

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
	$message .= "\t\t-S|--use-storable\tStore statistics object..\n";
	die($message."\n");
}


# FUNCTIONS #
sub condenseCigar($) {
	my $cig = $_[0];
	my $cigar = '';
	my $state = '';
	while( $cig =~ /([M]+|[D]+|[I]+|[H]+|[N]+|[S]+)/g ) {
		$state = $1;
		$cigar .= length($state);
		$cigar .= substr($state,0,1);
	}
	return $cigar;
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

$/ = ">"; 
$REF_LEN = $REF_NAME = '';
@REF_SEQ = ();	
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	@REF_SEQ = split('',uc($seq));
	$REF_LEN = length($seq);
	last;
}
close(REF);

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
@sam = <SAM>; chomp(@sam);
close(SAM);

%pairs = %insByIndex = ();
for($K=0;$K<scalar(@sam);$K++) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);
	if ( $bowtieFormat ) {
		$qMolID = $qname;
		if ( defined($pairs{$qMolID}) ) {
			$qSide = 2;
		} else {
			$qSide = 1;
		}
	} elsif ( $qname =~ $REgetMolID ) {
		$qMolID = $1;
		$qSide = $2;
	}

	if ( $REF_NAME eq $rn ) {
		@NTs = split('',uc($seq));
		@QCs = split('',$qual);
		@Qint = unpack("c* i*",$qual);
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
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$qAln .= ' ' x $inc;
				$aln .= '-' x $inc;
				for(1..$inc) {
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insByIndex{$K}{$rpos-1} = [substr($seq,$qpos,$inc),substr($qual,$qpos,$inc)];
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif( $op eq 'N' ) {
				$aln .= 'N' x $inc;
				$qAln .= ' ' x $inc;
				$rpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}
		$aln .= '.' x (($REF_LEN)-$rpos);
		$qAln .= ' ' x (($REF_LEN)-$rpos);
		$pairs{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
	}
}

%observations = ();
open(SAM,'>',$ARGV[2].'.sam') or die("Cannot open $ARGV[2].sam\n");
foreach $line ( @sam ) {
	if ( substr($line,0,1) eq '@' ) {
		print SAM $line,"\n";
	} else {
		last;
	}
}
$dmv = $obs = $fmv = $tmv = 0;
$insObs = $insErr = 0;
foreach $mID ( keys(%pairs) ) {
	@mPairs = keys(%{$pairs{$mID}});
	if ( scalar(@mPairs) == 2 ) {
		@a1 = @{$pairs{$mID}{'1'}};			
		@a2 = @{$pairs{$mID}{'2'}};
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
			$r = $REF_SEQ[$i];

			if ( $x ne '.' && $y ne '.' ) {
				$obs++;
				if ( $x eq $y ) {
					if ( $x eq '-' ) {
						$tmv++;
						$cigars .= 'D';
					} else {
						if ( $x ne $r ) { $tmv++; }
						$cigars .= 'M';
						$mSeq .= $x;	
						$qSeq .= chr(max($qx,$qy));
					}
				} elsif( $x eq $r ) {
					$fmv++;
					if ( $y eq '-' ) { $dmv++; }
					$mSeq .= $x;
					$qSeq .= chr($qx);
					$cigars .= 'M';
				} elsif( $y eq $r ) {
					$fmv++;
					if ( $x eq '-' ) { $dmv++; }
					$mSeq .= $y;
					$qSeq .= chr($qy);
					$cigars .= 'M';
				} else {
					$fmv++;
					if ( $x =~ $REisBase && $y !~ $REisBase ) {
						$cigars .= 'M';
						$mSeq .= $x;
						$qSeq .= chr($qx);
						if ( $y eq '-' ) { $dmv++; }
					} elsif ( $x !~ $REisBase && $y =~ $REisBase ) {
						$cigars .= 'M';
						$mSeq .= $y;
						$qSeq .= chr($qy);
						if ( $x eq '-' ) { $dmv++; }
					} elsif ( $qx > ($qy+4) ) {
						$cigars .= 'M';
						$mSeq .= $x;
						$qSeq .= chr($qx);
						if ( $y eq '-' ) { $dmv++; }
					} elsif ( $qy > ($qx+4) ) {
						$cigars .= 'M';
						$mSeq .= $y;
						$qSeq .= chr($qy);
						if ( $x eq '-' ) { $dmv++; }
					} else {
						$cigars .= 'M';
						$mSeq .= 'N';
						$qSeq .= chr(int(avg($qx,$qy)));
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


			if ( defined($insByIndex{$K1}{$i}) && defined($insByIndex{$K2}{$i}) ) {
				$ins1 = lc($insByIndex{$K1}{$i}[0]);
				$ins2 = lc($insByIndex{$K2}{$i}[0]);
				$insObs++;
				if ( $ins1 eq $ins2 ) {
					$mSeq .= $ins1;
					@qIns1 = split('',$insByIndex{$K1}{$i}[1]);	
					@qIns2 = split('',$insByIndex{$K2}{$i}[1]);
					$qSeqNew = '';
					for($qIndex=0;$qIndex<length($ins1);$qIndex++ ) {
						$qSeqNew .= chr(max(ord($qIns1[$qIndex]),ord($qIns2[$qIndex])));
					}
					$qSeq .= $qSeqNew;
					$cigars .= 'I' x length($ins1);
				} elsif ( $ins2 =~ /$ins1/ ) {
					# 1 in 2
					$mSeq .= $ins1;
					$qSeq .= $insByIndex{$K1}{$i}[1];
					$cigars .= 'I' x length($ins1);
					$insErr++;
				} elsif ( $ins1 =~ /$ins2/ ) {
					# 2 in 1
					$mSeq .= $ins2;
					$qSeq .= $insByIndex{$K2}{$i}[1];
					$cigars .= 'I' x length($ins2);
					$insErr++;
				} else {
					$insErr++;
				}
			} elsif ( defined($insByIndex{$K1}{$i}) ) {
				if ( $i != $end ) {
					$w = $bases2[$i+1]
				} else {
					$w = '.';
				}

				# TO-DO: can ssw permit hanging insertions?
				if ( $y ne '.' && $w ne '.' ) {
					$insObs++; $insErr++;
				} else {
					$ins1 = lc($insByIndex{$K1}{$i}[0]);

					$mSeq .= $ins1;
					$qSeq .= $insByIndex{$K1}{$i}[1];
					$cigars .= 'I' x length($ins1);
				}
			} elsif ( defined($insByIndex{$K2}{$i}) ) {
				if ( $i != $end ) {
					$v = $bases1[$i+1]
				} else {
					$v = '.';
				}

				if ( $x ne '.' && $v ne '.' ) {
					$insObs++; $insErr++;
				} else {
					$ins2 = lc($insByIndex{$K2}{$i}[0]);

					$mSeq .= $ins2;
					$qSeq .= $insByIndex{$K2}{$i}[1];
					$cigars .= 'I' x length($ins2);
				}
			}
		}
		
		$qname = $a1[5];
		$mapq = int(avg($a1[6],$a2[6]));

		if ( ! $bowtieFormat ) {
			$qname =~ s/(.+?[_ ])[12](:.+)/${1}3${2}/;
		}
		print SAM $qname,"\t",'0',"\t",$REF_NAME,"\t",($start+1),"\t",$mapq;
		print SAM "\t",condenseCigar($cigars),"\t*\t0\t0\t",$mSeq,"\t",$qSeq,"\n";
	} else {
		$K = $pairs{$mID}{$mPairs[0]}[2];
		print SAM $sam[$K],"\n";
	}
}

$observations{$REF_NAME}{'obs'} = $obs;
$observations{$REF_NAME}{'fmv'} = $fmv;
$observations{$REF_NAME}{'tmv'} = $tmv;
$observations{$REF_NAME}{'dmv'} = $dmv;
$observations{$REF_NAME}{'insObs'} = $insObs;
$observations{$REF_NAME}{'insErr'} = $insErr;

close(SAM);
if ( $useStorable ) {
	store(\%observations,$ARGV[2].'.sto');
}
