#!/usr/bin/env perl
# Sam Shepard -- mergeSAMpairs -- version 2
# 9.2014

use 5.16.1;
use strict;
use warnings;
use Storable;
use Getopt::Long;

use constant { false => 0, true => 1 };

my ($useStorable, $bowtieFormat) = (false,false);
GetOptions(
		'use-storable|S' => \$useStorable, 'bowtie-format|B' => \$bowtieFormat
	);

if ( scalar(@ARGV) != 3 ) {
	my $message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
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

my $REisBase = qr/[ATCG]/;
my $REgetMolID = qr/(.+?)[_ ]([12]):.+/;

$/ = ">"; 
my $REF_LEN = 0;
my $REF_NAME = '';
my @REF_SEQ = ();
my $REF;
open($REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while(my $record = <$REF>) {
	chomp($record);
	my @lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	my $seq = join('',@lines);
	if ( length($seq) < 1 ) {
		next;
	}
	@REF_SEQ = split('',uc($seq));
	$REF_LEN = length($seq);
	last;
}
close($REF);

$/ = "\n";
my $SAM;
open($SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
my @sam = <$SAM>; chomp(@sam);
close($SAM);

my %pairs = ();
my %insByIndex = ();
foreach my $K ( 0 .. $#sam ) {
	if ( substr($sam[$K],0,1) eq '@' ) {
		next;
	}

	my ($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$sam[$K]);

    my ($qMolID, $qSide);
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
		my @NTs = split('',uc($seq));
		my @QCs = split('',$qual);
		my @Qint = unpack("c* i*",$qual);
		my @cigars = split('',$cigar);
		my $rpos=$pos-1;
		my $qpos=0;
        my ($aln, $qAln) = ('','');
		
		#if ( $rpos > 0 ) {
		#	$aln = '.' x $rpos; 
		#	$qAln = ' ' x $rpos;
		#} else {
		#	$aln = '';
		#	$qAln = '';
		#}
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			my $inc=$1;
			my $op=$2;
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
		#$aln .= '.' x (($REF_LEN)-$rpos);
		#$qAln .= ' ' x (($REF_LEN)-$rpos);
		$pairs{$qMolID}{$qSide} = [$aln,$qAln,$K,($pos-1),($rpos-1),$qname,$mapq];
	}
}

my %observations = ();
my $OSAM;
open($OSAM,'>',$ARGV[2].'.sam') or die("Cannot open $ARGV[2].sam\n");
foreach my $line ( @sam ) {
	if ( substr($line,0,1) eq '@' ) {
		print $OSAM $line,"\n";
	} else {
		last;
	}
}
my ($dmv, $obs, $fmv, $tmv) = (0,0,0,0);
my ($insObs, $insErr) = (0,0);
foreach my $mID ( keys(%pairs) ) {
	my @mPairs = keys(%{$pairs{$mID}});
	if ( scalar(@mPairs) == 2 ) {
		my @a1 = @{$pairs{$mID}{'1'}};			
		my @a2 = @{$pairs{$mID}{'2'}};
		my ($s1,$e1) = ($a1[3],$a1[4]);
		my ($s2,$e2) = ($a2[3],$a2[4]);

		my $start = min($s1,$s2);
		my $end = max($e1,$e2);

		my $mSeq = '';
		my $cigars = '';
		my $qSeq = '';
		
		my $K1 = $a1[2];
		my $K2 = $a2[2];
		my @bases1 = split('',$a1[0]);
		my @bases2 = split('',$a2[0]);
		my @quals1 = unpack("c* i*",$a1[1]);
		my @quals2 = unpack("c* i*",$a2[1]);

        
        my $FB1 = sub {
            my $i = $_[0];
            if ( $i < $s1 || $i > $e1 ) {
                return '.';
            } else {
                if ( ($i - $s1) > $#bases1 ) { die("Bad: " . ($i-$s1) . " > " . $#bases1 ."\n"); }
                return $bases1[$i - $s1];
            }
        };

        my $FQ1 = sub {
            my $i = $_[0];
            if ( $i < $s1 || $i > $e1 ) {
                return ' ';
            } else {
                return $quals1[$i - $s1];
            }
        };

        my $FB2 = sub {
            my $i = $_[0];
            if ( $i < $s2 || $i > $e2 ) {
                return '.';
            } else {
                return $bases2[$i - $s2];
            }
        };

        my $FQ2 = sub {
            my $i = $_[0];
            if ( $i < $s2 || $i > $e2 ) {
                return ' ';
            } else {
                return $quals2[$i - $s2];
            }
        };

		foreach my $i ( $start .. $end ) {
			#my $x = $bases1[$i];
			#my $y = $bases2[$i];
			#my $qx = $quals1[$i];
			#my $qy = $quals2[$i];
			my $x   = $FB1->($i);
			my $y   = $FB2->($i);
			my $qx  = $FQ1->($i);
			my $qy  = $FQ2->($i);
			my $r   = $REF_SEQ[$i];

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
				my $ins1 = lc($insByIndex{$K1}{$i}[0]);
				my $ins2 = lc($insByIndex{$K2}{$i}[0]);
				$insObs++;
				if ( $ins1 eq $ins2 ) {
					$mSeq .= $ins1;
					my @qIns1 = split('',$insByIndex{$K1}{$i}[1]);	
					my @qIns2 = split('',$insByIndex{$K2}{$i}[1]);
					my $qSeqNew = '';
					foreach my $qIndex ( 0 .. (length($ins1)-1) ) {
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
                my $w = '';
				if ( $i != $end ) {
					#$w = $bases2[$i+1]
					$w = $FB2->($i+1);
				} else {
					$w = '.';
				}

				# TO-DO: can ssw permit hanging insertions?
				if ( $y ne '.' && $w ne '.' ) {
					$insObs++; $insErr++;
				} else {
					my $ins1 = lc($insByIndex{$K1}{$i}[0]);

					$mSeq .= $ins1;
					$qSeq .= $insByIndex{$K1}{$i}[1];
					$cigars .= 'I' x length($ins1);
				}
			} elsif ( defined($insByIndex{$K2}{$i}) ) {
                my $v = '';
				if ( $i != $end ) {
					#$v = $bases1[$i+1];
					$v = $FB1->($i+1);
				} else {
					$v = '.';
				}

				if ( $x ne '.' && $v ne '.' ) {
					$insObs++; $insErr++;
				} else {
					my $ins2 = lc($insByIndex{$K2}{$i}[0]);

					$mSeq .= $ins2;
					$qSeq .= $insByIndex{$K2}{$i}[1];
					$cigars .= 'I' x length($ins2);
				}
			}
		}
		
		my $qname = $a1[5];
		my $mapq = int(avg($a1[6],$a2[6]));

		if ( ! $bowtieFormat ) {
			$qname =~ s/(.+?[_ ])[12](:.+)/${1}3${2}/;
		}
		print $OSAM $qname,"\t",'0',"\t",$REF_NAME,"\t",($start+1),"\t",$mapq;
		print $OSAM "\t",condenseCigar($cigars),"\t*\t0\t0\t",$mSeq,"\t",$qSeq,"\n";
	} else {
		my $K = $pairs{$mID}{$mPairs[0]}[2];
		print $OSAM $sam[$K],"\n";
	}
}

$observations{$REF_NAME}{'obs'} = $obs;
$observations{$REF_NAME}{'fmv'} = $fmv;
$observations{$REF_NAME}{'tmv'} = $tmv;
$observations{$REF_NAME}{'dmv'} = $dmv;
$observations{$REF_NAME}{'insObs'} = $insObs;
$observations{$REF_NAME}{'insErr'} = $insErr;

close($OSAM);
if ( $useStorable ) {
	store(\%observations,$ARGV[2].'.sto');
}
