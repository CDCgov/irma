#!/usr/bin/env perl
# Sam Shepard -- merge SAM, call variants, create consensus, do coverage
# 9.2014


use Storable;
#use Getopt::Long;
#GetOptions(	'no-gap-allele|G' => \$noGap, );

if ( scalar(@ARGV) != 3 ) {
	$message = "Usage:\n\tperl $0 <ref> <sam> <prefix>\n";
#$message .= "\t\t-S|--sig-level <FLT>\t\t\tSignificance test (90, 95, 99, 99.9) variant is not machine error.\n";
	die($message."\n");
}

# FUNCTIONS #
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

sub sum($$$) {
        my $i;
        my $S = 0;
        for $i (0..($_[2]-1)) {
                $S += $_[0]->[$_[1]+$i];
        }
        return $S;
}
#############

$REisBase = qr/[ATCG]/;
$REgetMolID = qr/^(.+?)[_ ]([12]):.+/;

$/ = ">";
open(REF,'<',$ARGV[0]) or die("Cannot open $ARGV[0] for reading.\n");
while($record = <REF>) {
	chomp($record);
	@lines = split(/\r\n|\n|\r/, $record);
	$REF_NAME = shift(@lines);
	$REF_SEQ = join('',@lines);
	if ( length($REF_SEQ) < 1 ) {
		next;
	}
	$REF_LEN = length($REF_SEQ);
	last;
}
close(REF);
if ( !defined($REF_LEN) ) { die("No reference found.\n"); }

$/ = "\n";
open(SAM,'<',$ARGV[1]) or die("Cannot open $ARGV[1] for reading.\n");
while($line=<SAM>) {
	chomp($line);
	if ( $line eq '@' ) {
		next;
	}

	($qname,$flag,$rn,$pos,$mapq,$cigar,$mrnm,$mpos,$isize,$seq,$qual) = split("\t",$line);
	if ( $qname =~ $REgetMolID ) {
		$qMolID = $1;
		$qSide = $2;
	}

	if ( $REF_NAME eq $rn ) {
		$seq = uc($seq);
		@NTs = split('',$seq);
		@Qint = unpack("c* i*",$qual);
		@cigars = split('',$cigar);
		$rpos=$pos-1;
		$qpos=0;
		
		if ( $rpos > 0 ) {
			$aln = '.' x $rpos; 
		} else {
			$aln = '';
		}
		
		while($cigar =~ /(\d+)([MIDNSHP])/g ) {
			$inc=$1;
			$op=$2;
			if ( $op eq 'M' ) {
				for(1..$inc) {
					$aln .= $NTs[$qpos];
					$allele = substr($seq,$qpos,1);
					$table[0][$rpos]{$allele}++;
					$table[2][$rpos]{$allele} += $Qint[$qpos];
					$qpos++; $rpos++;
				}
			} elsif ( $op eq 'D' ) {
				$aln .= '-' x $inc;
				$table[4]{$rpos-1}{$inc}++;
				for(1..$inc) {
					$table[0][$rpos]{'-'}++;
					$rpos++;
				}
			} elsif( $op eq 'I' ) {
				$insert = lc(substr($seq,$qpos,$inc));
				$table[1]{$rpos-1}{$insert}++;
				$table[3]{$rpos-1}{$insert} += ( sum(\@Qint,$qpos,$inc) / length($insert) ) - 33;   #substr($qual,$qpos,$inc);
				$qpos += $inc;
			} elsif( $op eq 'S' ) {
				$qpos += $inc;
			} elsif( $op eq 'N' ) {
				# Analyzing merge pairs requires this change instead of the normal 'N' interpretation.
				$aln .= '.' x $inc;
				$rpos += $inc;
			} elsif ( $op eq 'H' ) {
				next;
			} else {
				die("Extended CIGAR ($op) not yet supported.\n");
			}
		}
		$aln .= '.' x (($REF_LEN)-$rpos);
		$table[5]{$aln}++;
	}
}
close(SAM);
store(\@table, $ARGV[2].'.sto');
