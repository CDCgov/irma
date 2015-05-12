#!/bin/bash

if [ $# -eq 2 ];then
	MODULE=$1
	RUN=$2
	PAIRED=1
elif [ $# -eq 3 ];then
	MODULE=$1
	RUN=$2			# name of the run
	PAIRED=$3
else
	echo -e "Usage:\n\t$0 <module-param> <run-dir> [paired:0/1]\n"
	exit 1
fi

# POSSIBLY RELEVANT
INS_T=0.15		# threshold for insertion refinement of references
DEL_T=0.75		# threshold for deletion refinement of references
MIN_C=2			# minimum count for variants
MIN_F=0.0075		# minimum frequency for variants
MIN_FI=0.0075		# minimum frequency for insertion variations
MIN_FD=0.0075		# minimum frequency for non-singleton deletion variations
MIN_AQ=24		# minimum average variant quality, does not apply to deletions
MIN_TCC=100		# minimum non-ambiguous column coverage
MIN_CONF=0.80		# minimum confidence not machine error
MIN_AMBIG=0.25		# minimum ambiguous nucleotides for amended consensus (with ambig codes)
SIG_LEVEL=0.999		# significance test level for variant calling (.90,.95,.99,.999). 
AUTO_F=1		# auto-adjust frequency threshold
NONSEGMENTED=0		# segmented versus non-segmented virus

# DO NOT ALTER
owd=$(pwd)
bpath=
if [ "$bpath" == "" ]; then
	bpath=$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
fi
rpath=$bpath/IRMA_RES	# IRMA resources
spath=$rpath/scripts	# IRMA specific scripts/binaries
LABEL=$(which LABEL)
Lspath=$(dirname $LABEL)/LABEL_RES/scripts
Lmpath=$(dirname $LABEL)/LABEL_RES/training_data/IRMA
PARALLEL="$spath/parallel --will-cite"

# check for module
MODULE2=$(echo $MODULE|cut -f1 -d'-')
if [ ! -d $rpath/modules/$MODULE2 ];then
	time_stamp "Error: $MODULE2 not found."
	exit 1
else
	mpath=$rpath/modules/$MODULE2	
fi

# load module parameter file
source $mpath/params-$MODULE.sh
ppath=$RUN/redo
if [ ! -d $ppath ];then
	mkdir $ppath
fi

# stage paths
if [ "$AUTO_F" -eq 1 ];then
	AUTO_F="-A"
else
	AUTO_F=""
fi

for bam in $RUN/*.bam;do
	gene=$(basename $bam .bam)
	ref=$RUN/$gene.fasta
	final=$(dirname $bam)/tables/$gene
	fix=$ppath/$gene

	if [ -r ${final}-insertions.txt -a -r ${final}-deletions.txt ];then
		callPhaseArgs="-C $MIN_C -F $MIN_F -I $MIN_FI -D $MIN_FD -Q $MIN_AQ -T $MIN_TCC -M $MIN_CONF -S $SIG_LEVEL $AUTO_F"
		$spath/vcfGenerator.pl $callPhaseArgs $ref ${final}-allAlleles.txt ${final}-insertions.txt ${final}-deletions.txt -V $fix-variants.txt > $fix.vcf
	fi
done
