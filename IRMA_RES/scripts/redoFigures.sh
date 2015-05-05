#!/bin/bash

# SET to install

bpath=~/Tools
spath=$bpath/IRMA_RES/scripts
if [ $# -eq "2" ];then
	PAIRED=$2
else
	PAIRED=0
fi
ppath=$1
fpath=$ppath/figures
apath=$ppath/intermediate/?-ASSEMBLE_*
logpath=$ppath/logs
tpath=$ppath/tables
sqmpath=$ppath/matrices

for i in $ppath/*fasta;do
	gene=$(basename $i .fasta)
	echo "DOING $gene"
	n=$(wc -l < $tpath/${gene}-variants.txt)
	if [ $n -gt 2 ]; then
		$spath/sqmHeatmap.R $sqmpath/$gene-EXPENRD.sqm $fpath/$gene-EXPENRD.pdf 2
		$spath/sqmHeatmap.R $sqmpath/$gene-JACCARD.sqm $fpath/$gene-JACCARD.pdf 2
		$spath/sqmHeatmap.R $sqmpath/$gene-MUTUALD.sqm $fpath/$gene-MUTUALD.pdf 2
		$spath/sqmHeatmap.R $sqmpath/$gene-NJOINTP.sqm $fpath/$gene-NJOINTP.pdf 2
	fi

	final=$tpath/$gene
	if [ $n -gt 1 ]; then
		$spath/coverageDiagram.R $final-coverage.txt $final-variants.txt $final-pairingStats.txt $fpath/$gene-coverageDiagram.pdf
	else
		$spath/simpleCoverageDiagram.R $final-coverage.txt $fpath/$gene-coverageDiagram.pdf
	fi
	$spath/heuristicDiagram.R ${final}-allAlleles.txt $fpath/$gene-heuristics.pdf
done
readSize=$(expr $(zcat $apath/reads.tar.gz |wc -l) / 4)
echo "Percents"
$spath/percentages.R $ppath/tables/READ_COUNTS.txt $ppath/figures/READ_PERCENTAGES.pdf $PAIRED
