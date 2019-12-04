#!/usr/bin/env Rscript
args = commandArgs(TRUE)
if (length(args) < 2) {
	cat("Usage:\n\tRscript ./sqmHeatmap.R <sqm> <pdf>\n")
	q()
}

if ( length(args) == 3 ) {
	s=args[3]
} else {
	s=3
}

sqmFile=args[1]
pdfFile=args[2]
sqm=read.table(sqmFile,header=FALSE,sep="\t")
mat=as.matrix(sqm[,s:length(sqm)])
rownames(mat)=sqm[,1]
colnames(mat)=sqm[,1]
pdf(pdfFile,width=8.5,height=8.5)
title=paste("Variant site clusters,",basename(args[1]))
heatmap(mat,symm=TRUE,revC=TRUE,main=title,margins=c(10,10),col=heat.colors(200))
cat("Saved '",pdfFile,"'\n",sep="")
