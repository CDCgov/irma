#!/usr/bin/env Rscript
args = commandArgs(TRUE)
if (length(args) != 4) {
	cat("Usage:\n\tRscript ./sqmHeatmap.R <RUN> <GENE> <COVG.txt> <out.pdf>\n")
	q()
}

run=args[1]
gene=args[2]
COVG=args[3]
pdfFile=args[4]

D=read.table(COVG,header=TRUE,sep="\t")
C=max(D$Coverage); Xrange=c(1,max(D$Position)); C2=C/2
Cols=vector()
Cols[["A"]] = "#1F77B4"
Cols[["C"]] = "#FF7F0E"
Cols[["G"]] = "#2CA02C"
Cols[["T"]] = "#D62728"
Cols[["-"]] = "#FFFFFF"

pdf(pdfFile,width=10.5,height=8)
plot(D$Position,D$Coverage,col="black",xlim=Xrange,ylab="Coverage depth",xlab=paste(sep="",gene," position (",run,")"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
segments(D$Position,0,D$Position,D$Coverage,col="gray")
dev.off()
