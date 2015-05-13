#!/usr/bin/env Rscript
args = commandArgs(TRUE)
if (length(args) != 2) {
	cat("Usage:\n\tRscript ./heuristicDiagram.R <ALL_ALLELES.txt> <out.pdf>\n")
	q()
}

D=read.table(args[1],sep="\t",header=TRUE)
pdf(args[2],width=10.5,height=8)
par(mfrow=c(4,1),mar=c(2,2,2,2))
plot(density(D$Average_Quality,from=min(D$Average_Quality,na.rm=TRUE),to=max(D$Average_Quality,na.rm=TRUE),na.rm=TRUE),main="Density of average allele quality"); abline(v=24,col="red")
plot(density(D$Frequency,from=min(D$Frequency,na.rm=TRUE),to=0.02,na.rm=TRUE),main="Density of observed frequency <= 2%"); abline(v=.005,col="red")
MX=quantile(D$Total,probs=(.20))
hist(D$Total[D$Total<=MX],breaks=50,main="Histogram of coverage (Depth <= 20% Quantile)",xlim=c(0,MX+1)); abline(v=100,col="red")
hist(D$ConfidenceNotMacErr[D$ConfidenceNotMacErr > 0],breaks=50,main="Histogram of confidence not machine error, non-zero"); abline(v=.5,col="red")
dev.off()
