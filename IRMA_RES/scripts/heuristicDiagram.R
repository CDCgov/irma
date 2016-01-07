#!/usr/bin/env Rscript
args = commandArgs(TRUE)
if (length(args) != 6) {
	cat("Usage:\n\tRscript ./heuristicDiagram.R <MIN_AQ> <MIN_F> <MIN_TCC> <MIN_CONF> <ALL_ALLELES.txt> <out.pdf>\n")
	q()
}

MIN_AQ=as.double(args[1])
MIN_F=as.double(args[2])
MIN_TCC=as.double(args[3])
MIN_CONF=as.double(args[4])

D=read.table(args[5],sep="\t",header=TRUE)
pdf(args[6],width=10.5,height=8)
par(mfrow=c(3,2),mar=c(2,2,2,2))

# Quality density plot
plot(density(D$Average_Quality,from=min(D$Average_Quality,na.rm=TRUE),to=max(D$Average_Quality,na.rm=TRUE),na.rm=TRUE),main="Density of average allele quality")
abline(v=MIN_AQ,col="red")
plot(density(D$Average_Quality,from=min(D$Average_Quality,na.rm=TRUE),to=MIN_AQ,na.rm=TRUE),main=paste(sep=""," to ",MIN_AQ))

# Frequency density plot
plot(density(D$Frequency,from=0.0,to=0.10,na.rm=TRUE),main="Density of observed frequency (to 10%)")
abline(v=MIN_F,col="red")
plot(density(D$Frequency,from=0.0,to=MIN_F,na.rm=TRUE),main=paste(sep=""," to ",MIN_F,"%"))

# Coverage depth histogram
MX=quantile(D$Total,probs=(.20))
hist(D$Total[D$Total<=MX],breaks=50,main="Histogram of coverage (Depth <= 20% Quantile)",xlim=c(0,MX+1))
abline(v=MIN_TCC,col="red")

# Confidence value histogram
hist(D$ConfidenceNotMacErr[D$ConfidenceNotMacErr > 0],breaks=50,main="Histogram of confidence not machine error, non-zero")
abline(v=MIN_CONF,col="red")
dev.off()
