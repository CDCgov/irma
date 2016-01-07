#!/usr/bin/env Rscript
args = commandArgs(TRUE)
if (length(args) != 6) {
	cat("Usage:\n\tRscript ./sqmHeatmap.R <run> <gene> <COVG.txt> <VARS.txt> <STATS.txt> <out.pdf>\n")
	q()
}

run=args[1]
gene=args[2]
COVG=args[3]
VARS=args[4]
STAT=args[5]
pdfFile=args[6]

D=read.table(COVG,header=TRUE,sep="\t")
V=read.table(VARS,header=TRUE,sep="\t")
V$Minority_Allele=substr(V$Minority_Allele,1,1)
V$Consensus_Allele=substr(V$Consensus_Allele,1,1)

C=max(D$Coverage); Xrange=c(1,max(D$Position)); Vlen=nrow(V); H=C/Vlen*as.integer(rownames(V)); C2=C/2
Cols=vector()
#Cols[["A"]] = "#1F77B4"
#Cols[["C"]] = "#FF7F0E"
#Cols[["G"]] = "#2CA02C"
#Cols[["T"]] = "#D62728"
#Cols[["-"]] = "#FFFFFF"

pdf(pdfFile,width=10.5,height=8)
par(mfrow=c(2,1),mar=c(4, 5, .1, 1))
plot(D$Position,D$Coverage,col="black",xlim=Xrange,ylab="Coverage depth",xlab=paste(sep="",gene," position (",run,")"))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "black")
segments(D$Position,0,D$Position,D$Coverage,col="gray")
variants=vector(length=Vlen)
for(i in 1:Vlen ) { 
	d=which(D$Position==V$Position[i])
	c=D$Coverage[d]
	a=V$Minority_Allele[i]

	if ( a == "A" ) {
		color="#1F77B4"
	} else if ( a == "C" ) {
		color="#FF7F0E"
	} else if ( a == "G" ) {
		color="#2CA02C"
	} else if ( a == "T" ) {
		color="#D62728"
	} else {
		color="#FFFFFF"
	}
	Cols[i] = color

#	color=Cols[[ V$Minor.Allele[i] ]]
	if ( c < C2) {
		segments(D$Position[d],c,D$Position[d],C,col=color)
	} else {
		segments(D$Position[d],0,D$Position[d],c,col=color)
	}
	variants[i] = paste(D$Consensus[d],'2',V$Minority_Allele[i],sep='')
}
Sizes=V$Minority_Frequency+1
points(x=V$Position,y=H,,ylab="Nth Variant",xlim=Xrange,pch=as.character(V$Minority_Allele),xlab="Position",cex=Sizes,col="white")


if ( file.exists(STAT) ) {
	S=read.table(STAT,header=FALSE,sep="\t")
	EE=S[S$V2=="ExpectedErrorRate",3]
	bp=barplot(c(EE,V$Minority_Frequency),beside=TRUE,names.arg=c('exp. err',variants),ylab="Observed frequency",xlab="minor variants",col=c('black',Cols),las=2)
	text(bp,c(0,V$Minority_Frequency*.5),labels=c('',V$Position),col="white",cex=.75)
	abline(h=EE,col="#282828",lty=2,lwd=.75)
} else {
	bp=barplot(V$Minority_Frequency,beside=TRUE,names.arg=variants,ylab="Observed frequency",xlab="minor variants",col=Cols,las=2)
	text(bp,(max(V$Minority_Frequency)*.03),labels=c(V$Position),col="white",cex=.75)
}
dev.off()
