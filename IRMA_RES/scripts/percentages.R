#!/usr/bin/env Rscript
# Sam Shepard

args = commandArgs(TRUE)
if (length(args) != 4) {
        cat("Usage:\n\tRscript ./percentages.R <run> <COUNTS.txt> <OUTPUT.pdf> <PAIRED:1/0>\n")
        q()
}

tableau10=c("#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")

run=args[1]
countsFile=args[2]
outputFile=args[3]
paired=as.numeric(args[4])

pdf(outputFile,width=12,height=11)
par(mfrow=c(2,2),mar=c(3, 3, 3, 3))

D=read.table(header=TRUE,sep="\t",file=countsFile)

total=D$Reads[D$Record =='1-initial']
#qual=D$Reads[D$Record =='2-passQC']
nonqual=D$Reads[D$Record =='2-failQC']
assembled=D$Reads[D$Record =='3-match']
therest=total-nonqual-assembled

vals=c(assembled,nonqual,therest)
S=vector(length=length(vals))
for( i in 1:length(vals)) {
	v=vals[i]
	if ( v >= 1000000) { 
		s=sprintf("%.1fM",v/1000000); 
	} else if ( v >= 1000) {
		s=sprintf("%.1fk",v/1000); 
	} else { 
		s=sprintf("%d",v); 
	} 
	S[i]=s
}
perc=sprintf("%.1f%% (%s)", vals/total*100,S)
grps=c("Assembled","QC filtered","Other")

cols=tableau10[c(1,2,4)]

if ( paired ) {
	pie(x=vals,labels=perc,col=cols,main="1. Percentages of total reads (R1 + R2)")
} else {
	pie(x=vals,labels=perc,col=cols,main="1. Percentages of total reads")
}
legend("bottomright", grps, fill=cols)


assembled=D$Patterns[D$Record =='3-match']
total=D$Patterns[D$Record =='2-passQC']
if ( sum(D$Record == '3-chimeric') > 0 ) {
	chimeric=sum(D$Patterns[D$Record =='3-chimeric'])
} else {
	chimeric=0
}

if ( sum(D$Record=='3-unrecognizable') > 0 ) {
	unused=D$Patterns[D$Record =='3-unrecognizable']
} else {
	unused=0
}

if ( sum(D$Record=='3-altmatch') > 0 ) {
	unused=unused+sum(D$Patterns[D$Record=='3-altmatch'])
}
unmatched=D$Patterns[D$Record =='3-nomatch']

cols=c("#1F77B4","#ed9a9a","#e05758","#d62728")
vals=c(assembled,unused,chimeric,unmatched)
S=vector(length=length(vals))
for( i in 1:length(vals)) {
	v=vals[i]
	if ( v >= 1000000) { 
		s=sprintf("%.1fM",v/1000000); 
	} else if ( v >= 1000) {
		s=sprintf("%.1fk",v/1000); 
	} else { 
		s=sprintf("%d",v); 
	} 
	S[i]=s
}
perc=sprintf("%.1f%% (%s)", vals/total*100,S)
grps=c("Assembled","Unusable","Chimeric","No match")
pie(x=vals,labels=perc,col=cols,main="2. Percentages of all read patterns passing QC")
legend("bottomright", grps, fill=cols)


G=D[grep(pattern="^4-",x=D$Record),]
G$Record=gsub(pattern="^4-",replacement="",x=G$Record)

o=order(G$Record)
grps=G$Record[o]
vals=G$PairsAndWidows[o]
S=vector(length=length(vals))
for( i in 1:length(vals)) {
	v=vals[i]
	if ( v >= 1000000) { 
		s=sprintf("%.1fM",v/1000000); 
	} else if ( v >= 1000) {
		s=sprintf("%.1fk",v/1000); 
	} else { 
		s=sprintf("%d",v); 
	} 
	S[i]=s
}
total=sum(vals)
cols=tableau10[1:8]

perc=vector(length=length(grps))
for(i in 1:length(grps) ) {
	percentage=vals[i]/total*100
	if ( percentage < 2 ) {
		perc[i] = grps[i]
	} else {
		perc[i] = sprintf("%.1f%% (%s)\n%s", vals[i]/total*100,S[i],grps[i])
	}

}

if ( paired ) {
	pie(x=vals,labels=perc,col=cols,main="3. Percentages of assembled, merged-pair reads")
} else {
	pie(x=vals,labels=perc,col=cols,main="3. Percentages of assembled reads")
}

plot.default(c(0,250),c(0,250), type="n", axes=FALSE,main=paste(sep="","SAMPLE \"",run,"\""),ylab="", xlab="")
if ( paired ) {
text(20, 30,
"READ PROPORTIONS.

1. Percentages of total read counts (R1 & R2)
    - ASSEMBLED: influenza reads in final assemblies.
    - QC FILTERED: didn't pass length/median quality thresholds.
    - OTHER: non-flu and contaminant/poor flu signal.

2. Percentages of all read patterns passing QC process
   - Patterns are clustered or non-redundant reads.
   - ASSEMBLED: excellent influenza read patterns.
   - UNUSABLE: poor or contaminant flu patterns.
   - CHIMERIC: flu patterns matching both strands.
   - NO MATCH: non-flu read patterns.

3. Percentages of assembled, merged-pair read counts
   - Shows the proportion of gene segments to the genome.
   - Paired-end reads have been merged into a single count
     unless not applicable: single-end reads have been used.", adj=c(0,0))
} else {
text(20, 30,
"READ PROPORTIONS.

1. Percentages of total read counts
    - ASSEMBLED: influenza reads in final assemblies.
    - QC FILTERED: didn't pass length/median quality thresholds.
    - OTHER: non-flu and contaminant/poor flu signal.

2. Percentages of all read patterns passing QC process
   - Patterns are clustered or non-redundant reads.
   - ASSEMBLED: excellent influenza read patterns.
   - UNUSABLE: poor or contaminant flu patterns.
   - CHIMERIC: flu patterns matching both strands.
   - NO MATCH: non-flu read patterns.

3. Percentages of assembled read counts
   - Shows the proportion of gene segments to the genome.
   - Paired-end reads have been merged into a single count
     unless not applicable: single-end reads have been used.", adj=c(0,0))

}
dev.off()
