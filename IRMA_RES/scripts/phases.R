#!/usr/bin/env Rscript
# Sam Shepard - 2017-06-06
# Use single-linkage clustering to get phase groups based on distance cutoff.

args = commandArgs(TRUE)
if (length(args) > 2 || length(args) < 1 ) {
        cat("Usage:\n\tRscript ./phases.R <SQM> [cutoff]\n")
        q()
}

matFromSQM = function(sqmFile, s=2) {
	sqm=read.table(sqmFile,header=FALSE,sep="\t")
	mat=as.matrix(sqm[,s:length(sqm)])
	rownames(mat)=sqm[,1]
	colnames(mat)=sqm[,1]
	return(mat)
}

groupsByTree = function(mat, cutoff = 0.78 ) {
	H=hclust(as.dist(mat),method="single")
	G=cutree(H,h=cutoff)
	return(G)
}

filename = args[1]
if ( file.exists(filename) ) {
	if ( length(args) > 1 ) {
		Groups = groupsByTree(matFromSQM(filename), as.numeric(args[2]) )
	} else {
		Groups = groupsByTree(matFromSQM(filename))
	}
	write.table(Groups, file="", col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")
}
