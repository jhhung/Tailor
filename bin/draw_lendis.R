#! /usr/bin/env Rscript
# this function is copied from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE)
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
}

pkgTest ("ggplot2")
pkgTest ("scales")

argv = commandArgs (T)
perfectMatchLendis = read.table (argv[1],F)
colnames (perfectMatchLendis) = c("length", "perfectMatchCounts")
tailedMatchLendis = read.table (argv[2],F)
colnames (tailedMatchLendis) = c("length", "prefixMatchCounts")

mergedLenDis = merge (perfectMatchLendis, tailedMatchLendis, all=T)
mergedLenDis[is.na (mergedLenDis)] = 0

main = argv[3]
pdf (paste (main, ".length_distribution.pdf", sep=''))
main=gsub ("\\."," ",main)
main=paste(strwrap(main, width = 50), collapse = "\n") 

# copied from http://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUp <- function(x, nice=c(1,2,4,5,6,8,10)) {
	if(length(x) != 1) stop("'x' must be of length 1")
	10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}

ru=roundUp( (max(mergedLenDis[,1])-min(mergedLenDis[,2]) )/10 ) 
barplot (t(as.matrix(mergedLenDis[,c(2,3)])), names = mergedLenDis[,1], legend=T, yaxt='n', cex.main=2, main=main, cex.names=1.5) 
+ axis (2, tck=0.01, lwd=3, at=seq(ru*-20, ru*20, ru), cex.axis=1.5)

dev.off ()