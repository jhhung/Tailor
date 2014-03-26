
# tailor pipeline
# Bo W Han (bowhan@me.com)

# this function is modified from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE, lib=paste(Sys.getenv ("PIPELINE_DIRECTORY"),"Rlib",sep="/"), repos='http://cran.us.r-project.org')
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
}

pkgTest("gplots")
pkgTest("multicore")

draw_microRNA_tailing_balloon = function (t1, name, outDir) {

	pdf (paste (outDir, '/', name, t1$V1[1], ".miRNATailingBalloonPlot.pdf", sep=''), family="Helvetica")
	par(mfrow=c(3,3), cex.main=0.75, mar=c(0.5,0.5,0.5,0.5), oma=c(1.5,2,1,1))

	normalized_signal = sum(t1[,4]/t1[,5])
	if (normalized_signal > 0) {

		t1[,4] = 100*(t1[,4]/t1[,5])/normalized_signal ;

		five_three = xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V3)
		balloonplot (five_three, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
		five_tail = xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V9)
		balloonplot (five_tail, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
		three_tail = xtabs (as.integer (t1$V4) ~ t1$V3 + t1$V9)
		balloonplot (three_tail, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )

		t2 = t1[t1$V9==0,]
		if (nrow (t2) > 0) {
			five_three = xtabs (as.integer (t2$V4) ~ t2$V2 + t2$V3)
			balloonplot (five_three, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
			five_tail = xtabs (as.integer (t2$V4) ~ t2$V2 + t2$V9)
			balloonplot (five_tail, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
			three_tail = xtabs (as.integer (t2$V4) ~ t2$V3 + t2$V9)
			balloonplot (three_tail, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )
		}

		t3 = t1[t1$V9!=0,]
		if (nrow (t3) > 0) {
			five_three = xtabs (as.integer (t3$V4) ~ t3$V2 + t3$V3)
			balloonplot (five_three, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
			five_tail = xtabs (as.integer (t3$V4) ~ t3$V2 + t3$V9)
			balloonplot (five_tail, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
			three_tail = xtabs (as.integer (t3$V4) ~ t3$V3 + t3$V9)
			balloonplot (three_tail, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="white",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )
		}
	}
	invisible(dev.off())
}

argv = commandArgs (TRUE)
mirRelativePos = read.table (argv[1], F, sep="\t", stringsAsFactors=F)
numOfCore = argv[2]
name = argv[3]
outDir = argv[4]
mirRelativePosSplitted = split (mirRelativePos, mirRelativePos$V1)
mclapply (mirRelativePosSplitted, draw_microRNA_tailing_balloon, mc.cores=numOfCore, name = name, outDir = outDir)


