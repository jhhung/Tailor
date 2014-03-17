
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
	sum4th = sum(t1[,4])	
	fivePrimeArm_het =  xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V3)	
	pdf (paste (outDir, '/', name, t1$V1[1], ".miRNATailingBalloonPlot.pdf", sep=''), family="Helvetica")
#	par (mfrow=c(2,2),mar=c(5,2,2,1))
	balloonplot (threePrimeArm_mut, main=main,xlab="5'",ylab="3'",sorted=T,label.size=.6,text.size=.5,rowmar=1,show.zeros=T, dotcolor="lightgreen")
	invisible(dev.off())
}

argv = commandArgs (TRUE)
mirRelativePos = read.table (argv[1], F, sep="\t", stringsAsFactors=F)
numOfCore = argv[2]
name = argv[3]
outDir = argv[4]
final_pdf = argv[5]
mirRelativePosSplitted = split (mirRelativePos, mirRelativePos$V1)
mclapply (mirRelativePosSplitted, draw_microRNA_tailing_balloon, mc.cores=numOfCore, name = name, outDir = outDir)



