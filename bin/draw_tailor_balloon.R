# Tailor, a BWT-based aligner for non-templated RNA tailing
# Copyright (C) 2014 Min-Te Chou, Bo W Han, Jui-Hung Hung
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# this function is modified from http://stackoverflow.com/questions/9341635/how-can-i-check-for-installed-r-packages-before-running-install-packages
source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/Tailor.R",sep=""))
pkgTest("gplots")
pkgTest("parallel")

draw_microRNA_tailing_balloon = function (t1, name, outDir) {

	pdf (paste (outDir, '/', name, t1$V1[1], ".miRNATailingBalloonPlot.pdf", sep=''), family="Helvetica")
	par(mfrow=c(3,3), cex.main=0.75, mar=c(0.5,0.5,0.5,0.5), oma=c(1.5,2,1,1))

	normalized_signal = sum(t1[,4]/t1[,5])
	if (normalized_signal > 0) {

		t1[,4] = 100*(t1[,4]/t1[,5])/normalized_signal ;

		five_three = xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V3)
		balloonplot (five_three, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
		five_tail = xtabs (as.integer (t1$V4) ~ t1$V2 + t1$V9)
		balloonplot (five_tail, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="5'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
		three_tail = xtabs (as.integer (t1$V4) ~ t1$V3 + t1$V9)
		balloonplot (three_tail, main= paste (name, t1$V1[1], "perfect+tailed", sep=" "), xlab="3'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )

		t2 = t1[t1$V9==0,]
		if (nrow (t2) > 0) {
			five_three = xtabs (as.integer (t2$V4) ~ t2$V2 + t2$V3)
			balloonplot (five_three, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
			five_tail = xtabs (as.integer (t2$V4) ~ t2$V2 + t2$V9)
			balloonplot (five_tail, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="5'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
			three_tail = xtabs (as.integer (t2$V4) ~ t2$V3 + t2$V9)
			balloonplot (three_tail, main= paste (name, t1$V1[1], "perfect", sep=" "), xlab="3'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )
		}

		t3 = t1[t1$V9!=0,]
		if (nrow (t3) > 0) {
			five_three = xtabs (as.integer (t3$V4) ~ t3$V2 + t3$V3)
			balloonplot (five_three, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="5'",ylab="3'",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="red",label.lines=T, )
			five_tail = xtabs (as.integer (t3$V4) ~ t3$V2 + t3$V9)
			balloonplot (five_tail, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="5'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="blue",label.lines=T, )
			three_tail = xtabs (as.integer (t3$V4) ~ t3$V3 + t3$V9)
			balloonplot (three_tail, main= paste (name, t1$V1[1], "tailed", sep=" "), xlab="3'",ylab="tail len",sorted=T,label.size=1, label.color="darkgrey",text.color="black",text.size=1.5,rowmar=0.5,colmar=0.5,show.zeros=T, dotcolor="green",label.lines=T, )
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


