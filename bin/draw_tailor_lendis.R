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

# part of tailor package
# used to draw length distribution
# from bed2 format given by tailor package
source (paste (Sys.getenv ("PIPELINE_DIRECTORY"),"/bin/Tailor.R",sep=""))
pkgTest ("ggplot2")
pkgTest ("gridExtra")
library ("grid")

argv  = commandArgs (TRUE)
table = read.table (argv[1], F)
name1 = argv[3]
name2 = argv[4]
colnames (table) = c("length", "tail", "count")
pdf (argv[2], onefile=TRUE, width=8.5, height=11)
table = transform (table, order = factor (tail, levels=c('*','A','C','G','T','the_others'), ordered=TRUE))
g1 = ggplot (table, aes (x=length, y=count, order=order)) + geom_bar (stat="identity", aes (fill=order)) + scale_fill_manual(values=c("black", "blue", "yellow", "darkgreen", "red", "darkgrey")) + theme_minimal()
g2 = ggplot (table, aes (x=length, y=count, fill=tail)) + geom_bar(stat="identity", position=position_dodge(width=0.5)) + scale_fill_manual(values=c("black", "blue", "yellow", "darkgreen", "red", "darkgrey")) + theme_minimal() + theme(legend.title=element_blank())
g3 = ggplot (table, aes (x=length, y=count)) + geom_bar(stat='identity') + facet_wrap(~tail)
grid.arrange(g1, g2, g3, ncol=1, as.table=TRUE, main = textGrob(paste (name1, name2), vjust = 1, gp = gpar(fontface = "bold", cex = 1) ))
gc = dev.off ()
