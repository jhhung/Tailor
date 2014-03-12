# part of tailor package
# used to draw length distribution
# from bed2 format given by tailor package
library (ggplot2)
library (gridExtra)
argv  = commandArgs (TRUE)
table = read.table (argv[1], F)
colnames (table) = c("length", "tail", "count")
pdf (argv[2], onefile=TRUE, width=8.5, height=11)
table$tail = factor ( table$tail, levels = c("*", "A","C","G","T","the_others"))
g1 = ggplot (table, aes (x=length, y=count, fill=tail)) + geom_bar (stat="identity") + scale_fill_manual(values=c("black", "blue", "yellow", "darkgreen", "red", "darkgrey")) + theme_minimal() 
g2 = ggplot (table, aes (x=length, y=count, fill=tail)) + geom_bar(stat="identity", position=position_dodge(width=0.5)) + scale_fill_manual(values=c("black", "blue", "yellow", "darkgreen", "red", "darkgrey")) + theme_minimal() + theme(legend.title=element_blank())
g3 = ggplot (table, aes (x=length, y=count)) + geom_bar(stat='identity') + facet_wrap(~tail)
grid.arrange(g1, g2, g3, ncol=1, as.table=TRUE)
gc = dev.off ()
