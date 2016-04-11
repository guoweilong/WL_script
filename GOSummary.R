#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2015-08-29

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-t", "--title"), dest = "title", default = "", 
              help = "[opt] text shown on title"),
  make_option(c("-w", "--width"), dest = "png_width", default = 800, 
              help = "[opt] width (in px) of png"),
  make_option(c("--height"), dest = "png_height", default = 800, 
              help = "[opt] height (in px) of png"),
  make_option(c("-s", "--site"), dest = "sitefile", default = "", 
              help = "[opt] file of site to be marked")
)
#
parser <- OptionParser(usage = "%prog [options] file", 
   option_list=option_list, description = "Author: Guo, Weilong; guoweilong@126.com; 2015-11-23 \
Description: \
Input Ex: \
	cut -f2,3,4,5,10,13 GavidGO.result \
	Term	Count	%	PValue	Fold Enrichment	FDR \
	membrane	98	5.730994152	7.00E-10	1.7081798142	8.63E-07
   "
)
#
arguments <- parse_args(parser)
opt <- arguments$options
infile <- arguments$infile
#
if(infile == "") {
  print_help(parser)
  stop(sprintf("input file is not specified"))
}  
if(file.access(infile) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", infile))
}
#
outfile = arguments$outfile
if( outfile == "") {
  outfile = paste( infile, "png", sep=".")
}
#
#infile = "Yellow_module_David_ZhuList_simple.txt"
T = read.table(infile, header=TRUE, check.names = FALSE, sep = "\t")
NR = nrow(T)

#plot( -log(T[,4], base = 10), c(1:NR),
#      pwd = log(T[,3], base = 2),
#      xlab = "-log2(p-value)", ylab="", yaxt="n")


library("gplots")
colpool = rev(heat.colors(101))
MinusLogPv = -log(T[,6]/100, base=10)
MAX = max(MinusLogPv)
MAX_X = max(T[,5])*1.1
MAX_Y = NR*1.1
MinusLogPv = as.integer(MinusLogPv/MAX*100)
png(outfile, width=arguments$png_width, height=arguments$png_height)
#postscript("YellowModule_GOsummary.eps", width=8, height=9)
par(mai=c(1, 8, 1, 0.5) )
plot(NA, NA, ylim=c(0, MAX_Y), xlim=c(0, MAX_X),
     ylab="", xlab="Fold enrichment",  yaxt="n",
     cex=1.5, cex.axis=1.5, cex.lab=1.5)
symbols(x=T[,5], y=c(1:NR), 
        #circles=as.integer(log(T[,2], base=2) ), 
        circles=log(T[,2], base=2)/5, 
        #cex = 0.2,
        bg = colpool[MinusLogPv],
        inches=1/3, ann=F, fg=NULL, add=TRUE,
        xlab = "", ylab="", yaxt="n")

#mtext("enrichment (%)", line=2, side=1, cex=1.5)
axis(2, at=c(1:NR), labels=T[,1], 
     #col.axis="red",
     cex = 1.5, cex.axis = 1.5,
     las=2)
# Draw the color legends
seq(1, ceiling(MAX))
text = rep("",11)
text[c(1,6,11)] = rev(as.character(c(0, ceiling(MAX)/2.0, ceiling(MAX)) ))
pos=legend("topright", fill = rev(colpool)[seq(1,101,10)],
       border = NA, text, bty = "n", cex = 1.5,
       y.intersp=0.5, x.intersp = 0.8)
text(x=pos$rect$left+pos$rect$w/2-1, 
     y=pos$rect$top-pos$rect$h-0.5,
     "-log10(FDR)", cex = 1.5)
# text the size of group
text(T[,5], c(1:NR), labels = t(T[,2]), col="black", cex = 1.2)
dev.off()


