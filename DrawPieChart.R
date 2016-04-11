#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2015-11-20

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-t", "--title"), dest = "title", default = "", 
              help = "[opt] text shown on title"),
  make_option(c("-w", "--width"), dest = "png_width", default = 600, 
              help = "[opt] width (in px) of png"),
  make_option(c("--height"), dest = "png_height", default = 600, 
              help = "[opt] height (in px) of png"),
  make_option(c("-c", "--col"), dest = "col4value", default = 2, 
              help = "The column for values"),
  make_option(c("-n", "--name"), dest = "col4name", default = 1, 
              help = "The column for names")
)
#
parser <- OptionParser(usage = "%prog [options] file", 
   option_list=option_list, description = "Author: Guo, Weilong; guoweilong@126.com; 2015-11-20 \
Description: Draw the pie chart for the value\
Input examples: \
  chr1	3279\
  chr2	1980"
)
#
arguments <- parse_args(parser)
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
col4value = arguments$col4value
col4name = arguments$col4name
#
#infile = "Yellow_module_David_ZhuList_simple.txt"
T = read.table(infile, header=FALSE, check.names = FALSE, sep = "\t")
NR = nrow(T)
Data=T[,col4value]
names(Data)=T[,col4name]
png(outfile, width = arguments$png_width, height = arguments$png_height)
par( mai=c(0, 0.5, 0.5, 0.5) )
pie(Data, init.angle = 90, 
    main = arguments$title,
    cex = 1.5, cex.lab = 1.5, cex.main = 2,
    col = terrain.colors(NR) 
)
dev.off()
  


