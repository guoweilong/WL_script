#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2016-04-11

# Argument
suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="[opt] input file, use STDIN if omitted"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file, use STDOUT if omitted")
)
#
parser <- OptionParser(usage = "%prog [options] file", 
     option_list=option_list, description = "Author: Guo, Weilong; guoweilong@126.com; 2015-05-20 \
Description: \
  Transpose an input matrix by row-to-column.\
Example: \
  Transpos.R -i input.xls -o output.xls \
-Input matrix: \
    M rows X N columns ==> N rows X M colums\
              "
)
#
arguments <- parse_args(parser)
opt <- arguments$options
# infile
infile <- arguments$infile
if(infile == "") {
  #print_help(parser)
  #stop(sprintf("input file is not specified"))
  infile = file("stdin")
}  
#if(file.access(infile) == -1) {
#  stop(sprintf("Specified file ( %s ) does not exist", infile))
#}
# outfile
outfile = arguments$outfile
if( outfile == "") {
  outfile = stdout()
}
# Read the input file
T = read.table( infile, header=FALSE, sep = "\t", check.names = FALSE)
write.table(t(T), file = outfile, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = FALSE)
#
