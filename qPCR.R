#!/usr/bin/env Rscript
# Guo, Weilong; guoweilong@126.com; 2015-07-30

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
  make_option(c("-i", "--infile"), dest = "infile", default = "",
              help="input file"),
  make_option(c("-o", "--outfile"), dest = "outfile", default = "",
              help = "[opt] output file"),
  make_option(c("-w", "--width"), dest = "png_width", default = 800, 
              help = "[opt] width (in px) of png"),
  make_option(c("--height"), dest = "png_height", default = 800, 
              help = "[opt] height (in px) of png")
)

parser <- OptionParser(usage = "%prog [options] file", 
                       option_list=option_list, description = "Author: Guo, Weilong; guoweilong@126.com; 2015-07-30 \
                       Description: \
                        Show the statistic for qPRC \
                       Example: \
                        qPCR -i input -o output.png 
                       "
)

arguments <- parse_args(parser)
opt <- arguments$options
infile <- arguments$infile

if(infile == "") {
  print_help(parser)
  stop(sprintf("input file is not specified"))
}  
if(file.access(infile) == -1) {
  stop(sprintf("Specified file ( %s ) does not exist", infile))
}

outfile = arguments$outfile

if( outfile == "") {
  outfile = paste( infile, "png", sep=".")
}

png(outfile, 
    height=as.integer(arguments$png_height),
    width=as.integer(arguments$png_width) )

#infile="P VALUE_miR-16-5p control.txt"
infile="test.txt"
T = read.table(infile, header=TRUE, sep="\t", check.names=FALSE)

N_Sample = ncol(T)-2
Sample = T[,2]
SampleList = levels(Sample)
NSample = length(SampleList)

qPCR = 2^T[,c(-1,-2)]

EachMean <- function ( vec ){
  x1 = mean(vec[1:3])
  x2 = mean(vec[4:6])
  #x3 = mean(vec[7:9])
  #c(x1, x2, x3)
  c(x1, x2)
}


EachSd <- function ( vec ){
  x1 = sd(vec[1:3])
  x2 = sd(vec[4:6])
  #x3 = sd(vec[7:9])
 # c(x1, x2, x3)
  c(x1, x2)
}

EachPv <- function ( vec ){
  x12 = t.test(vec[1:3], vec[4:6])$p.value
  #x23 = t.test(vec[4:6], vec[7:9])$p.value
  #x13 = t.test(vec[1:3], vec[7:9])$p.value
  #c(x12, x23, x13)
  x12
}

MEANs = apply( qPCR, 2, EachMean )
SDs = apply( qPCR, 2, EachSd )
PVs = apply( qPCR, 2, EachPv )

# Draw the bars and error bars
par(las=2, mai=c(2.5,1,0,0))
YMAX = max(MEANs+SDs, na.rm = TRUE)*1.3
x = barplot(MEANs, beside=TRUE, 
            ylim=c(0, YMAX),
            col = grey.colors(2),
            cex.lab=1.5, cex.axis=1.5, cex=1.5)
arrows(x, MEANs+SDs, x, MEANs-SDs,
       code=3, angle=90, length=0.05)

legend("topleft", SampleList, fill=grey.colors(3))

for(n in c(1:N_Sample)) {
  points( rep(x[,n], each = 3), qPCR[,n],
          col = rgb(1, 0, 0, 0.3), pch=19 )
}

# Draw the p-values
for(m in c(1:3)) {
  for(n in c(1:N_Sample)) {
  #
    pv = PVs[m, n]
    if(pv < 0.05) {
      if( m==1 ){
        x1 = x[1,n]
        x2 = x[2,n]
        y_max = max((MEANs+SDs)[,n])
        y1 = y_max * 1.1
        y2 = y_max * 1.12
        lines( c(x1, x1, x2, x2), c(y1, y2, y2, y1) )
        if(pv<=0.05&&pv>0.01) {
          text( (x1+x2)/2, y_max*1.15, "*" )
        } 
        if(pv<=0.01&&pv>0.001) {
          text( (x1+x2)/2, y_max*1.15, "**" )
        } 
        if(pv<=0.001) {
          text( (x1+x2)/2, y_max*1.15, "***" )
        } 
      }
      if( m==2 ){
        x1 = x[2,n]
        x2 = x[3,n]
        y_max = max((MEANs+SDs)[,n])
        y1 = y_max * 1.1
        y2 = y_max * 1.12
        lines( c(x1, x1, x2, x2), c(y1, y2, y2, y1) )
        if(pv<=0.05&&pv>0.01) {
          text( (x1+x2)/2, y_max*1.15, "*" )
        } 
        if(pv<=0.01&&pv>0.001) {
          text( (x1+x2)/2, y_max*1.15, "**" )
        } 
        if(pv<=0.001) {
          text( (x1+x2)/2, y_max*1.15, "***" )
        }       }
      if( m==3 ){
        x1 = x[1,n]
        x2 = x[3,n]
        y_max = max((MEANs+SDs)[,n])
        y1 = y_max * 1.2
        y2 = y_max * 1.22
        lines( c(x1, x1, x2, x2), c(y1, y2, y2, y1) )
        if(pv<=0.05&&pv>0.01) {
          text( (x1+x2)/2, y_max*1.25, "*" )
        } 
        if(pv<=0.01&&pv>0.001) {
          text( (x1+x2)/2, y_max*1.25, "**" )
        } 
        if(pv<=0.001) {
          text( (x1+x2)/2, y_max*1.25, "***" )
        }       }
    }
  #
  } #for n
} #for m

dev.off()
