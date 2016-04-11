#!/usr/bin/env python
import os
import gzip
import sys
#import re
import numpy as np

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)



def Quantile( fn, col = 1, pct = 0.5 ):
    try:
        if fn :
            if fn.endswith(".gz") :
                IN = gzip.open(fn, 'rb')
            else :
                IN = open(fn, 'r')
        else :
            IN = sys.stdin
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", fn
        exit(-1)
    #
    col = int(col)-1
    pct = float(pct)
    PCT = pct*100
    if PCT < 0 :
        PCT = 0
    elif PCT >100 :
        PCT =100
    #
    LST = []
    for line in IN :
        line = line.strip()
        tokens = line.split()
        if col >= len(tokens) or col < 0 :
            sys.stderr.write("line with wrong number of column\n  \"%s\"\n" % line)
        else :
            LST.append( float(tokens[col]) )
    #
    if LST != [] :
        print "%.2f" % np.percentile(LST, PCT)
    if fn :
        IN.close()


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog [-i <input.txt>] -c 1 -p 0.5 \n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2015-09-30\n" \
            "Description: Get the quantile number for certain column\n" \
            "Example:\n" \
            "Input file:\n" \
            "  1	0.1\n" \
            "  2	0.5\n" \
            "  3	0.2\n"\
            "Output:\n" \
            "   2 (by default)"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file, or STDIN if omitted", metavar="FILE")
    parser.add_option("-c", dest="col",
                  help="The column in input file [default: %default]", default = 1, metavar="INT")
    parser.add_option("-p", dest="pct",
                  help="the percentage, range from 0~1 [default: %default]", default = 0.5, metavar="double")
    (options, args) = parser.parse_args()

    Quantile (options.infile, options.col, options.pct)


# ===========================================
if __name__ == "__main__":
    main()
