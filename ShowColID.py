#!/usr/bin/env python
import os
import gzip
import sys
import re

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)



def ShowColID( fn, sep="\t" ):
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

    usage = {}
    line = IN.readline();
    tokens = line.strip().split(sep)
    for i in xrange( len(tokens) ) :
        print "%d\t%s" % (i+1, tokens[i])
    if fn :
        IN.close()


from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog -i <input.txt> [-o <output>]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2015-05-19\n" \
            "Description: Show the columns of 1st line in input by column numbers.\n" \
            "Example:\n" \
            "Input file:\n" \
            "  Wrnip1	NM_030215	chr13\n" \
            "Output:\n" \
            "  1	Wrnip1\n" \
            "  2	NM_030215\n" \
            "  3	chr13"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file", metavar="FILE")
    (options, args) = parser.parse_args()

    ShowColID (options.infile)


# ===========================================
if __name__ == "__main__":
    main()
