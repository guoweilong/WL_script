#!/usr/bin/env python
import os
import gzip
import sys
import re

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)



def get_char_counts( fn ):
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
    for c in IN.read():
        if c not in usage:
            usage.update({c:1})
        else:
            usage[c] += 1

    for i in sorted(usage) :
        if i not in ["\n", " ", "\t"]  :
            print "%c\t%d" % (i, usage[i])

    if fn :
        IN.close()



from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog -i <genome.fa> [-o <output>]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2014-04-23\n" \
            "Description: Count the occurrences of chars (excluding \"\\n\"," \
            "\"\\t\", \" \").\n" \
            "Example:\n" \
            "\t$echo \"guogguogguo\" | %prog\n" \
            "\tg\t5\n" \
            "\to\t3\n" \
            "\tu\t3\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input file", metavar="FILE")
    (options, args) = parser.parse_args()

    get_char_counts (options.infile)


# ===========================================
if __name__ == "__main__":
    main()
