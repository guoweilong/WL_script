#!/usr/bin/python

import sys
import os
import os.path
import re

#from collections import Counter

BAM_MATCH = 0
BAM_INS = 1
BAM_DEL = 2
BAM_SOFTCLIP = 4

CIGAR_OPS = {'M' : BAM_MATCH, 'I' : BAM_INS, 'D' : BAM_DEL, 'S' : BAM_SOFTCLIP}
def parse_cigar(cigar_string):
    i = 0
    prev_i = 0
    cigar = []
    while i < len(cigar_string):
        if cigar_string[i] in CIGAR_OPS:
            cigar.append((CIGAR_OPS[cigar_string[i]], int(cigar_string[prev_i:i])))
            prev_i = i + 1
        i += 1
    return cigar


"""
SRR299053.104125        0       chr1    174644084       255     72M3I25M        *       0       0       TTTTTTTATTTTTTGTGTGTGTATAAAAGATGGGTTGATTTATTGGTGATTAGCGTTTAATGGAATTTTTTTTTATTTTTTTTAAATCGAAGTGTAAATA  *       XO:Z:+FW      XS:i:0  NM:i:3  XM:Z:--z-zz---z-z----------------------y-------y----------X---z----------z---------------X----------z-        XG:Z:AT_TTCTCCTATCTCTTGTGTGTGTATAAAAGATGGGCTGATTTACTGGTGATTAGCGTTCAATGGAATTTCTTTTTTTTTTTAAATCGAAGTGTAAACA_TA
SRR299053.104431        0       chr1    77429024        255     27M73S  *       0       0       CATTGGAATAGATGTTTTTTTTTTTTTCGTTTTTTTTTGTTGATTTTTTATGTTTTTAGAAGAGTTTCGTTTTTAAGTTTTTTTTCTTGGGTTTGTTTGT  *       XO:Z:+FW     XS:i:0   NM:i:0  XM:Z:Z-y-------------z---z-zzzyx        XG:Z:AC_CACTGGAATAGATGTTCTTTCTCCCCC_GT
"""

from operator import itemgetter
from collections import defaultdict

def process_sam_for_cigar( samfn ):
    try:
        samF = open(samfn, 'r')
    except IOError:
        print "\n[Error]:\n\t File cannot be open: ", samfn
        exit(-1)

    FiveClip = 0
    ThreeClip = 0
    BothClip = 0
    TotalRead = 0
    FiveClip_list = []
    ThreeClip_list = []
    for line in samF:
        TotalRead += 1
        line = line.strip().split()
        cigar = parse_cigar(line[5])
        if line[1] == "0" :
            #print "+", cigar
            if cigar[0][0] == BAM_SOFTCLIP :
                FiveClip += 1
                FiveClip_list.append(cigar[0][1])
            if cigar[-1][0] == BAM_SOFTCLIP :
                ThreeClip +=1
                ThreeClip_list.append(cigar[-1][1])
            if cigar[0][0] == BAM_SOFTCLIP and cigar[-1][0] == BAM_SOFTCLIP :
                BothClip += 1
        elif line[1] == "16" :
            #print "-", cigar
            if cigar[0][0] == BAM_SOFTCLIP :
                ThreeClip += 1
                ThreeClip_list.append(cigar[0][1])
            if cigar[-1][0] == BAM_SOFTCLIP :
                FiveClip +=1
                FiveClip_list.append(cigar[-1][1])
            if cigar[0][0] == BAM_SOFTCLIP and cigar[-1][0] == BAM_SOFTCLIP :
                BothClip += 1
        else :
            print "Bug"

    tolist = lambda A : [(i, A[i]) for i in A]

    print "# of total reads :", TotalRead
    print "# of reads with 5' clipped : ", FiveClip
    print "# of reads with 3' clipped : ", ThreeClip
    print "# of reads with both ends clipped : ", BothClip

    FiveClip_freq = defaultdict(int)
    for item in FiveClip_list:
        FiveClip_freq[item] += 1

    ThreeClip_freq = defaultdict(int)
    for item in ThreeClip_list:
        ThreeClip_freq[item] += 1


    print "5' end clip:\t", sorted(tolist(FiveClip_freq))
    print "3' end clip:\t", sorted(tolist(ThreeClip_freq))

#    print "5' end clip:\t", sorted(tolist(Counter(FiveClip_list)))
#    print "3' end clip:\t", sorted(tolist(Counter(ThreeClip_list)))

    samF.close()



from optparse import OptionParser
# ===========================================
def main():
    usage = "%prog -s <sam_file>\n" \
            "Description : Check the CIGAR section to count 5'/3'-end trimming\n"\
            "Date : 2013-04-24"
    parser = OptionParser(usage)
    parser.add_option("-s", dest="samfile", type="string", help="Name of the input SAM file", metavar="FILE")
    (options, args) = parser.parse_args()

    if (options.samfile is None) :
        print parser.print_help()
        exit(-1)

    process_sam_for_cigar(options.samfile)


# ===========================================
if __name__ == "__main__":
    main()
