#!/usr/bin/env python
import os
import gzip
import sys
import re

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)

from optparse import OptionParser

# ===========================================
def main():
    usage = "Usage: %prog -i <input> [-E -I -F -T -O -G -P --up=500 --down=100]\n" \
            "Author : Guo, Weilong; guoweilong@gmail.com; 2013-09-19\n" \
            "Last Update: 2015-11-22\n" \
            "Description: Generate sub-region in BED format from RefFlat file.\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Input RefFlat", metavar="FILE")
    parser.add_option("-E", dest="Exon", help="File name for exon\n", metavar="FILE", default=None)
    parser.add_option("--1stExon", dest="FirstExon", help="File name for 1st exon\n", metavar="FILE", default=None)
    parser.add_option("--posteriorExon", dest="PosteriorExon", help="File name for posterior exons, 2nd to last\n", metavar="FILE", default=None)
    parser.add_option("-I", dest="Intron", help="File name for intron\n", metavar="FILE", default=None)
    parser.add_option("--1stIntron", dest="FirstIntron", help="File name for 1st intron\n", metavar="FILE", default=None)
    parser.add_option("--posteriorIntron", dest="PosteriorIntron", help="File name for posterior introns, 2nd to last\n", metavar="FILE", default=None)
    parser.add_option("-L", "--LeftIntron", dest="LeftIntron", help="File name for (5'SS) left intron\n", metavar="FILE", default=None)
    parser.add_option("-M", "--MiddleIntron", dest="MiddleIntron", help="File name for (MI) middle intron\n", metavar="FILE", default=None)
    parser.add_option("-R", "--RightIntron", dest="RightIntron", help="File name for (3'SS) right intron\n", metavar="FILE", default=None)

    parser.add_option("--1stLeftIntron", dest="FirstLeftIntron", help="File name for 1st (5'SS) left intron\n", metavar="FILE", default=None)
    parser.add_option("--1stMiddleIntron", dest="FirstMiddleIntron", help="File name for 1st (MI) middle intron\n", metavar="FILE", default=None)
    parser.add_option("--1stRightIntron", dest="FirstRightIntron", help="File name for 1st (3'SS) right intron\n", metavar="FILE", default=None)

    parser.add_option("--posteriorLeftIntron", dest="PosteriorLeftIntron", help="File name for posterior (5'SS) left intron\n", metavar="FILE", default=None)
    parser.add_option("--posteriorMiddleIntron", dest="PosteriorMiddleIntron", help="File name for posterior (MI) middle intron\n", metavar="FILE", default=None)
    parser.add_option("--posteriorRightIntron", dest="PosteriorRightIntron", help="File name for posterior (3'SS) right intron\n", metavar="FILE", default=None)

    parser.add_option("-F", dest="FiveUTR", help="File name for 5\'UTR\n", metavar="FILE", default=None)
    parser.add_option("-T", dest="ThreeUTR", help="File name for 3\'UTR\n", metavar="FILE", default=None)
    parser.add_option("-O", dest="ORF", help="File name for ORF\n", metavar="FILE", default=None)
    parser.add_option("-G", dest="Genebody", help="File name for genebody\n", metavar="FILE", default=None)
    parser.add_option("-P", dest="Promoter", help="File name for promoter\n", metavar="FILE", default=None)
    parser.add_option("--up", dest="upstream", help="length of upstream for promoter [%default]\n", metavar="INT", default=500)
    parser.add_option("--down", dest="downstream", help="length of downstream for promoter [%default]\n", metavar="INT", default=100)
    parser.add_option("--Intergenic", dest="Intergenic", help="File name for Intergenic region\n", metavar="FILE", default=None)
    parser.add_option("--IntergenicFrom", dest="IntergenicFrom", help="From the position upstream of TSS [%default]\n", metavar="INT", default=10000)
    parser.add_option("--IntergenicTo", dest="IntergenicTo", help="To the position upstream of TSS [%default]\n", metavar="INT", default=9000)
    (options, args) = parser.parse_args()

    if (options.infile is None) :
    #    print parser.print_help()
    #    exit(-1)
        INPUT = sys.stdin
    else :
        try :
            INPUT =  open(options.infile)
        except IOError:
            print "[Error] Cannot find input file : %s !" % options.infile
            exit(-1)
    #
    Exon_lst=[]
    FirstExon_lst=[]
    PosteriorExon_lst=[]
    Intron_lst=[]
    FirstIntron_lst=[]
    PosteriorIntron_lst=[]
    LeftIntron_lst=[]
    MiddleIntron_lst=[]
    RightIntron_lst=[]
    #
    FirstLeftIntron_lst=[]
    FirstMiddleIntron_lst=[]
    FirstRightIntron_lst=[]
    #
    PosteriorLeftIntron_lst=[]
    PosteriorMiddleIntron_lst=[]
    PosteriorRightIntron_lst=[]
    #
    FiveUTR_lst=[]
    ThreeUTR_lst=[]
    ORF_lst=[]
    Genebody_lst=[]
    #
    Promoter_lst=[]
    upstream = int(options.upstream)
    downstream = int(options.downstream)
    #
    Intergenic_lst=[]
    IntergenicFrom = int(options.IntergenicFrom)
    IntergenicTo   = int(options.IntergenicTo)
    #
    #LINC00336       NR_027908       chr6    -       33661860        33669093        33669093        33669093   33661860,33668751,       33663625,33669093,
    for line in INPUT :
        #print line.strip().split()
        [_, _, chr, strand, Left, Right, LeftCoding, RightCoding, _, ExonLeft, ExonRight]=line.strip().split()
        ExonLefts=ExonLeft[:-1].split(",")
        ExonRights=ExonRight[:-1].split(",")
        # Exon
        for i in xrange(len(ExonLefts)) :
            Exon_lst.append( "\t".join([chr, ExonLefts[i], ExonRights[i], strand ]) )
        #
        # First Exon
        if strand == "+" :
            FirstExon_lst.append( "\t".join([chr, ExonLefts[0], ExonRights[0], strand ]) )
        else :
            FirstExon_lst.append( "\t".join([chr, ExonLefts[-1], ExonRights[-1], strand ]) )
        #
        # Posterior Exons
        if strand == "+" :
            for i in xrange(1, len(ExonLefts)) :
                PosteriorExon_lst.append( "\t".join([chr, ExonLefts[i], ExonRights[i], strand ]) )
        else :
            for i in xrange(0, len(ExonLefts)-1) :
                PosteriorExon_lst.append( "\t".join([chr, ExonLefts[i], ExonRights[i], strand ]) )
        #
        # Intron
        for i in xrange(len(ExonLefts)-1) :
            Intron_lst.append( "\t".join( [chr, str( int(ExonRights[i])+1), str(int(ExonLefts[i+1])-1), strand ] ) )
        #
        # First Intron
        if len(ExonLefts) > 1 :
            if strand == "+" :
                FirstIntron_lst.append( "\t".join( [chr, str( int(ExonRights[0])+1), str(int(ExonLefts[1])-1), strand ] ) )
            else :
                FirstIntron_lst.append( "\t".join( [chr, str( int(ExonRights[-2])+1), str(int(ExonLefts[-1])-1), strand ] ) )
        #
        # Posterior Intron
        if len(ExonLefts) > 2 :
            if strand == "+" :
                for i in xrange(1, len(ExonLefts)-1) :
                    PosteriorIntron_lst.append( "\t".join( [chr, str( int(ExonRights[i])+1), str(int(ExonLefts[i+1])-1), strand ] ) )
            else :
                for i in xrange(0, len(ExonLefts)-2) :
                    PosteriorIntron_lst.append( "\t".join( [chr, str( int(ExonRights[i])+1), str(int(ExonLefts[i+1])-1), strand ] ) )
        #
        # LeftIntron
        # MiddleIntron
        # RightIntron
        if strand == "+" :
            for i in xrange(len(ExonLefts)-1) :
                LeftPos = int(ExonRights[i])
                RightPos = int(ExonLefts[i+1])
                if  (RightPos - LeftPos) >= 500 :
                    LeftIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                    MiddlePos = (LeftPos + RightPos) / 2
                    MiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                    RightIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        else :
            for i in xrange(len(ExonLefts)-1) :
                LeftPos = int(ExonRights[i])
                RightPos = int(ExonLefts[i+1])
                if  (RightPos - LeftPos) >= 500 :
                    RightIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                    MiddlePos = (LeftPos + RightPos) / 2
                    MiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                    LeftIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        #
        # FirstLeftIntron
        # FirstMiddleIntron
        # FirstRightIntron
        if strand == "+" :
            if len(ExonLefts) > 1 :
                LeftPos = int(ExonRights[0])
                RightPos = int(ExonLefts[1])
                if  (RightPos - LeftPos) >= 500 :
                    FirstLeftIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                    MiddlePos = (LeftPos + RightPos) / 2
                    FirstMiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                    FirstRightIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        else :
            if len(ExonLefts) > 1 :
                LeftPos = int(ExonRights[-2])
                RightPos = int(ExonLefts[-1])
                if  (RightPos - LeftPos) >= 500 :
                    FirstRightIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                    MiddlePos = (LeftPos + RightPos) / 2
                    FirstMiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                    FirstLeftIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        #
        # PosteriorLeftIntron
        # PosteriorMiddleIntron
        # PosteriorRightIntron
        if strand == "+" :
            if len(ExonLefts) > 2 :
                for i in xrange(1, len(ExonLefts)-1) :
                    LeftPos = int(ExonRights[i])
                    RightPos = int(ExonLefts[i+1])
                    if  (RightPos - LeftPos) >= 500 :
                        PosteriorLeftIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                        MiddlePos = (LeftPos + RightPos) / 2
                        PosteriorMiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                        PosteriorRightIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        else :
            if len(ExonLefts) > 2 :
                for i in xrange(0, len(ExonLefts)-2) :
                    LeftPos = int(ExonRights[i])
                    RightPos = int(ExonLefts[i+1])
                    if  (RightPos - LeftPos) >= 500 :
                        PosteriorRightIntron_lst.append( "\t".join( [chr, str(LeftPos+5), str(LeftPos+125), strand ] ) )
                        MiddlePos = (LeftPos + RightPos) / 2
                        PosteriorMiddleIntron_lst.append( "\t".join( [chr, str(MiddlePos-60), str(MiddlePos+60), strand ] ) )
                        PosteriorLeftIntron_lst.append( "\t".join( [chr, str(RightPos-125), str(RightPos-5), strand ] ) )
        #
        # FiveUTR
        if strand == "+" :
            FiveUTR_lst.append( "\t".join( [chr, Left, LeftCoding, strand ] ) )
        else :
            FiveUTR_lst.append( "\t".join( [chr, RightCoding, Right, strand ] ) )
        #
        # ThreeUTR
        if strand == "+" :
            ThreeUTR_lst.append( "\t".join( [chr, RightCoding, Right, strand ] ) )
        else :
            ThreeUTR_lst.append( "\t".join( [chr, Left, LeftCoding, strand ] ) )
        #
        # ORF
        ORF_lst.append( "\t".join( [chr, LeftCoding, RightCoding, strand ] ) )
        #
        # Genebody
        Genebody_lst.append( "\t".join( [chr, Left, Right, strand ] ) )
        #
        # Promoter
        if strand=="+" :
            if int(Left) - upstream > 0 :
                Promoter_lst.append( "\t".join( [chr, str(int(Left)-upstream),
                                               str(int(Left)+downstream), strand] ) )
            else :
                Promoter_lst.append( "\t".join( [chr, "0", str(int(Left)+downstream), strand] ) )
        else :
            if int(Right) - downstream > 0 :
                Promoter_lst.append( "\t".join( [chr, str(int(Right)-downstream),
                                               str(int(Right)+upstream), strand] ) )
            else :
                Promoter_lst.append( "\t".join( [chr, "0", str(int(Right)+upstream), strand] ) )
        #
        # Intergenic
        if strand=="+" :
            if int(Left) - IntergenicFrom > 0 :
                Intergenic_lst.append( "\t".join( [chr, str(int(Left)-IntergenicFrom),
                                               str(int(Left)-IntergenicTo), strand] ) )
            else :
                Intergenic_lst.append( "\t".join( [chr, "0", str(int(Left)-IntergenicTo), strand] ) )
        else :
            if int(Right) + IntergenicTo > 0 :
                Intergenic_lst.append( "\t".join( [chr, str(int(Right)+IntergenicTo),
                                               str(int(Right)+IntergenicFrom), strand] ) )
            else :
                Intergenic_lst.append( "\t".join( [chr, "0", str(int(Right)+IntergenicFrom), strand] ) )
        #
    #
    if options.infile is not None :
        INPUT.close()
    #
    for [fn, lst]  in [[options.Exon, Exon_lst],
                       [options.FirstExon, FirstExon_lst],
                       [options.PosteriorExon, PosteriorExon_lst],
                       [options.Intron, Intron_lst],
                       [options.FirstIntron, FirstIntron_lst],
                       [options.PosteriorIntron, PosteriorIntron_lst],
                       [options.LeftIntron, LeftIntron_lst],
                       [options.MiddleIntron, MiddleIntron_lst],
                       [options.RightIntron, RightIntron_lst],
                       [options.FirstLeftIntron, FirstLeftIntron_lst],
                       [options.FirstMiddleIntron, FirstMiddleIntron_lst],
                       [options.FirstRightIntron, FirstRightIntron_lst],
                       [options.PosteriorLeftIntron, PosteriorLeftIntron_lst],
                       [options.PosteriorMiddleIntron, PosteriorMiddleIntron_lst],
                       [options.PosteriorRightIntron, PosteriorRightIntron_lst],
                       [options.FiveUTR, FiveUTR_lst],
                       [options.ThreeUTR, ThreeUTR_lst],
                       [options.ORF, ORF_lst],
                       [options.Genebody, Genebody_lst],
                       [options.Promoter, Promoter_lst],
                       [options.Intergenic, Intergenic_lst]] :
        if fn is not None :
            OF = open(fn, 'w')
            for i in sorted( set(lst) ) :
                [chr, start, end, strand] = i.split()
                if start < end :
                    OF.writelines(i+"\n")
            OF.close()
        #
    #
#
# ===========================================
if __name__ == "__main__":
    main()


"""
Command for Testing:
../SubRegionFromRefFlat.py -i test.refFlat -E exon.tmp --1stExon 1stExon.tmp --posteriorExon=posteriorExon.tmp \
                           -I intron.tmp --1stIntron=1stIntron.tmp --posteriorIntron=posteriorIntron.tmp \
                           -L leftIntron.tmp -M middleIntron.tmp -R rightIntron.tmp --1stLeftIntron=1stLeftIntron.tmp \
                           --1stMiddleIntron=1stMiddleIntron.tmp --1stRightIntron=1stRightIntron.tmp \
                           --posteriorLeftIntron=posteriorLeftIntron.tmp \
                           --posteriorMiddleIntron=posteriorMiddleIntron.tmp \
                           --posteriorRightIntron=posteriorRightIntron.tmp -F fiveUTR.tmp -T threeUTR.tmp -O orf.tmp \
                           -G genebody.tmp -P promoter.tmp

../SubRegionFromRefFlat.py -i test.refFlat -E exon.tmp --1stExon 1stExon.tmp --posteriorExon=posteriorExon.tmp \
                           -I intron.tmp --1stIntron=1stIntron.tmp --posteriorIntron=posteriorIntron.tmp \
                           -L leftIntron.tmp -M middleIntron.tmp -R rightIntron.tmp --1stLeftIntron=1stLeftIntron.tmp \
                           --1stMiddleIntron=1stMiddleIntron.tmp --1stRightIntron=1stRightIntron.tmp \
                           --posteriorLeftIntron=posteriorLeftIntron.tmp \
                           --posteriorMiddleIntron=posteriorMiddleIntron.tmp \
                           --posteriorRightIntron=posteriorRightIntron.tmp -F fiveUTR.tmp -T threeUTR.tmp -O orf.tmp \
                           -G genebody.tmp -P promoter.tmp

"""

