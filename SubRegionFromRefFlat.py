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
            "Last Update: 2013-09-19\n" \
            "Description: Generate sub-region in BED format from RefFlat file.\n"
    parser = OptionParser(usage)
    parser.add_option("-i", dest="infile",
                  help="Genome sequence file in Fasta format", metavar="FILE")
    parser.add_option("-E", dest="Exon", help="File name for exon\n", metavar="FILE", default=None)
    parser.add_option("-I", dest="Intron", help="File name for intron\n", metavar="FILE", default=None)
    parser.add_option("-F", dest="FiveUTR", help="File name for 5\'UTR\n", metavar="FILE", default=None)
    parser.add_option("-T", dest="ThreeUTR", help="File name for 3\'UTR\n", metavar="FILE", default=None)
    parser.add_option("-O", dest="ORF", help="File name for ORF\n", metavar="FILE", default=None)
    parser.add_option("-G", dest="Genebody", help="File name for genebody\n", metavar="FILE", default=None)
    parser.add_option("-P", dest="Promoter", help="File name for promoter\n", metavar="FILE", default=None)
    parser.add_option("--up", dest="upstream", help="length of upstream for promoter [%default]\n", metavar="INT", default=500)
    parser.add_option("--down", dest="downstream", help="length of downstream for promoter [%default]\n", metavar="INT", default=100)
    (options, args) = parser.parse_args()

    if (options.infile is None) :
        print parser.print_help()
        exit(-1)

    try :
        INPUT =  open(options.infile)
    except IOError:
        print "[Error] Cannot find input file : %s !" % options.infile
        exit(-1)

    Exon_lst=[]
    Intron_lst=[]
    FiveUTR_lst=[]
    ThreeUTR_lst=[]
    ORF_lst=[]
    Genebody_lst=[]
    Promoter_lst=[]

    #LINC00336       NR_027908       chr6    -       33661860        33669093        33669093        33669093   33661860,33668751,       33663625,33669093,
    for line in INPUT :
        #print line.strip().split()
        [_, _, chr, strand, Left, Right, LeftCoding, RightCoding, _, ExonLeft, ExonRight]=line.strip().split()
        ExonLefts=ExonLeft[:-1].split(",")
        ExonRights=ExonRight[:-1].split(",")
        # Exon
        for i in range(len(ExonLefts)) :
            Exon_lst.append( "\t".join([chr, ExonLefts[i], ExonRights[i], strand ]) )

        # Intron
        for i in range(len(ExonLefts)-1) :
            Intron_lst.append( "\t".join( [chr, str( int(ExonRights[i])+1), str(int(ExonLefts[i+1])-1), strand ] ) )

        # FiveUTR
        if strand == "+" :
            FiveUTR_lst.append( "\t".join( [chr, Left, LeftCoding, strand ] ) )
        else :
            FiveUTR_lst.append( "\t".join( [chr, RightCoding, Right, strand ] ) )

        # ThreeUTR
        if strand == "+" :
            ThreeUTR_lst.append( "\t".join( [chr, RightCoding, Right, strand ] ) )
        else :
            ThreeUTR_lst.append( "\t".join( [chr, Left, LeftCoding, strand ] ) )

        # ORF
        ORF_lst.append( "\t".join( [chr, LeftCoding, RightCoding, strand ] ) )

        # Genebody
        Genebody_lst.append( "\t".join( [chr, Left, Right, strand ] ) )

        # Promoter
        if strand=="+" :
            if int(Left) - options.upstream > 0 :
                Promoter_lst.append( "\t".join( [chr, str(int(Left)-options.upstream),
                                               str(int(Left)+options.downstream), strand] ) )
            else :
                Promoter_lst.append( "\t".join( [chr, "0", str(int(Left)+options.downstream), strand] ) )

        else :
            if int(Right) - options.downstream > 0 :
                Promoter_lst.append( "\t".join( [chr, str(int(Right)-options.downstream),
                                               str(int(Right)+options.upstream), strand] ) )
            else :
                Promoter_lst.append( "\t".join( [chr, "0", str(int(Right)+options.upstream), strand] ) )

    INPUT.close()

    for [fn, lst]  in [[options.Exon, Exon_lst],
                       [options.Intron, Intron_lst],
                       [options.FiveUTR, FiveUTR_lst],
                       [options.ThreeUTR, ThreeUTR_lst],
                       [options.ORF, ORF_lst],
                       [options.Genebody, Genebody_lst],
                       [options.Promoter, Promoter_lst] ] :
        if fn is not None :
            OF = open(fn, 'w')
            for i in sorted( set(lst) ) :
                [chr, start, end, strand] = i.split()
                if start < end :
                    OF.writelines(i+"\n")
            OF.close()


# ===========================================
if __name__ == "__main__":
    main()
