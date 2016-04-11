#!/usr/bin/env python
import os
import gzip
import sys
import re

def error(msg):
    print >> sys.stderr, 'ERROR: %s' % msg
    exit(1)

from optparse import OptionParser

# RefFlat format
# LINC00336	NR_027908	chr6	-	33661860	33669093	33669093	33669093	2	33661860,33668751,	33663625,33669093,

# GTF format
"""
chr1	mm9_refGene	start_codon	134212807	134212809	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134212807	134213049	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134212702	134213049	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134221530	134221650	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134221530	134221650	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134224274	134224425	0.000000	+	2	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134224274	134224425	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134224708	134224773	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134224708	134224773	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134226535	134226654	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
guoweilong@Downloads$ zcat mm9_refseq.gtf.gz | head -20
chr1	mm9_refGene	start_codon	134212807	134212809	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134212807	134213049	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134212702	134213049	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134221530	134221650	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134221530	134221650	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134224274	134224425	0.000000	+	2	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134224274	134224425	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134224708	134224773	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134224708	134224773	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134226535	134226654	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134226535	134226654	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134227136	134227268	0.000000	+	0	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134227136	134227268	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	CDS	134227898	134228955	0.000000	+	2	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	stop_codon	134228956	134228958	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	exon	134227898	134230065	0.000000	+	.	gene_id "NM_028778"; transcript_id "NM_028778";
chr1	mm9_refGene	start_codon	134212807	134212809	0.000000	+	.	gene_id "NM_001195025"; transcript_id "NM_001195025";
chr1	mm9_refGene	CDS	134212807	134213049	0.000000	+	0	gene_id "NM_001195025"; transcript_id "NM_001195025";
chr1	mm9_refGene	exon	134212702	134213049	0.000000	+	.	gene_id "NM_001195025"; transcript_id "NM_001195025";
chr1	mm9_refGene	CDS	134221530	134221650	0.000000	+	0	gene_id "NM_001195025"; transcript_id "NM_001195025";

"""

"""
url: http://mblab.wustl.edu/GTF2.html

Fields

Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

    seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix.
    source - name of the program that generated this feature, or the data source (database or project name)
    feature - feature type name, e.g. Gene, Variation, Similarity
    start - Start position of the feature, with sequence numbering starting at 1.
    end - End position of the feature, with sequence numbering starting at 1.
    score - A floating point value.
    strand - defined as + (forward) or - (reverse).
    frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
    attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.

"""


# ===========================================
def main():
    usage = "Usage: %prog -r <refFlat>\n" \
            "Author : Guo, Weilong; guoweilong@126.com; 2014-06-24\n" \
            "Description: Generate GTF file from RefFlat file.\n" \
            "Note: The generated GTF format can be used in Cufflinks."
    parser = OptionParser(usage)
    parser.add_option("-r", dest="refFlat",
                  help="Input file, refFlat format. If not specified, "
                       "STDIN will be used.", metavar="FILE", default=None)
    parser.add_option("-o", dest="outfile",
                  help="Output file, GTF format", metavar="FILE", default=None)
    (options, args) = parser.parse_args()

    if options.refFlat :
        try :
            IN = open(options.refFlat)
        except IOError:
            print "[Error] Cannot find input file : %s !" % options.refFlat
            exit(-1)
    else :
        IN = sys.stdin

    if (options.outfile is not None) :
        if options.outfile.endswith('.gz') :
            sys.stdout = gzip.open(options.outfile, 'wb')
        else :
            sys.stdout = open(options.outfile, 'w')

    source = "RefFlat"
    score = 0
    frame = "."

    #LINC00336       NR_027908       chr6    -       33661860        33669093        33669093        33669093   33661860,33668751,       33663625,33669093,
    GTF_list = []
    for line in IN :
        #print line.strip().split()
        [ Gene_ID, Trans_ID, chr, strand, Left, Right, LeftCodon, RightCodon, _, ExonLeft, ExonRight]=line.strip().split()
        LeftCodon = int(LeftCodon)
        RightCodon = int(RightCodon)
        # GTF format
        # seqname, source, feature, start, end, score, strand, frame, attribute
        attribute ="gene_id \"" + Gene_ID + "\"; transcript_id \"" + Trans_ID + "\""
        if strand == "+" :
            for [exon_left, exon_right] in zip(ExonLeft.strip(",").split(","), ExonRight.strip(",").split(",")) :
                # exon
                GTF_list.append( [chr, source, "exon", exon_left, exon_right, score, strand, frame, attribute] )
                # CDS
                CDS_left = exon_left if exon_left > LeftCodon else LeftCodon
                CDS_right = exon_right if exon_right < RightCodon else RightCodon
                if CDS_left < CDS_right :
                    GTF_list.append( [chr, source, "CDS", CDS_left, CDS_right, score, strand, frame, attribute] )
            if LeftCodon < RightCodon :
                # start_codon
                GTF_list.append( [chr, source, "start_codon", LeftCodon, LeftCodon+2, score, strand, frame, attribute] )
                # stop_codon
                GTF_list.append( [chr, source, "stop_codon", RightCodon-2, RightCodon, score, strand, frame, attribute] )

        else :
            for [exon_left, exon_right] in zip(ExonLeft.strip(",").split(","), ExonRight.strip(",").split(",")) :
                # exon
                GTF_list.append( [chr, source, "exon", exon_left, exon_right, score, strand, frame, attribute] )
                # CDS
                CDS_left = exon_left if exon_left > LeftCodon else LeftCodon
                CDS_right = exon_right if exon_right < RightCodon else RightCodon
                if CDS_left < CDS_right :
                    GTF_list.append( [chr, source, "CDS", CDS_left, CDS_right, score, strand, frame, attribute] )
            if LeftCodon < RightCodon :
                # stop_codon
                GTF_list.append( [chr, source, "stop_codon", LeftCodon, LeftCodon+2, score, strand, frame, attribute] )
                # start_codon
                GTF_list.append( [chr, source, "start_codon", RightCodon-2, RightCodon, score, strand, frame, attribute] )

    # End of "for line in INPUT"
    IN.close()

    GTF_list.sort(key=lambda ele : [ele[0], int(ele[3]), int(ele[4]), ele[2] ])
    for ele in GTF_list :
        print "\t".join( str(i) for i in ele )



# ===========================================
if __name__ == "__main__":
    main()


# todo
# the <frame> field of GTF files was not fixed.
"""
Here is an example of a gene on the negative strand. Larger coordinates are 5' of smaller coordinates. Thus, the start codon is 3 bp with largest coordinates among all those bp that fall within the CDS regions. Similarly, the stop codon is the 3 bp with coordinates just less than the smallest coordinates within the CDS regions.

AB000123    Twinscan     CDS    193817    194022    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    199645    199752    .    -    2    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    200369    200508    .    -    1    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     CDS    215991    216028    .    -    0    gene_id "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     start_codon   216026    216028    .    -    .    gene_id    "AB000123.1"; transcript_id "AB00123.1.2";
AB000123    Twinscan     stop_codon    193814    193816    .    -    .    gene_id    "AB000123.1"; transcript_id "AB00123.1.2";

Note the frames of the coding exons. For example:

    The first CDS (from 216028 to 215991) always has frame zero.
    Frame of the 1st CDS =0, length =38.  (frame - length) % 3  = 1, the frame of the 2nd CDS.
    Frame of the 2nd CDS=1, length=140. (frame - length) % 3  = 2, the frame of the 3rd CDS.
    Frame of the 3rd CDS=2, length=108. (frame - length) % 3  =  2, the frame of the terminal CDS.
    Alternatively, the frame of terminal CDS can be calculated without the rest of the gene. Length of the terminal CDS=206. length % 3 =2, the frame of the terminal CDS.
"""


