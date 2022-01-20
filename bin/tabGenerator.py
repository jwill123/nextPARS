#!/usr/bin/python

import os
import sys
import os.path
import argparse

## this is just wrapper script
## create TAB files from proifle fiels

#~ /users/tg/jrodriguez/lncRNA_PROJECT/__bin_2013_10_24/RNA_Get_PARS_Position_Exons.py
#~ -i /users/tg/jrodriguez/lncRNA_PROJECT/PARS_POSITIONS/OUTPUTS_ALL_2013_10_28/tophat2/Dig33_V1/Dig33_V1.bam 
#~ -r /users/tg/jrodriguez/lncRNA_PROJECT/__RAW_DATA/__REFERENCE/__PAPER/HOTAIR2.bed 
#~ -o /users/tg/jrodriguez/lncRNA_PROJECT/PARS_POSITIONS/OUTPUTS_ALL_2013_10_28/tophat2/Dig33_V1/pars/HOTAIR2_Dig33_V1.tab

#bed = "/home/dloska/Desktop/RNA/tab/tabs/mRNAsaccharomycsControls.bed" # Saccharomyces cerevisiae  * Controls: TETp4, TETp9, HOT2, SRA, B2, U1  
#~ bed = "/home/dloska/Desktop/RNA/tab/tabs/test_mRNA.bed"
#~ bed = "/home/dloska/Desktop/RNA/tab/tabs/ester-trimmed.bed"
#~ bed = "/home/dloska/Desktop/RNA/tab/tabs/reference-Hall.bed"


#~ suffix = "X"
#~ suffix = "NoDig"
#~ suffix = "S1"
#~ suffix = "V1"

#home = '/home/jwillis/users/tg'
#home = '/users/tg'
#src  = '%s/jwillis/lncRNA_DATA/RNA2D' %home

#~ outDir = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/candida/"
#outDir = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/"
## "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/refR61_bedR61/"
#outDir = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data' %home



home = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
data = '%s/data' %home
src = '%s/bin' %home
outDir = '%s/tabGenerator_outputs' %data
########################################################################
########################################################################
def repeatIt(bamfolder, name, bed, ref, date, mapQ, minC, v_lev):

    for dirpath, dirnames, filenames in os.walk(bamfolder):
        #~ for filename in [f for f in filenames if  f.endswith(".bam") and f.startswith("D") ]:
        for filename in [f for f in filenames if  f.endswith(".bam")  ]:
            if "unmapped.bam" in filename: continue
            if "accepted_hits.bam" in filename: continue
            if "accepted_hits.sorted.bam" in filename: continue
            if "sorted_reads.bam" in filename: continue
            if "tmp" in dirpath: continue

            #command = "python RNA_Get_PARS_Position_Exons.py -i " + os.path.join(dirpath, filename) + " -r " + bed +" -o /home/dloska/Desktop/RNA/tab/tabs/"+ name +".tab" ;
            #~ print os.path.join(dirpath, filename)
            #print os.path.join(command)
            #os.system(command)

            if "_A_" in bamfolder:
                suffix = "A"
            elif "control" in bamfolder:
                suffix = "NoDig"
            elif "Nodig" in bamfolder:
                suffix = "NoDig"
            elif "NoDigest" in bamfolder:
                suffix = "NoDig"
            elif "NO_dig" in bamfolder:
                suffix = "NoDig"
            elif "_S1" in bamfolder or "S1_" in bamfolder:
                suffix = "S1"
            elif "_V1" in bamfolder or "V1_" in bamfolder:
                suffix = "V1"
            elif "_T1" in bamfolder:
                suffix = "T1"
            elif "_MIX" in bamfolder:
                suffix = "MIX"
            
            if not os.path.exists("%s/%s" %(outDir, ref)):
                os.makedirs("%s/%s" %(outDir, ref))
            if not os.path.exists("%s/%s/%s" %(outDir, ref, date)):
                os.makedirs("%s/%s/%s" %(outDir, ref, date))
        
            t2c = "%s/transcript2count.py" %src
#            if filename.endswith('Aligned.sortedByCoord.out.bam'):
            command = "python %s -v %s -a %s -b %s -q %s -m %s | cut -f1,3 | sed \"s/$/;/g\" >  %s/%s/%s/%s_%s.tab" %(t2c, v_lev, 
                                                                                                                os.path.join(dirpath, filename),
                                                                                                                bed, mapQ, minC, outDir, ref, date, 
                                                                                                                name, suffix)
#            else:
#                command = "python %s -v %s -a %s -b %s -q %s | cut -f1,3 | sed \"s/$/;/g\" >  %s/%s/%s/%s_%s.tab" %(t2c, v_lev, 
#                                                                                                                    os.path.join(dirpath, filename), 
#                                                                                                                    bed, mapQ, outDir, ref, date, 
#                                                                                                                    name, suffix)
            #~ command = "python transcript2count.py -v -a " + os.path.join(dirpath, filename) + " -b " + bed + " | cut -f1,3 | sed \"s/$/;/g\" >  /home/dloska/Desktop/RNA/tab/tabs/" + name +".tab" ; 
#            print command

            os.system(command)

#python /users/tg/lpryszcz/cluster/pars/genome/tophat2/src/transcript2count.py -v -a /users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2A_S1_5437_TGACCA/H2A_S1_5437_TGACCA.bam -b /home/dloska/Desktop/RNA/tab/tabs/reference-Hall.bed -q 20
#src/transcript2count.py -v -a $f -b paper/mRNA.bed -q $q | gzip > $f.q$q.counts.gz
#skrypt jest w /users/tg/lpryszcz/cluster/pars/genome/tophat2/src
#zerknij tez w /users/tg/lpryszcz/cluster/pars/genome/tophat2/log.txt
# zgrep YAL044C /users/tg/lpryszcz/cluster/pars/genome/tophat2/Dig26_S1.bam.q0.counts.gz
			
               #~
#    if not os.path.exists("%s/%s/%s" %(outDir, ref, name)):
#        os.makedirs("%s/%s/%s" %(outDir, ref, name))

    # read the temporary file with list of 
    fin = open("%s/%s/%s/%s_%s.tab" %(outDir, ref, date, name, suffix), 'r')
    for line in fin:

        line2 = line.strip()

        molecule = line2.split("\t")

        if "_A_" in bamfolder:
            suffix = "A"
        elif "control" in bamfolder:
            suffix = "NoDig"
        elif "Nodig" in bamfolder:
            suffix = "NoDig"
        elif "NoDigest" in bamfolder:
            suffix = "NoDig"
        elif "NO_dig" in bamfolder:
            suffix = "NoDig"
        elif "_S1" in bamfolder or "S1_" in bamfolder:
            suffix = "S1"
        elif "_V1" in bamfolder or "V1_" in bamfolder:
            suffix = "V1"
        elif "_T1" in bamfolder:
            suffix = "T1"
        elif "_MIX" in bamfolder:
            suffix = "MIX"
#        else:
#            continue




        out = open("%s/%s/%s/%s_%s_%s.tab" %(outDir, ref, date, molecule[0], name, suffix), 'w')
        #~ out = open("/home/dloska/Desktop/RNA/tab/tabsSingle/"  + arr[0] + "_" + name + "_" + suffix + ".tab", 'w')
        out.write(line2 + "\n")
        out.close()
    fin.close()
    
    
#    os.system('rm %s/%s/%s.tab' %(outDir, ref, name))
    print
########################################################################
########################################################################
########################################################################

parser = argparse.ArgumentParser(description="Produce tab files from bam inputs")
parser.add_argument("-e", "--experiment", dest="exp", action="store", required=True, help="Experiment name")
parser.add_argument("-b", "--bed", dest="bedfile", action="store", required=True, help="Bedfile with chromosome, position of gene (start-end) and the name of the gene")
parser.add_argument("--bam", dest="bamfolder", action="store", required=True, help="Folder containing bam files")
parser.add_argument("-r", "--ref", dest="ref", action="store", required=True, help="Name of reference organism")
parser.add_argument("-d", "--date", dest="date", action="store", help="date of experiment")
parser.add_argument("-q", "--mapQ", dest="mapQ", action="store",help="Minimum mapping quality for reads")
parser.add_argument("-m", dest="minCount", default=5.0, type=float, help="min average counts for given transcript  [%(default)s]")
parser.add_argument("-v", dest="verbose_level", default=2, action="store", help="Indicate verbose level for transcript2count.py (None, 1, 2)")

args = parser.parse_args()


repeatIt(args.bamfolder, args.exp, args.bedfile, args.ref, args.date, args.mapQ, args.minCount, args.verbose_level)







##############################
# For c. elegans Total from 2016-09-23
#name = "Total-2016-09-23a"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/Total_1_S1_16022_ACAGTG" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "Total-2016-09-23b"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/Total_1_V1_16021_ATCACG" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "Total-2016-09-23c"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/Total_2_S1_16024_GTGGCC" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "Total-2016-09-23d"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/Total_2_V1_16023_TAGCTT" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
##############################


##############################
# For c. elegans PolyA from 2016-09-23
#name = "PolyA-2016-09-23a"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/S1_PolyA_1_16017_ACTTGA" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "PolyA-2016-09-23b"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/S1_PolyA_2_16019_CGTACG" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "PolyA-2016-09-23c"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/V1_PolyA_1_16016_GCCAAT" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "PolyA-2016-09-23d"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/V1_PolyA_2_16018_GGCTAC" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
##############################


##############################
# For c. elegans Enriched from 2016-09-23
#name = "ENRICH-2016-09-23a"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/ENRICH_PARS_C_elegans_16020_ACAGTG_S1" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "ENRICH-2016-09-23b"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/ENRICH_PARS_C_elegans_16020_ATCACG_V1" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "ENRICH-2016-09-23c"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/ENRICH_PARS_C_elegans_16020_GTGGCC_S1" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
#
#name = "ENRICH-2016-09-23d"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/C_ELEGANS/WS256/c_elegans.PRJNA13758.WS256.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-09-23/ENRICH_PARS_C_elegans_16020_TAGCTT_V1" %home
#ref = "c_elegans"
#repeatIt(bamfolder, name, bed, ref)
##############################










#name = "2016-03-18a"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HT29.bed' %home
##bed = 'test.bed'
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-03-18/aa_1_V1_REP1_12905_AGTCAA" %home
#ref = "human"
#repeatIt(bamfolder, name, bed)
#
#name = "2016-03-18b"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HT29.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-03-18/aa_2_S1_REP1_12906_ACAGTG" %home
#ref = "human"
#repeatIt(bamfolder, name, bed)
#
#name = "2016-03-18c"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HT29.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-03-18/aa_3_V1_REP2_12907_TAGCTT" %home
#ref = "human"
#repeatIt(bamfolder, name, bed)

#name = "2016-03-18d"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HT29.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-03-18/aa_4_S1_REP2_12908_GTGGCC" %home
#ref = "human"
#repeatIt(bamfolder, name, bed)
###############################


#
###############################
#name = "2016-05-18a"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HCT116.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-05-18/HCT_116_S1_REP1_13369_TGACCA_MERGED" %home #***
#ref = "human"
#repeatIt(bamfolder, name, bed)
#
#name = "2016-05-18b"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HCT116.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-05-18/HCT_116_V1_REP1_13368_CGATGT_MERGED" %home #***
#ref = "human"
#repeatIt(bamfolder, name, bed)
##*** a and b must use these "MERGED" bams because the data comes from two runs that were merged after rep1 was repeated
#
#name = "2016-04-19c"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HCT116.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-04-19/4_S1_REP2_13371_CCGTCC" %home #name is different bc of library prep
#ref = "human"
#repeatIt(bamfolder, name, bed)
#
#name = "2016-04-19d"
#bed = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/HUMAN/HCT116.bed' %home
#bamfolder = "%s/jwillis/lncRNA_PROJECT/2016-04-19/HCT_116_V1_REP2_13370_CAGATC" %home
#ref = "human"
#repeatIt(bamfolder, name, bed)
###############################























###############################
##~ bed = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/CANDIDA/refs/cglabrata_controls.bed"
#bed = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/CANDIDA/refs/SNCA.bed"
#name = "2015-01-09a"
##~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09/aa_C_gla_V1_rep1_8962_ACAGTG" 
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09_SNCA/aa_C_gla_V1_rep1_8962_ACAGTG.fa" 
#repeatIt(bamfolder, name, bed)
#
#name = "2015-01-09b"
##~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09/aa_C_gla_S1_rep1_8963_GCCAAT" 
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09_SNCA/aa_C_gla_S1_rep1_8963_GCCAAT.fa" 
#repeatIt(bamfolder, name, bed)
##~ 
#name = "2015-01-09c"
##~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09/aa_C_gla_V1_rep2_8964_CTTGTA" 
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09_SNCA/aa_C_gla_V1_rep2_8964_CTTGTA.fa" 
#repeatIt(bamfolder, name, bed)
##~ 
#name = "2015-01-09d"
##~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09/aa_C_gla_S1_rep2_8965_GTGAAA" 
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/2015-01-09_SNCA/aa_C_gla_S1_rep2_8965_GTGAAA.fa" 
#repeatIt(bamfolder, name, bed)
###############################












#~ name = "2014-04-29a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2014-04-29/TotalYeast_S1_rep1_6740_GGCTAC"
#~ repeatIt(bamfolder, name, bed)









###############################
#~ bed = "/users/tg/dloska/scripts/yeastR62_controls.bed"
#~ bed = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/YEAST/scerevisiae_R61/yeastR61_controls.bed"

#~ name = "2014-04-29a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2014-04-29/TotalYeast_S1_rep1_6740_GGCTAC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2014-04-29b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2014-04-29/TotalYeast_S1_rep2_6742_GAGTGG"
#~ repeatIt(bamfolder, name, bed)

#~ name = "2014-04-29c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2014-04-29/TotalYeast_V1_rep1_6739_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2014-04-29d"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2014-04-29/TotalYeast_V1_rep2_6741_CGTACG"
#~ repeatIt(bamfolder, name, bed)
###############################











#~ name = "test"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig41_V1_5679/unmappedResults/Dig41_V1_5679"
#~ repeatIt(bamfolder, name, bed)




########################################bed
########################################bed
########################################bed
#~ bed = "/home/dloska/Desktop/RNA/tab/tabs/reference-Hall.bed"
#~ 
#~ 
#~ name = "2013-11-18a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1A_V1_5434_ATCACG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1A_S1_5435_CGATGT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2A_V1_5436_TTAGGC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2A_S1_5437_TGACCA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3A_V1_5438_ACAGTG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3A_S1_5439_GCCAAT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18d"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4A_V1_5440_CAGATC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4A_S1_5441_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18e"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5A_V1_5442_GATCAG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5A_S1_5443_TAGCTT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18f"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1B_V1_5444_GGCTAC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1B_S1_5445_CTTGTA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2B_V1_5446_AGTCAA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2B_S1_5447_AGTTCC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18h"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3B_V1_5448_ATGTCA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3B_S1_5449_CCGTCC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18i"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4B_V1_5450_GTAGAG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4B_S1_5451_GTCCGC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18j"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5B_V1_5452_GTGAAA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5B_S1_5453_GTGGCC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ 
#~ 
#~ 
#~ 
#~ 
#~ 
########################################bed
########################################bed
########################################bed


#~ name = "2014-01-22a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/1_V1_0_03_A_6010_GTTTCG"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2014-01-22b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/2_V1_0_03_B_6011_CGTACG"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2014-01-22c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/3_V1_0_03_C_6012_GAGTGG"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/5_T1_1_A_6013_GGTAGC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/7_T1_1_B_6014_ATCACG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/8_T1_1_C_6015_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22g"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/13_A_0_05_A_6016_TTAGGC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/15_A_0_05_B_6017_TGACCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22i"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/16_A_0_05_C_6018_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22j"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/18_MIX_A_6019_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22k"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/19_MIX_B_6020_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22l"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/ester-2014-01-22/20_MIX_C_6021_ACTTGA"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ 
					#~ 
					#~ name = "2013-11-25a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig40_NoDigest_5675/unmappedResults/Dig40_NoDigest_5675"
					#~ repeatIt(bamfolder, name, bed)
#~ name = "2013-11-26a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig40_V1_5676/unmappedResults/Dig40_V1_5676"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2013-11-29a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig40_S1_5677/unmappedResults/Dig40_S1_5677"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "2013-12-10a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig41_NoDigest_5678/unmappedResults/Dig41_NoDigest_5678"
			
					#~ repeatIt(bamfolder, name, bed)
#~ name = "2013-12-03a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig41_V1_5679/unmappedResults/Dig41_V1_5679"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2013-12-10b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig41_S1_5680/unmappedResults/Dig41_S1_5680"  #bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig41_S1_5680/unmappedResults/Dig41_S1_5680/"
#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ 
					#~ name = "0000-00-02a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_Nodig_4768_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_V1_0_004_4769_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_V1_0_01_4770_CTTGTA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_V1_0_02_4771_ATGTCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_S1_20_4772_CCGTCC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_S1_50_4773_GTAGAG"
					#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-02g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_S1_200_4774_GTGAAA"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig39_S1_500_4775_GTGGCC"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ 
					#~ name = "0000-00-01a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_Nodig_4712_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_V1_0_004_4713_TGACCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_V1_0_01_4714_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_V1_0_02_4715_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_S1_20_4716_AGTCAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_S1_50_4717_AGTTCC"
					#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-01g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_S1_200_4718_GTCCGC"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig38_S1_500_4719_GAGTGG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-00l"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig_37_S1_500_4647_CAAAAG"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig_37_V1_0_02_4646_CGTACG"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00k"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig_36_S1_20_4645_GTGAAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig_36_V1_0_004_4644_GGCTAC"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
#~ name = "0000-00-00j"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig34_NO_dig_4642_TTAGGC"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig_34_S1_200_4643_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00i"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig33_Nodig_4573/"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig33_S1_4575/"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig33_V1_4574/"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig31_V1_4479_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig31_S1_4480_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00g"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig29_V1_4050_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig29_S1_4051_GTGAAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-29-39/Dig29_Nodig_4049_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig29_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig29_V1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig29_S1"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
#~ name = "0000-00-00e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig26_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig26_V1"
#~ repeatIt(bamfolder, name, bed)
					
#~ name = "0000-00-00d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig25_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig25_V1"
#~ repeatIt(bamfolder, name, bed)
					
#~ name = "0000-00-00c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig20_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig20_V1"
#~ repeatIt(bamfolder, name, bed)

					#~ name = "0000-00-00b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17B_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17B_S1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17B_V1"
					#~ repeatIt(bamfolder, name, bed)

					#~ name = "0000-00-00a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17_S1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/Dig-17-29/Dig17_V1"
					#~ repeatIt(bamfolder, name, bed)
					
					
					
					
					
					
					
########################################################################	
########################################################################	
########################################################################	
########################################################################	
############## ref R61, bed R61 ########################################	
########################################################################	
########################################################################	
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################


#~ bed = "/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/GENOMES/YEAST/scerevisiae_R61/yeastR61_controls.bed"
	#~ 
					#~ 
#~ name = "2014-04-29a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/TotalYeast_S1_rep1_6740_GGCTAC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2014-04-29b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/TotalYeast_S1_rep2_6742_GAGTGG"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2014-04-29c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/TotalYeast_V1_rep1_6739_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2014-04-29d"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/TotalYeast_V1_rep2_6741_CGTACG"
#~ repeatIt(bamfolder, name, bed)



#~ name = "2013-11-18a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1A_V1_5434_ATCACG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1A_S1_5435_CGATGT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2A_V1_5436_TTAGGC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2A_S1_5437_TGACCA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3A_V1_5438_ACAGTG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3A_S1_5439_GCCAAT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18d"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4A_V1_5440_CAGATC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4A_S1_5441_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18e"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5A_V1_5442_GATCAG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5A_S1_5443_TAGCTT"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18f"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1B_V1_5444_GGCTAC"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H1B_S1_5445_CTTGTA"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2B_V1_5446_AGTCAA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H2B_S1_5447_AGTTCC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18h"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3B_V1_5448_ATGTCA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H3B_S1_5449_CCGTCC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18i"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4B_V1_5450_GTAGAG"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H4B_S1_5451_GTCCGC"
#~ repeatIt(bamfolder, name, bed)
#~ 
#~ name = "2013-11-18j"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5B_V1_5452_GTGAAA"
#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R62/2013-11-18/H5B_S1_5453_GTGGCC"
#~ repeatIt(bamfolder, name, bed)
#~ 






#~ name = "2014-01-22a"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_1_V1_0_03_A_6010_GTTTCG"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2014-01-22b"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_2_V1_0_03_B_6011_CGTACG"
#~ repeatIt(bamfolder, name, bed)
#~ name = "2014-01-22c"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_3_V1_0_03_C_6012_GAGTGG"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_5_T1_1_A_6013_GGTAGC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_7_T1_1_B_6014_ATCACG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_8_T1_1_C_6015_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22g"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_13_A_0_05_A_6016_TTAGGC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_15_A_0_05_B_6017_TGACCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22i"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_16_A_0_05_C_6018_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22j"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_18_MIX_A_6019_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22k"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_19_MIX_B_6020_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "2014-01-22l"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/aa_20_MIX_C_6021_ACTTGA"
					#~ repeatIt(bamfolder, name, bed)
					
					
					
					#name = "2013-11-25a"
					#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig40_NoDigest_5675/unmappedResults/Dig40_NoDigest_5675"
					#repeatIt(bamfolder, name, bed)
#name = "2013-11-26a"
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig40_V1_5676/unmappedResults/Dig40_V1_5676"
#repeatIt(bamfolder, name, bed)
#name = "2013-11-29a"
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig40_S1_5677/unmappedResults/Dig40_S1_5677"
#repeatIt(bamfolder, name, bed)
					#name = "2013-12-10a"
					#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig41_NoDigest_5678/unmappedResults/Dig41_NoDigest_5678"
					#repeatIt(bamfolder, name, bed)
#name = "2013-12-03a"
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig41_V1_5679/unmappedResults/Dig41_V1_5679"
#repeatIt(bamfolder, name, bed)
#name = "2013-12-10b"
#bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/R62/ester-40-41-sh/Dig41_S1_5680/unmappedResults/Dig41_S1_5680"  #bamfolder = "/users/tg/dloska/lncRNA_PROJECT/shfiles/ester-40-41-sh/Dig41_S1_5680/unmappedResults/Dig41_S1_5680/"
#~ repeatIt(bamfolder, name, bed)
					
					
					#~ name = "0000-00-02a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_Nodig_4768_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_V1_0_004_4769_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_V1_0_01_4770_CTTGTA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_V1_0_02_4771_ATGTCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_S1_20_4772_CCGTCC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_S1_50_4773_GTAGAG"
					#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-02g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_S1_200_4774_GTGAAA"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-02h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig39_S1_500_4775_GTGGCC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_Nodig_4712_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_V1_0_004_4713_TGACCA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_V1_0_01_4714_ACAGTG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_V1_0_02_4715_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_S1_20_4716_AGTCAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_S1_50_4717_AGTTCC"
					#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-01g"
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_S1_200_4718_GTCCGC"
#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-01h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig38_S1_500_4719_GAGTGG"
					#~ repeatIt(bamfolder, name, bed)
					#~ name = "0000-00-00l"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig_37_S1_500_4647_CAAAAG"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig_37_V1_0_02_4646_CGTACG"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00k"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig_36_S1_20_4645_GTGAAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig_36_V1_0_004_4644_GGCTAC"
					#~ repeatIt(bamfolder, name, bed)
					
#~ name = "0000-00-00j"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig34_NO_dig_4642_TTAGGC"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig_34_S1_200_4643_ACTTGA"
#~ repeatIt(bamfolder, name, bed)
					
					#~ name = "0000-00-00i"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig33_Nodig_4573/"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig33_S1_4575/"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig33_V1_4574/"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00h"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig31_V1_4479_GCCAAT"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig31_S1_4480_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00g"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_V1_4050_CAGATC"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_S1_4051_GTGAAA"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_Nodig_4049_CGATGT"
					#~ repeatIt(bamfolder, name, bed)
					#~ 
					#~ name = "0000-00-00f"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_V1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig29_S1"
					#~ repeatIt(bamfolder, name, bed)
					
#~ name = "0000-00-00e"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig26_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig26_V1"
#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-00d"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig25_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig25_V1"
#~ repeatIt(bamfolder, name, bed)
#~ name = "0000-00-00c"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig20_S1"
					#~ repeatIt(bamfolder, name, bed)
#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig20_V1"
#~ repeatIt(bamfolder, name, bed)
#~ 
					#~ name = "0000-00-00b"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17B_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17B_S1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17B_V1"
					#~ repeatIt(bamfolder, name, bed)
#~ 
					#~ name = "0000-00-00a"
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17_control"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17_S1"
					#~ repeatIt(bamfolder, name, bed)
					#~ bamfolder = "/users/tg/dloska/lncRNA_PROJECT/results/results_R61/Dig17_V1"
					#~ repeatIt(bamfolder, name, bed)
					
					
					

