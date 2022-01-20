#!/usr/bin/env python
desc="""Count number of 5'-end reads for all transcripts at given load.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Dublin/Barcelona, 4/07/2012

versions:
- 0.3 - unstranded option added (to be done!)
- 0.2 - solved issue of including reverse reads for + transcripts (thanks to Jose) 
- 0.1 
"""

import argparse, os, sys
import pysam
import numpy as np
from datetime import datetime
from genome_annotation import load_transcripts_bed
########################################################################
def get_reads_5ends( samfile,ref,start,end,mapq,reverse ):
    """Return number of reads starting at each position
    of defined region.
    """
    #define empty counts list
    counts = [ 0 for x in range(end-start) ]
    #process all reads in region
    try:
        for sam in samfile.fetch( ref,start,end+1 ):
            #skip if low quality
            if sam.mapq < mapq:
                continue
            ##skip if read comes from opposite strand
            #if the transcript is reverse and read is not reverse
            #if the transcript is forward and the read is reverse
            if sam.is_reverse and not reverse or not sam.is_reverse and reverse: 
                continue
            #get 5' position of read and 1bp upstream
            if reverse:
                pos = sam.aend - 1 
            else:
                pos = sam.pos
            #check if start in transcript range
            tpos = pos-start
            if 0 <= tpos < len(counts):
                #update counts
                counts[tpos] += 1
        return counts, 'good', None
    except ValueError as e:
        #For when pysam.calignmentfile throws ValueError: start out of range (-1)
#        for sam in samfile.fetch( ref,start+1,end+1 ):
#            if sam.mapping_quality < mapq:
#                continue
#            if sam.is_reverse and not reverse or not sam.is_reverse and reverse: 
#                continue
#            if reverse:
#                pos = sam.aend - 1 
#            else:
#                pos = sam.pos
#            tpos = pos-start
#            if 0 <= tpos < len(counts):
#                counts[tpos] += 1
        return counts, 'bad', e
    except IOError as e:
        #For when pysam.calignmentfile throws IOError: truncated file
        return counts, 'truncated', e
########################################################################
def process_exons():
    
    
    return counts
########################################################################
def process_transcripts( bed,bam,mapq,minCount,offset,verbose_level ):
    """
    """
    #first load transcripts
    if verbose_level:
        sys.stderr.write("Loading BED file...\n")
    tbed = datetime.now()
    transcripts = load_transcripts_bed( bed, PARS=True )#oneoff=True ; print transcripts
    if verbose_level:
        sys.stderr.write("## Time to load %s unique transcripts: %s\n" %(len(transcripts), datetime.now() - tbed ))
        #then parse bam file for every transcript
        sys.stderr.write("Parsing BAM file for %s BED entries...\n\n" % len(transcripts) )
    
    ### Process transcripts
    samfile = pysam.Samfile( bam,"rb" )
    refs    = set( samfile.references )
#    sys.stderr.write("%s\t%s\n" %(refs, len(refs)))
    bad_start_ranges = open( '%s/start_range_min1_transcripts.txt' %bam.rsplit('/',1)[0], 'w' )
    truncated_transcripts = open( '%s/truncated_transcripts.txt' %bam.rsplit('/',1)[0], 'w' )
    low_exp_transcripts = open( '%s/low_exp_transcripts.txt' %bam.rsplit('/',1)[0], 'w' )

    i = 0
    k = 0
    le = 0 #for lowly expressed transcripts
    ge = 0 #for transcripts with average expression above minCount
    le_list, ge_list = [], []
    t_process_begin = datetime.now()
    percents = [10,20,30,40,50,60,70,80,90,100]
    for transcript in transcripts:
        #if transcript not in ("YAL001C","YAR002W","YAR002C-A","YAR003W","YAL003W"): continue
        i += 1
        
        #TODO: make list of already printed trqnscripts to avoid searching
#        if transcript in already_printed:
        #TODO: will need to change output to an argparse option in order to get the output directory to see whats been printed already 
        
        
        if int(verbose_level) > 1:
            sys.stderr.write(" %s %s    \r" % (i,transcript) )
        elif verbose_level:
            if round((float(i)/len(transcripts))*100, 0) in percents:
                sys.stderr.write("%3s%% of transcripts processed\tTime: %s...\n" % (percents.pop(0), datetime.now() - t_process_begin) )
        #ref
        ref    = transcripts[transcript]["chromosome"]
        strand = transcripts[transcript]["strand"]
        #check if ref indeed in bam
        if ref not in refs:
            if verbose_level:
                sys.stderr.write( " Warning: %s not in BAM file\n" % ( ref, ) )
            continue
        #define strand
        reverse = False
        if strand == "-":
            reverse = True            
        #process all exons (and UTRs)
        counts = []
        exonCount = len(transcripts[transcript]["intervals"])
        for ii in range(exonCount):
            #get start and end
            start,end,score = transcripts[transcript]["intervals"][ii]
            #add 1bp for transcript start if reverse
            if ii   == 0 and reverse:
                start -= offset
            #add 1bp for transcript end if not reverse
            if ii+1 == exonCount and not reverse:
                end   += offset
            #get counts
            r5e, start_range, error_message = get_reads_5ends( samfile,ref,start,end,mapq,reverse )
            if start_range == 'bad':
                #write a list of the transcripts that threw ValueError: start out of range (-1)
                # these have been adjusted for the function get_reads_5ends() above
                bad_start_ranges.write('%s\t%s\texon %s\t%s\t%s | %s\n' %(transcript, ref, ii, start, end, error_message))
                continue
            elif start_range == 'truncated':
                #write a list of the transcripts that threw IOError: truncated file
                truncated_transcripts.write('%s\t%s\texon %s\t%s\t%s | %s\n' %(transcript, ref, ii, start, end, error_message))
                continue
#            sys.stderr.write("%s exon #%s good\n" %(transcript, ii))
            counts += r5e
                
        #reverse (or not) and add tailing base count
        if reverse:
            counts.reverse()
        #move all by 1bp - as 1bp added to transcript start)
        counts = counts[offset:]
        #output line
        if len(counts) > 0 and np.mean(counts) > minCount:
            #only if actual data is present for the given transcript
            ge += 1
            ge_list.append(np.mean(counts))
            line   = "%s\t%s\t%s\n" % ( transcript,len(counts),";".join(str(x) for x in counts) ) 
            sys.stdout.write( line )
        else:
            le += 1
            if len(counts) > 0:
                le_list.append(np.mean(counts))
                low_exp_transcripts.write('%s\t%s\n' %(transcript, np.mean(counts)))
            else:
                le_list.append(0)
                low_exp_transcripts.write('%s\t%s\n' %(transcript, '0.0'))
        
        k += 1

    if int(verbose_level) == 1:
        sys.stderr.write("%s transcripts found in BAM file and processed.\n\n" %k)
        sys.stderr.write("%s of those transcripts had an average count per site above the minimum (%s).\n" %(ge, minCount))
        sys.stderr.write("The mean counts per site for these was %s, with a max of %s and a min of %s\n\n" %(np.mean(ge_list),max(ge_list),min(ge_list)))
        sys.stderr.write("%s transcripts had an average count per site below the minimum (%s).\n" %(le, minCount))
        sys.stderr.write("The mean counts per site for these was %s, with a max of %s and a min of %s\n\n" %(np.mean(le_list),max(le_list),min(le_list)))
    
    bad_start_ranges.close()
    truncated_transcripts.close()
    low_exp_transcripts.close()
########################################################################
def main():

    usage  = "%(prog)s [options]"
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )
    
    parser.add_argument("-v", dest="verbose_level", default=None, action="store", help="Indicate verbose level (None, 1, 2)")
    parser.add_argument('--version', action='version', version='0.2')
    parser.add_argument("-a", dest="bam",  type=file,
                        help="bam file (sorted)     [%(default)s]")
    parser.add_argument("-b", dest="bed",  type=file,
                        help="bed file              [%(default)s]")
    parser.add_argument("-q", dest="mapq", default=0, type=int,
                        help="min mapping quality   [%(default)s]")
    parser.add_argument("-m", dest="minCount", default=5.0, type=float,
                        help="min average counts for given transcript  [%(default)s]")
    parser.add_argument('--offset', dest='offset', default=1, type=int,
                        help="reads are counted for base upstream [%(default)s]")
  
    o = parser.parse_args()
    if o.verbose_level:
        sys.stderr.write( "Options: %s\n\n" % str(o) )

    #process all transcripts
    process_transcripts( o.bed.name,o.bam.name,o.mapq,o.minCount,o.offset,o.verbose_level )
########################################################################
########################################################################
if __name__=='__main__': 
    t0 = datetime.now()
    main()
    dt = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
  
