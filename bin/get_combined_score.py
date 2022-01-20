#! usr/bin/python
# (c) Ton Gabaldon 2014

#Example of usage:
#$ python get_combined_score_v1.1.py -i TETp4p6 <options>

import argparse, os, sys#, glob
import numpy as np
import tabscorelib as tsl
from termcolor import colored


# global variables
home = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
data = '%s/data' %home
src = '%s/bin' %home

#####################################################################################
#####################################################################################

def bp_nucs(calls, mol_seq):
    '''Count number of times each nucleotide is in a paired position and return list of them with their position in sequence'''
    mol_seq = ['U' if n=='T' else n for n in mol_seq]
    nucs = []
    for i in xrange(len(mol_seq)):
        if calls[i] == 1:
            nucs.append((mol_seq[i], i+1))
#    print 'Nucs in pairs:', nucs, len(nucs)
    g = 0
    c = 0
    a = 0
    u = 0
    for n in nucs:
        if n[0] == 'G': g += 1
        elif n[0] == 'C':   c += 1
        elif n[0] == 'A':   a += 1
        elif n[0] == 'U':   u += 1
    print 'G = %i' %g
    print 'C = %i' %c
    print 'A = %i' %a
    print 'U = %i' %u
    print 'Total = %i' %(g+c+a+u)
    print '\n'
    return nucs
#####################################################################################
def pairs_from_ref_ct(molname):
    '''Create tuples of info about pairs from reference structure'''
    if molname == 'TETp4p6':
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.ct' %(data,molname), 'r')
#        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s-PDB.ct' %(data,molname), 'r')
    elif molname == 'SRA' or molname == 'B2' or molname == 'U1':
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s-ref.ct' %(data,molname), 'r')
    else:
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.ct' %(data,molname), 'r')
    
    pairs_ct = []
    nucs = set(['A', 'C', 'G', 'U', 'T', 'X', 'a', 'c', 'g', 'u', 't', 'x'])  #To use only lines with the position and pairing info
    allnuc = []
    
    for line in ct_file:
        if not line.startswith('#'):     #Ignore commented lines at beginning
            line = line.split()
            if len(line) > 1 and line[1] in nucs:    #Ignore empty lines and first line with number of positions and additional info
                allnuc.append(line[1].upper())
                if line[4] != '0':                #Pairing is found in 5th column of ct files
                    pairs_ct.append((line[1].upper(), int(line[5]), int(line[4])))  #(NUC, pos, pos of pair)
                    #TODO: add an 'and' statement to ensure pairing nuc is within range of sequence length (no pseudoknots)
    ct_file.close()
    
#    print pairs_ct, len(pairs_ct)
#    print allnuc, len(allnuc)
    pairs_ct = [('U',n[1],n[2],allnuc[n[2]-1]) if n[0]=='T' else (n[0],n[1],n[2],allnuc[n[2]-1]) for n in pairs_ct] #(Nuc, pos, pair pos, pair nuc)
    #pairs_ct[0] = Nucleotide
    #pairs_ct[1] = Position of nucleotide in reference sequence
    #pairs_ct[2] = Position of pairing nucleotide
    #pairs_ct[3] = Paired nucleotide
#    print pairs_ct, len(pairs_ct), '\n'
    
    gc = ['GC', 0]
    au = ['AU', 0]
    gu = ['GU', 0]
    ac = ['AC', 0]
    ga = ['GA', 0]
    cu = ['CU', 0]
    aa = ['AA', 0]
    gg = ['GG', 0]
    uu = ['UU', 0]
    cc = ['CC', 0]
    
    for p in pairs_ct:
        if (p[0]=='G' and p[3]=='C'):# or (p[0]=='C' and p[3]=='G'):
            gc[1] += 1
        elif (p[0]=='A' and p[3]=='U'):# or (p[0]=='U' and p[3]=='A'):
            au[1] += 1
        elif (p[0]=='G' and p[3]=='U'):# or (p[0]=='U' and p[3]=='G'):
            gu[1] += 1
        elif (p[0]=='A' and p[3]=='C'):# or (p[0]=='C' and p[3]=='A'):
            ac[1] += 1
        elif (p[0]=='G' and p[3]=='A'):# or (p[0]=='A' and p[3]=='G'):
            ga[1] += 1
        elif (p[0]=='C' and p[3]=='U'):# or (p[0]=='U' and p[3]=='C'):
            cu[1] += 1
        elif (p[0]=='A' and p[3]=='A'):
            aa[1] += 1
        elif (p[0]=='G' and p[3]=='G'):
            gg[1] += 1
        elif (p[0]=='U' and p[3]=='U'):
            uu[1] += 1
        elif (p[0]=='C' and p[3]=='C'):
            cc[1] += 1
#        else:
#            print p
    aa[1] /= 2
    gg[1] /= 2
    uu[1] /= 2
    cc[1] /= 2
    
    plist = sorted([gc, au, gu, ac, ga, cu, aa, gg, uu, cc], key=lambda x: x[1], reverse=True)
    print colored('Pairs in reference:', 'yellow')
    i = 0
    for pair in plist:
        if pair[1] != 0:
            print pair[0], '=', pair[1], 'pairs'
            i += 1
    if i == 0:
        print colored('No pairs', 'red')
    print '\n'
    return pairs_ct
#####################################################################################
def correct_bps(paired_nucs, pairs_ct):
    '''Check how many full pairs are predicted correctly from combined scores'''
    index = 0
    ind2 = 0
    corresp = []
    pred_pairs = []
    for pn in paired_nucs:
        for pc in pairs_ct:
            if (pn[0], pn[1]) == (pc[0], pc[1]):
                
                for n_p in paired_nucs:
                    if pc[2] == n_p[1]:
                        corresp.append(n_p)
                        ind2 += 1
                index += 1
                if len(corresp) > 0 and pc[2] == corresp[ind2-1][1]:
#                    print index, pc, corresp[ind2-1], ind2
                    pred_pairs.append(pc)
#                else:
#                    print index, pc, 'No match', ind2
                  
    gc = ['GC', 0]
    au = ['AU', 0]
    gu = ['GU', 0]
    ac = ['AC', 0]
    ga = ['GA', 0]
    cu = ['CU', 0]
    aa = ['AA', 0]
    gg = ['GG', 0]
    uu = ['UU', 0]
    cc = ['CC', 0]
    
    for p in pred_pairs:
        if (p[0]=='G' and p[3]=='C'):# or (p[0]=='C' and p[3]=='G'):
            gc[1] += 1
        elif (p[0]=='A' and p[3]=='U'):# or (p[0]=='U' and p[3]=='A'):
            au[1] += 1
        elif (p[0]=='G' and p[3]=='U'):# or (p[0]=='U' and p[3]=='G'):
            gu[1] += 1
        elif (p[0]=='A' and p[3]=='C'):# or (p[0]=='C' and p[3]=='A'):
            ac[1] += 1
        elif (p[0]=='G' and p[3]=='A'):# or (p[0]=='A' and p[3]=='G'):
            ga[1] += 1
        elif (p[0]=='C' and p[3]=='U'):# or (p[0]=='U' and p[3]=='C'):
            cu[1] += 1
        elif (p[0]=='A' and p[3]=='A'):
            aa[1] += 1
        elif (p[0]=='G' and p[3]=='G'):
            gg[1] += 1
        elif (p[0]=='U' and p[3]=='U'):
            uu[1] += 1
        elif (p[0]=='C' and p[3]=='C'):
            cc[1] += 1
                
    aa[1] /= 2
    gg[1] /= 2
    uu[1] /= 2
    cc[1] /= 2
    
    plist = sorted([gc, au, gu, ac, ga, cu, aa, gg, uu, cc], key=lambda x: x[1], reverse=True)
    print colored('Pairs correct in predicted:', 'yellow')
    i = 0
    for pair in plist:
        if pair[1] != 0:
            print pair[0], '=', pair[1], 'pairs'
            i += 1
    if i == 0:
        print colored('No pairs', 'red')
    print '\n'                

#####################################################################################
#####################################################################################

def removepos(combined, correct_x, incorrect_x, corfile, incfile):
    cor_pos = []
    inc_pos = []
    for i in xrange(len(combined)):
        if i in correct_x:
            cor_pos.append(i)   #Create list of indices of correct positions
        elif i in incorrect_x:
            inc_pos.append(i)   #Create list of indices of incorrect positions
    
#    print cor_pos, len(cor_pos)
#    print inc_pos, len(inc_pos)
    
    corfile.write(str(cor_pos)+'\n')
    incfile.write(str(inc_pos)+'\n')
#####################################################################################
#####################################################################################

def get_tabfilelist(molname):
    '''Write list of tab file names to file to be read by function tsl.fileList() and 
       return molname_full and real_path (path used for opening file with function tsl.readtabfile()'''
    
    PATHS = ['/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/refR61_bedR61',
             '/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/candida',
             '/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/human',
             '/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/c_elegans',
             '/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/refR62_bedR62',
             '/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/refR62_bedR61']

    if molname.startswith('RDN'):
        final_concs = ['2014-04-29a','2014-04-29b','2014-04-29c','2014-04-29d']
        molname_full = molname
    elif molname == 'smHotair_2015-01-09':
        final_concs = [ '2015-01-09a','2015-01-09b','2015-01-09c','2015-01-09d' ]
        molname_full = molname
        molname = 'smHotair' #because will need just the base name in many places
    elif molname == 'smHotair_2016-03-18':
        final_concs = [ '2016-03-18a','2016-03-18b','2016-03-18c','2016-03-18d' ]
        molname_full = molname
        molname = 'smHotair' #because will need just the base name in many places
    elif molname == 'smHotair':
        sm_date = raw_input("Which smHotair do you want, (a) 2015-01-09 or (b) 2016-03-18? ")
        if sm_date in ['a','A', '(a)', '2015', '2015-01-09']:
            final_concs = [ '2015-01-09a','2015-01-09b','2015-01-09c','2015-01-09d' ]
            molname_full = 'smHotair_2015-01-09'
        elif sm_date in ['b','B', '(b)', '2016', '2016-03-18']:
            final_concs = [ '2016-03-18a','2016-03-18b','2016-03-18c','2016-03-18d' ]
            molname_full = 'smHotair_2016-03-18'
    else:
        final_concs = ['0000-00-00c','0000-00-00d','0000-00-00e','0000-00-00j','0000-00-01g','0000-00-02g',
                   '2013-11-18a','2013-11-18b','2013-11-18c','2013-11-18d','2013-11-18e',
                   '2013-11-18f','2013-11-18g','2013-11-18h','2013-11-18i','2013-11-18j',
                   '2014-01-22a','2014-01-22b','2014-01-22c',
                   '2014-04-29a','2014-04-29b','2014-04-29c','2014-04-29d',
                   '2015-01-09a','2015-01-09b','2015-01-09c','2015-01-09d',
                   '2016-03-18a','2016-03-18b','2016-03-18c','2016-03-18d',
                   '2016-05-18a','2016-05-18b','2016-05-18c','2016-05-18d',
                   'ENRICH-2016-09-23a','ENRICH-2016-09-23b','ENRICH-2016-09-23c','ENRICH-2016-09-23d',
                   'PolyA-2016-09-23a','PolyA-2016-09-23b','PolyA-2016-09-23c','PolyA-2016-09-23d',
                   'Total-2016-09-23a','Total-2016-09-23b','Total-2016-09-23c','Total-2016-09-23d']
        molname_full = molname
        
    tabs = []
    append_tabs = tabs.append
    for PATH in PATHS:
        if len(tabs) == 0:  #Check to see if tabs were already found in one of the other directories
            real_path = PATH  #path to be used to open file with readtabfile() function
            for exp in sorted(os.listdir(PATH)):
                if exp in final_concs:
#                print PATH, exp
                    if os.path.isdir(PATH+'/'+exp):
                        if os.path.exists('%s/%s/%s_%s_V1.tab' %(PATH,exp,molname,exp)):
#                            print '%s/%s/%s_%s_V1.tab' %(PATH,exp,molname,exp)
                            append_tabs('%s_%s_V1.tab' %(molname, exp))
                        if os.path.exists('%s/%s/%s_%s_S1.tab' %(PATH,exp,molname,exp)):
#                            print '%s/%s/%s_%s_S1.tab' %(PATH,exp,molname,exp)
                            append_tabs('%s_%s_S1.tab' %(molname, exp))
        else:
            break

    tabfilelist = open(molname+'.txt', 'w')
    tabs = sorted(tabs)
    for t in tabs:
        tabfilelist.write(t+'\n')
    tabfilelist.close()
    
    return molname_full, real_path, None, None #None's are for tabfile names..wont be used in this case
#####################################################################################
def get_tabfilelist_nextPARS(molname, exp_dir=None):
    '''Write list of tab file names to file to be read by function tsl.fileList() and 
       return molname_full and real_path (path used for opening file with function tsl.readtabfile()'''
    
    molname_full = molname
    
    if molname in ['TETp4p6', 'RDN5-1', 'RDN18-1', 'RDN25-1', 'RDN58-1']:
        tabs = []
        for exp in os.listdir('%s/%s' %(data, molname)):
            tabs.append(exp)

#        tabOutputFolder = glob.glob('%s/%s/%s*.tab' %(data, molname, molname))
#        if len(tabOutputFolder) > 0:
#            for exp in tabOutputFolder:
#                tabs.append(exp.split('/')[-1])
        
        tabfilelist = open(molname+'.txt', 'w')
        tabs = sorted(tabs)
        for t in tabs:
            tabfilelist.write(t+'\n')
        tabfilelist.close()
        
        return molname_full, '%s/%s' %(data, molname), None, None #None's are for tabfile names..wont be used in this case
    
    else:
        V1_tabs, S1_tabs = [], []
        v_exps, s_exps = [], []
        
        # in this case, exp_dir argument should be included and must be a directory containing tab files, so first must check:
        if os.path.isdir(exp_dir):
            for exp in sorted(os.listdir(exp_dir)):
                if exp.endswith('.tab'):
                    with open('%s/%s' %(exp_dir, exp)) as opentab:
                        for l in opentab:
                            if molname in l:
                                if "V1" in exp:
                                    V1_tabs.append(l)
                                    v_exps.append(exp)
                                elif "S1" in exp:
                                    S1_tabs.append(l)
                                    s_exps.append(exp)
                                else:
                                    print "The tab files should indicate which enzyme was used (V1/S1)"
                                    print "Please rename them accordingly"
                                    exit()
                                break
            
            return V1_tabs, S1_tabs, molname, v_exps, s_exps, None # return None in place of hte variable real_path, which wont be used in this case
                
        else:
            print "If not using one of the provided test molecules as input, please include a path containing files with the given input molecule"
            print "EX:  $ python get_combined_score.py -i <input molecule> -indir <directory_of_input_files>"
            exit()
#####################################################################################
def readtabfile_nextPARS(tabFile):
    '''Reads a tab file with the digestion information'''
    
    name, counts = tabFile.split()
    counts = [ float(x) for x in counts.strip(';').split(';')]
    
    return counts, name

#####################################################################################
def get_enzyme_profiles(flist, oldFileOrganization, tabnames, misalign, exp_indices, PARScalc, capto, last50, 
                        last50_matches, verbose, molname_full, low_coverage_flag, enz_name, molname, real_path=None):
    enz=[]
    included = 0
    tabcounter = 0
    for f in flist:
        if oldFileOrganization:
            data, filename = tsl.readtabfile(f, real_path, oldFileOrganization)
        else:
            if molname in ['TETp4p6', 'RDN5-1', 'RDN18-1', 'RDN25-1', 'RDN58-1']:
                data, filename = tsl.readtabfile(f, real_path)
            else:
                data, filename = readtabfile_nextPARS(f)
                filename = tabnames[tabcounter] #changing it so that it keeps the enzyme suffix
                tabcounter += 1
            
        if misalign:
            if exp_indices == 'all_indices':
                data = data
            else:
                data = [data[i] for i in exp_indices]
#        print data
        if PARScalc:
            a = data
            included += 1
        else:
            a, enzmean = tsl.normList_to_average(tsl.cap_to_percentile(data, capto), data, filename, 
                                               last50, last50_matches, verbose, molname_full)
    #        a, enzmean = tsl.exclude_underexpressed(tsl.cap_to_percentile(data, capto), filename)
            if len(a) > 1: #exclude empty lists returned when avg count per site too low
                included += 1
        #a=tsl.normList_to_average_log(tsl.cap_to_percentile(data, 90), filename)
        if len(a)==0:      ##For those files that had mean=0 and thus are discarded.
            enz = enz
        elif len(enz)==0:
            enz = a
        else:
            enz = [l+m for l,m in zip(enz,a)]   #Add corresponding values from each enz tab file in list file to be averaged in next step.
#        print a
    enz = [x/included for x in enz]   #Average of values from all enz tab files in list file
#    print colored(enz_name,'red'), enz
    try:
        enz_norm = tsl.normList_to1(enz)
    except ValueError:
        print colored('Insufficient average coverage in %s experiments' %enz_name,'red',attrs=['bold','underline'])
        low_coverage_flag = 1
        enz_norm = []
    
    return enz, enz_norm, low_coverage_flag, included
    
#####################################################################################
#####################################################################################

def main():

    #####################################################################################
    #####################################################################################
    parser = argparse.ArgumentParser(description="Get combined scores from tab files.\n",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-a', '--allthresh', dest = 'allthresh', action = 'store', default = None,
                        help = "Generate files with calculations to plot with plot_method_differences.py.")
                        
    parser.add_argument('-b', '--bpnucs', dest = 'bpnucs', action = 'store_true', default = False,
                        help = "Determine the particular nucleotides in base pairs.")
    
    parser.add_argument('-c', '--capto', dest = 'capto', action = 'store', type = int, default = 95,
                        help = "Cap to given percentile.")
    
    parser.add_argument('-f', '--fasta', dest = 'fasta', action = 'store', default = None,
                        help = "Path to fasta file for input molecule.")
    
    parser.add_argument('--full', dest = 'full', action = 'store_false', default = True,
                        help = "Use full length of molecules and real calls for plots and calculations.")
                        
    parser.add_argument('-i', '--input', dest = 'input', action = 'store', default = None,
                        help = "Input molecule name (EX: TETp4p6)")
    
    parser.add_argument('-inDir', '--inDir', dest = 'inDir', action = 'store', default = None,
                        help = "Directory containing input files for the given molecule")
    
    parser.add_argument('--ignore', dest = 'ignore', action = 'store_true', default = False,
                        help = "Ignore positions that are missing structural info or are part of pseudoknots.")
    
    parser.add_argument('-k', '--knots', dest = 'knots', action = 'store', default = 'remove_none',
                        help = "Collect and remove positions in pseudoknots or without structure info.")
    
    parser.add_argument('-l', '--last50', dest = 'last50', action = 'store_false', default = True,
                        help = "Remove last 50 positions because reads not available.")
    
    parser.add_argument('-m', '--misalign', dest = 'misalign', action = 'store_false', default = True,
                        help = "Read alignment file and exclude positions that do not match.")
    #                    action = 'store',
    #                    type = int,
    #                    nargs = '*',
    #                    default = None,
    #                    help = '''Remove positions that do not match reference in alignment. 
    #                        Indicate multiple positions by separating with a space ($ python ... -m 1 12 57)''')
    
    parser.add_argument('-n', '--norm', dest = 'norm', action = 'store', default = 'D',
                        help = "Normalization method.")
    
    parser.add_argument('-nr', '--normreads', dest = 'normreads', action = 'store', default = None,
                        help = "Tab file of normalized S1 (-nr {s/s1}) or V1 (-nr {v/v1}) read counts.")
    
    parser.add_argument('--nP_only', dest = 'nP_only', action = 'store', default = None,
                        help = 'Output file with pre-RNN nextPARS scores in tab file format.')
    
    parser.add_argument('-o', '--output', dest = 'output', action = 'store', default = None,
                        help = 'Output file with scores in tab file format.')
    
    parser.add_argument('-old', '--oldFileOrganization', dest = 'oldFileOrganization', action = 'store_true', default = False,
                        help = 'Use old organization of files, in which each transcript is stored in a single file.')
    
    parser.add_argument('-p', '--ppv', dest = 'ppv', action = 'store', default = None,
                        help = "Perform ACC calculations and write to Excel file.")
    
    parser.add_argument('-P', '--PARScalc', dest = 'PARScalc',action = 'store_true', default = False,
                        help = "Calculate scores as in original PARS paper (using log2(V-S)).")
    
    parser.add_argument('-r', '--removepos', dest = 'removepos', action = 'store', type = int, default = None,
                        help = "Create lists of cor and inc to remove consistent inc.")
    
    parser.add_argument('-s', '--spp', dest = 'spp', action = 'store_true', default = False,
                        help = 'Generate Structure Prefernece Profile file.')
    
    parser.add_argument('-ss', '--sppstats', dest = 'sppstats', action = 'store_true', default = False,
                        help = 'Output the stats about the spp files. ')
    
    parser.add_argument('-t', '--thresh', dest = 'thresh', action = 'store', type = float, default = 0.8,
                        help = 'Threshold for strict calls.')
    
    parser.add_argument('-T', '--threshvarna', dest = 'threshvarna', action = 'store', type = float, default = None,
                        help = 'Threshold for Varna output.')
    
    parser.add_argument('-v', '--verbose', dest = 'verbose', action = 'store_true', default = False,
                        help = "Verbose.")
    
    parser.add_argument('-V', '--varna', dest = 'varna', action = 'store', default = '',
                        help = 'Create colormap file from scores values. {nextPARS, spp, snorm, vnorm}')
    
    parser.add_argument('-w', '--window', dest = 'window', action = 'store', default = None,
                        help = "Summarize scores by average of windows around each position.")
    
    options = parser.parse_args()
    
    #####################################################################################
    #####################################################################################
    
#################### Obtain and manipulate data from V1 and S1 digestion files and real calls file (ref structure) ####################

    molname = options.input
    if options.inDir:
        exp_dir = options.inDir.rstrip('/')
    
    ###############################
    #collect digestion data
    try:
        real_calls = tsl.get_real_calls(molname)
    except IOError:
        real_calls = []  #if no ref structure available
    
    if options.oldFileOrganization:
        molname_full, real_path, v_exps, s_exps = get_tabfilelist(molname)
        V1,S1 = tsl.fileList(molname+'.txt')
        os.system('rm '+molname+'.txt')
    else:
        if molname in ['TETp4p6', 'RDN5-1', 'RDN18-1', 'RDN25-1', 'RDN58-1']:
            molname_full, real_path, v_exps, s_exps = get_tabfilelist_nextPARS(molname)
            V1,S1 = tsl.fileList(molname+'.txt')
            os.system('rm '+molname+'.txt')
        else:
            V1, S1, molname_full, v_exps, s_exps, real_path = get_tabfilelist_nextPARS(molname, exp_dir)
            
    
    
    ###############################
    

    ###############################
    #account for mismatches, gaps, etc
    if options.misalign: #FIXME: remove this option because must perform these steps always...
        if molname in ['RDN5-1', 'RDN25-1', 'RDN58-1']:
            superstruc = '3U5H'
        elif molname == 'RDN18-1':
            superstruc = '3U5F'
        elif molname == 'TETp4p6':
            superstruc = '1GID'
        else:
            superstruc = molname
            
        if options.ignore and os.path.exists('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.dp' %(data,superstruc)):
            #Use PDB output and ct to determine positions missing info or in internal/external pseudoknots
            dots = tsl.get_dot_string(molname)
            exts = tsl.get_external_knots(molname)
            bad_indices = tsl.get_bad_indices(dots, exts)
            ignore = tsl.indices_to_ignore(bad_indices, options.knots)
            exp_indices, match_in_ref, exp_fill_in, full_ref_seq, full_exp_seq, last50_matches = tsl.align_calls_ignore(molname_full, options.bpnucs, options.verbose, options.last50, knots=ignore)
            #FIXME: this section may now be deprecated, can erase??****
        else:
            exp_indices, match_in_ref, exp_fill_in, full_ref_seq, full_exp_seq, last50_matches = tsl.align_calls(molname_full, options.verbose, options.last50)

    else:
        #TODO: not this...also deprecated??******
        if options.last50:    #Removes last 50 positions by default to avoid bias due to lack of reads at these positions
            real_calls = real_calls[:-50]    #Retain by using command line option '-l' or '--last50'

    ###############################

    ####################
    low_coverage_flag = 0
    
    # get V1 and S1 profiles
    if options.verbose:
        print
        print colored("Tabfile \t Avg. counts per site", "cyan")
    V, V_norm, low_coverage_flag, included_v = get_enzyme_profiles(V1, options.oldFileOrganization, v_exps, options.misalign, exp_indices, 
                                                                   options.PARScalc, options.capto, options.last50, last50_matches, 
                                                                   options.verbose, molname_full, low_coverage_flag, "V1", molname, real_path)
    
    S, S_norm, low_coverage_flag, included_s = get_enzyme_profiles(S1, options.oldFileOrganization, s_exps, options.misalign, exp_indices, 
                                                                   options.PARScalc, options.capto, options.last50, last50_matches, 
                                                                   options.verbose, molname_full, low_coverage_flag, "S1", molname, real_path)

    ####################
    # exit if coverage is too low
    if low_coverage_flag == 1:
        print 'Try again later with better PARS data...Now exiting program'
        exit()
        
    ####################
    # compute combined score as the difference between normalized values of S1 and V1
    if options.PARScalc:
        combined = []
        for y,z in zip(S_norm,V_norm):
            if z == 0.0 and y == 0.0:
                combined.append(0.0)
            elif z == 0.0:
                combined.append('min') #no cuts by V1 => suggests unpaired
            elif y == 0.0:
                combined.append('max')  #no cuts by S1 => suggests paired
            else:
                combined.append(np.log2(z/y))
        combined_num = [x for x in combined if type(x) != str]
        combined = [max(combined_num) if x == 'max' else x for x in combined]
        combined = [min(combined_num) if x == 'min' else x for x in combined]
    else:
        # computes combined score as the difference between normalized values (after capping to percentile 90) of S1 and V1
        if options.norm.upper() == 'A':
            combined = tsl.normList_to1([z-y for y, z in zip(S,V)])
        elif options.norm.upper() == 'B':
            combined = [z-y for y,z in zip(S_norm, V_norm)]
        elif options.norm.upper() == 'C':
            combined = tsl.norm_pos_neg([z-y for y,z in zip(S,V)])
        elif options.norm.upper() == 'D':
            combined = [z-y for y,z in zip(S_norm, V_norm)]
            combined_norm = tsl.norm_pos_neg(combined)
            combined=combined_norm
    #combined=[safe_number((np.log10(a)/np.log10(b))) for a, b in zip(a,b)]
    #
#    print "Combined", combined, len(combined)
    
    
    if options.full and options.misalign:
#        print V_len, len(combined), combined
#        print len(exp_indices), len([x for x in exp_indices if x < V_len-50])
        combined_full = []
        append_cf = combined_full.append
        zer = []
        f = 0
        z = 0
        if full_ref_seq == 'all_indices':
            seq_len = len(combined)
            
        else:
            seq_len = len(full_ref_seq)
        
        for i in xrange(seq_len):
            if (type(match_in_ref)==list and i in match_in_ref) or match_in_ref == 'all_indices':
                append_cf(combined[f])
                f += 1
            elif i in exp_fill_in:
                append_cf(0.0)
                z+=1
                zer.append(i)
            else:
                print colored('bad %s' %i, 'red')
#        print combined_full, len(combined_full), z, f
#        print colored('zer', 'red'), zer, len(zer)
        combined = combined_full
    
    # to keep the nextPARS scores without incorporating the RNN score
    pre_combined = combined
        

#################### Incorporate the RNN model to finalize the scores ##############################    
    # First get fasta file for input molecule
    try:
        fas = "%s/SEQS/PROBES/%s.fa" %(data, molname)
        fas_temp = "%s_tmp" %fas
        fo = open(fas)
        fo.close()
        os.system("cp %s %s" %(fas, fas_temp))
    except IOError:
        fas = options.fasta
        fas_temp = "%s_tmp" %fas
        fo = open(fas)
        fw = open(fas_temp, "w")
        
        seq = ''
        for line in fo:
            l = line.strip()
            if l.startswith(">"):
                name = l.split(">")[1] 
            if name == molname and not l.startswith(">"):
                seq += l
        fw.write(">%s\n" %molname)
        fw.write(seq)
        
        fo.close()
        fw.close()
    
    # then get file with initial nextPARS scores
    score_tab = []
    for a in combined:
        #this is to avoid problems with exponentials => very low order values become 0
        if a<=0.0001 and a>=0.0:
            a=0.0
        elif a>=-0.0001 and a<=0.0:
            a=0.0
        score_tab.append(str(a))
    score_tab = '%s\t%s;' %(molname, ';'.join(score_tab))
#        print score_tab
    out = open("%s.tmp.out" %molname, 'w')
    out.write(score_tab)      ##Creates file with scores in tab format
    out.close()
    
    # Finally, run the RNN model calculation of scores
    os.system("python predict2.py -f %s -p %s.tmp.out -o %s.RNN.tab_tmp 2> /dev/null" %(fas_temp, molname, molname))
    os.system("rm %s.tmp.out" %molname)
    os.system("rm %s" %fas_temp)
    
    # Now can use those scores as the combined scores
    rnn = open("%s.RNN.tab_tmp" %molname)
    rnn_write = open("%s.RNN.tab" %molname, "w")
    combined = [ float(x) for x in rnn.readline().split()[-1].rstrip(';').split(';') ]
    combined_norm = tsl.norm_pos_neg(combined)
    combined=combined_norm
    rnn_write.write("%s\t%s;" %(molname, ';'.join([str(x) for x in combined])))
    rnn.close()
    rnn_write.close()
    os.system("rm %s.RNN.tab_tmp" %molname)
    
    
    if options.verbose:
        print
        print '%s V1 files and %s S1 files included' %(colored(str(included_v),'cyan'), colored(str(included_s),'cyan'))
        print colored("Combined", 'green'), combined, colored(len(combined), 'cyan')
        if full_ref_seq != "all_indices": #since these values will all be 0s because no CT file was available
            print colored("Real Calls", 'green'), real_calls, colored(len(real_calls), 'cyan')
            
        
#################### Strict calls and percentages ##############################
    if options.allthresh:   #For use with the script plot_best.py
        ACCs = []
        cors = []
        incs = []
        for t in np.linspace(0.0, 1.0, 101):
            strict_calls = tsl.only_reliable_calls(combined, t)
            ACC, PPV, sen, fpr, cor, inc, und, correct_x, correct_y, incorrect_x, incorrect_y, undet_x, undet_y = tsl.calculations(real_calls, strict_calls)
            ACCs.append(ACC)
            cors.append(cor)
            incs.append(inc)
            
#        print ACCs, len(ACCs)
#        print cors, len(cors)
#        print incs, len(incs)
        if options.allthresh == 'norm':
            facc = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_accs.txt' %(data,molname,molname,options.capto), 'a')
            fcor = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_cors.txt' %(data,molname,molname,options.capto), 'a')
            finc = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_incs.txt' %(data,molname,molname,options.capto), 'a')
        elif options.allthresh == 'cap':
            facc = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_accs.txt' %(data,molname,molname,options.norm), 'a')
            fcor = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_cors.txt' %(data,molname,molname,options.norm), 'a')
            finc = open('%s/TAB_FILES/REF_molecule_tab_files/ACC_plots/%s/%s_%s_incs.txt' %(data,molname,molname,options.norm), 'a')
            

        for a in ACCs:
            facc.write(str(a)+' ')
        facc.write('\n')
        facc.close()
        
        for c in cors:
            fcor.write(str(c)+' ')
        fcor.write('\n')
        fcor.close()
        
        for i in incs:
            finc.write(str(i)+' ')
        finc.write('\n')
        finc.close()
    
    else:
        strict_calls=tsl.only_reliable_calls(combined, options.thresh)
#        print strict_calls, len(strict_calls), len(real_calls)
#        c = 0
#        n = 0
#        u = 0
        
        if options.removepos == 2:    #For use after pipeline_removepos.py to remove bad positions
            corfile = open('../corpos_all.txt', 'r')    #TODO: figure out a use for this, if any
            incfile = open('../incpos_all.txt', 'r')
            for line in incfile:
                if line.startswith('>'+molname):
                    Is = line.split('\t')[1].strip()
                    Is = Is.split(' ')
            corfile.close()
            incfile.close()
            good_stricts = []
            good_reals = []
            for x in xrange(len(strict_calls)):
                if str(x) not in Is:
                    good_stricts.append(strict_calls[x])    #Positions excluding those incorrect under 80% of conditions
                    good_reals.append(real_calls[x])
            
            ACC, PPV, sen, fpr, cor, inc, und, correct_x, correct_y, incorrect_x, incorrect_y, undet_x, undet_y = tsl.calculations(good_reals, good_stricts)
        else:
            ACC, PPV, sen, fpr, cor, inc, und, correct_x, correct_y, incorrect_x, incorrect_y, undet_x, undet_y = tsl.calculations(real_calls, strict_calls)
            

        if options.verbose and full_ref_seq != "all_indices": #since these values will all be 0s because no CT file was available
            print colored('ACC', 'green'), '= %f%% of predicted values that are correct' %(ACC)
            print colored('PPV', 'green'), '= %f%% of predicted pairs that are actually in real calls' %(PPV)
            print colored('Correct', 'green'), '= %f%% of predicted values match real calls' %(cor)
#            print correct_y, colored("Correct " + str(len(correct_y)), 'cyan')
            print colored('  correct indices', 'cyan'), correct_x, colored("%s/%s correct" %(len(correct_x),len(strict_calls)), 'cyan')
            print colored('Incorrect', 'green'), '= %f%% of predicted values match real calls' %(inc)
#            print incorrect_y, colored("Incorrect " + str(len(incorrect_y)), 'cyan')
            print colored('  incorrect indices', 'cyan'), incorrect_x, colored("%s/%s incorrect" %(len(incorrect_x),len(strict_calls)), 'cyan')
            print colored('Undetermined', 'green'), '= %f%% of predicted values do not meet given threshold' %(und)
#            print undet_y, colored("Undetermined " + str(len(undet_y)), 'cyan')
            print colored('  undetermined indices', 'cyan'), undet_x, colored("%s/%s undetermined" %(len(undet_x),len(strict_calls)), 'cyan')
            print colored('Sensitivity', 'green'), '= %f%%' %sen
            print colored('FPR', 'green'), '= %.2f%%' %fpr

######################  Outputs  ####################################
    
    if options.ppv:
        from xlrd import open_workbook
        from xlutils.copy import copy
        
        if options.norm == 'A':
            if options.capto == 90: col = 1
            elif options.capto == 85:   col = 2
            elif options.capto == 75:   col = 3
        elif options.norm == 'B':
            if options.capto == 90: col = 5
            elif options.capto == 85:   col = 6
            elif options.capto == 75:   col = 7
        elif options.norm == 'C':
            if options.capto == 90: col = 9
            elif options.capto == 85:   col = 10
            elif options.capto == 75:   col = 11
        elif options.norm == 'D':
            if options.capto == 90: col = 13
            elif options.capto == 85:   col = 14
            elif options.capto == 75:   col = 15
        
        rbook = open_workbook(options.ppv, formatting_info=True)
        rppv = rbook.sheet_by_name('PPV')
        rcor = rbook.sheet_by_name('Correct')
#        rinc = rbook.sheet_by_name('Incorrect')
#        rund = rbook.sheet_by_name('Undetermined')
#        rsen = rbook.sheet_by_name('Sensitivity')
                
        wbook = copy(rbook)
        wppv = wbook.get_sheet(0)
        wcor = wbook.get_sheet(1)
        winc = wbook.get_sheet(2)
        wund = wbook.get_sheet(3)
        wsen = wbook.get_sheet(4)
        
        for t in np.linspace(0.1,0.8,15):
            if '%.2f' %t == '%.2f' %options.thresh:
                for row in xrange(rppv.nrows):
                    if rppv.cell_value(row,0) == 'T=%.2f' %options.thresh:
                        for row in xrange(row+1, row+6):
                            if rcor.cell_value(row, 0).startswith(molname):
                                wppv.write(row, col, '%.2f%%' %(ACC))
                                wcor.write(row, col, '%i (%.2f%%)' %(len(correct_y), cor))
                                winc.write(row, col, '%i (%.2f%%)' %(len(incorrect_y), inc))
                                wund.write(row, col, '%i (%.2f%%)' %(len(undet_y), und))
                                wsen.write(row, col, '%.2f%%' %sen)
        
                    
        if options.ppv.endswith('xlsx'):
            wbook.save(options.ppv[:-1])
        else:
            wbook.save(options.ppv)


    
    ##############################
    #printing to files
    ##############################   
    if options.spp:
        score_tab=[]
        append_st = score_tab.append
        for c in strict_calls:
            if c == -1.0:   #SPP files use 1 for SS, 0 for DS, NA for undetermined
                s = 1
            elif c == 1.0:
                s = 0
            elif c == 0.0:
                s = 'NA'
            append_st(str(s))
        score_tab = '%s\t%s;' %(molname_full, ';'.join(score_tab))
#        print score_tab
        spp_file = open(molname_full+'.spp','w')
        spp_file.write(score_tab)
        spp_file.close()
        
    elif options.sppstats:
#        statsfile = open('/home/jwillis/users/tg/jwillis/lncRNA_DATA/RNA2D/sppstats.txt','a')
        p = 0.0
        for c in strict_calls:
            if c in [-1.0, 1.0]:
                p += 1
#        print molname, p, len(strict_calls), p/len(strict_calls)
        sys.stdout.write('%s' %(100*p/len(strict_calls)))
#        statsfile.write('%s\t%s\t%s\n' %(molname, options.thresh, 100*p/len(strict_calls)))
#        statsfile.close()
        
    elif options.output:
        score_tab = []
        for a in combined:
            #this is to avoid problems with exponentials => very low order values become 0
            if a<=0.0001 and a>=0.0:
                a=0.0
            elif a>=-0.0001 and a<=0.0:
                a=0.0
            score_tab.append(str(a))
        score_tab = '%s\t%s;' %(molname_full, ';'.join(score_tab))
#        print score_tab
        out = open(options.output, 'w')
        out.write(score_tab)      ##Creates file with scores in tab format
        out.close()
    
    elif options.nP_only:
        score_tab = []
        for a in pre_combined:
            #this is to avoid problems with exponentials => very low order values become 0
            if a<=0.0001 and a>=0.0:
                a=0.0
            elif a>=-0.0001 and a<=0.0:
                a=0.0
            score_tab.append(str(a))
        score_tab = '%s\t%s;' %(molname_full, ';'.join(score_tab))
#        print score_tab
        out = open(options.nP_only, 'w')
        out.write(score_tab)      ##Creates file with scores in tab format
        out.close()
    
    elif options.normreads:
        if options.normreads.upper() in ['S','S1']:
            enz_norm = S_norm
            enz = 'S1'
        elif options.normreads.upper() in ['V','V1']:
            enz_norm = V_norm
            enz = 'V1'
        enz_norm = [str(round(x,4)) for x in enz_norm]
#        score_tab = molname+'\t'+';'.join(enz_norm)+';'
#        nr = open('/home/jwillis/lncRNA_DATA/RNA2D/%s_%s_norm_reads.tab' %(molname, enz),'w' )
#        nr.write(score_tab)
#        nr.close()
        tsl.colormap('%s/%s_%s_norm_reads.colormap' %(src, molname, enz), enz_norm, threshold=None)
    
    ##Create file with scores in colormap format
    if options.varna.upper() in ['C','COMBINED','NEXTPARS']:
        tsl.colormap('%s/%s_nextPARS.colormap' %(src,molname), combined, options.threshvarna)
    elif options.varna.upper() in ['S','SPP']:
        score_tab=[]
        append_st = score_tab.append
        for c in strict_calls:
            if c == -1.0:   #SPP files use 1 for SS, 0 for DS, NA for undetermined
                s = 1
            elif c == 1.0:
                s = 0
            elif c == 0.0:
                s = 'NA'
            append_st(str(s))
        print score_tab, molname
        tsl.colormap('%s/%s_spp.colormap' %(src,molname), score_tab, options.threshvarna)
    elif options.varna.upper() in ['SNORM','SN','S_NORM','S1']:
        enz_norm = [str(round(x,4)) for x in S_norm]
        tsl.colormap('%s/%s_S1_norm_reads.colormap' %(src,molname), enz_norm, threshold=None)
    elif options.varna.upper() in ['VNORM','VN','V_NORM','V1']:
        enz_norm = [str(round(x,4)) for x in V_norm]
        tsl.colormap('%s/%s_V1_norm_reads.colormap' %(src,molname), enz_norm, threshold=None)
    
    
    
    
    ##############################
    if options.removepos == 1:  #For use with script pipeline_removepos.py
        corfile = open('corpos.txt', 'a')
        incfile = open('incpos.txt', 'a')
        removepos(combined, correct_x, incorrect_x, corfile, incfile)
        corfile.close()
        incfile.close()
    
    elif options.removepos == 2:    #For use after pipeline_removepos.py to remove bad positions
        corfile = open('../corpos_all.txt', 'r')
        incfile = open('../incpos_all.txt', 'r')
        for line in incfile:
            if line.startswith('>'+molname):
                Is = line.split('\t')[1].strip()
                Is = Is.split(' ')
        print 'Is', Is, type(Is), type(Is[0])
        corfile.close()
        incfile.close()
        #TODO: use Is to ignore these positions when calculating percents => implement in align_calls function(??)
        #Will also have a options.removepos == 3 for completely removing before calculating combined scores
    ##############################
    
    
    ##############################
    if options.window:
        V_win = tsl.windows(V_norm, int(options.window))
        V_win_norm = tsl.normList_to1(V_win)
#        print 'V_win:', V_win, '\n'
        S_win = tsl.windows(S_norm, int(options.window))
        S_win_norm = tsl.normList_to1(S_win)
#        print 'S_win:', S_win, '\n'
#        win = tsl.windows(combined, int(options.window))
#        for x,y,z in zip(V_win, S_win, combined):
        for x,y,z in zip(V_win_norm, S_win_norm, combined):
            if z >= options.thresh or z <= options.thresh*(-1.0):
                print (x,y,z)
            else:
                print 'Not above threshold'


######################  Determine nucleotides  ####################################
    if options.bpnucs:
        print colored('Predicted paired nucs:', 'yellow')
        strict_nucs = bp_nucs(strict_calls, full_exp_seq)
        
        print colored('Real paired nucs:', 'yellow')
        real_nucs = bp_nucs(real_calls, full_exp_seq)
        
        paired_nucs = []
        for n in strict_nucs:
            if n in real_nucs:
                paired_nucs.append(n)
        print colored('In both','green'), paired_nucs, len(paired_nucs), '\n'
        
        pairs_ct = pairs_from_ref_ct(molname)
        
        correct_bps(paired_nucs, pairs_ct)


#####################################################################################
####################### Main function call ##########################################
#####################################################################################

if __name__ == "__main__":
    main()

