import os
import numpy as np
from Bio import AlignIO
from termcolor import colored


# global variables
#home = '/home/jwillis/users/tg'
#home = '/users/tg'
#data = '%s/jwillis/lncRNA_DATA' %home
#src  = '%s/jwillis/lncRNA_DATA/RNA2D' %home


home = '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-1])
data = '%s/data' %home
sbin = '%s/bin' %home


#####################################################################################
#####################################################################################
def get_real_calls(molname):
    '''Read ct file to return list of real calls'''
    
#    ct_path = '%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/STRUCTURES/REFERENCE_STRUCTURES/' %home
    ct_path = '%s/STRUCTURES/REFERENCE_STRUCTURES/' %data
    if molname in ['RDN18-1']:
        ct = open(ct_path+molname.replace('-1','')+'-PDB.ct', 'r') #replace since file called RDN18-PDB in folder
    elif molname in ['B2', 'U1', 'SRA']:
        ct = open(ct_path+molname+'-ref.ct', 'r')
    else:
        ct = open(ct_path+molname+'.ct', 'r')
        
    reals = []
    append_r = reals.append
    nucs = ['A', 'C', 'G', 'U', 'T', 'X', 'a', 'c', 'g', 'u', 't', 'x']  #To ensure use of only those lines containing the position and pairing info
    
    for line in ct:
        if not line.startswith('#'):     #Ignore commented lines at beginning
            pair = line.split()
            if len(pair) > 1 and pair[1] in nucs:    #Ignore empty lines and first line with number of positions and additional info
                if pair[4] == '0':                #Pairing is found in 5th column of ct files
                    append_r('0')          #When position is SS
                else:
                    append_r('1')          #When position is DS
    return [float(r) for r in reals]
#####################################################################################
def fileList(filename, T1_A=False):
    '''Reads a file with the tabfiles list'''
    listfile=open(filename,'r')
    lines=listfile.readlines()
    listfile.close()
    
    if T1_A:
        enz = []
        for line in lines:
            enz.append(line.strip())
        return enz
        
    else:
        V1=[]
        S1=[]
        for line in lines:
    #        test=string.split(line,'_')
    #        #print test
    #        
    #        # this adds the path to where we store all data in the cluster
    #        # note, shall we make this optional?
    #        #path='/users/tg/tgabaldon/PROJECTS/NONCODEVOL/DATA/TAB_FILES/our_data/'+test[1]+'/'
    #        #path = '/Users/owner1/lncRNA_DATA/Scripts/'+test[1]+'/'  #Temporary for my computer
    #        path = '/Users/owner1/lncRNA_DATA/Scripts/'
    #        #path = '/Users/owner1/lncRNA_DATA/TAB_FILES/our_data/0000-00-00c/'
    #        #if options.infile:     #For if it is necessary to indicate.......???? 
    #        
            if 'V1' in line:
                V1.append(line.strip())
            if 'S1' in line:
                S1.append(line.strip())
        return V1, S1
    
#####################################################################################
def readtabfile(filename, path, oldFileOrganization=False):
    '''Reads a tab file with the digestion information'''

    if oldFileOrganization:
        experiment = filename.split('_')[-2] #experiment name, EX: 2014-04-29c
        tabfile=open(path+'/'+experiment+'/'+filename,'r')
    else:
        tabfile = open('%s/%s' %(path,filename), 'r')
    lines=tabfile.readlines()
    tabfile.close()
#    print lines, len(lines)
    dat=[]
    test=lines[0].strip().split('\t')   #strip() removes '\n'
#    print 'test', test
    dat=map(float,test[1].strip(';').split(';'))
#    print 'dat', dat
    return dat, filename
#####################################################################################
#####################################################################################

def get_dot_string(molname):
    '''Extract dot-bracket sequence from dp file'''
    
    if molname in ['RDN5-1', 'RDN25-1', 'RDN58-1']:
        dp_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/3U5H.dp' %data, 'r')
    elif molname == 'RDN18-1':
        dp_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/3U5F.dp' %data, 'r')
    elif molname == 'TETp4p6':
        dp_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/TETp4p6.dp' %data, 'r')
    else:
        dp_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.dp' %(data,molname), 'r')
        
    dots = ['.', '(', ')', '[', ']', ' {', '}', '<', '-'] #exclude '>' because header lines begin with it
    dot_list = []
    for line in dp_file:
        for line in dp_file:
            if 1 in [d in line for d in dots]:  #take only lines with dot-bracket symbols
                dot_list.append(line.strip())
    dp_file.close()

    seq = ''.join(dot_list)
    if molname == 'RDN25-1':
        dot_seq = seq[0:3395]
    elif molname == 'RDN5-1':
        dot_seq = seq[3396:3516]
    elif molname == 'RDN58-1':
        dot_seq = seq[3517:]
    elif molname == 'TETp4p6':
        dot_seq = seq[0:158]
    else:
        dot_seq = seq
    return dot_seq
#####################################################################
def get_external_knots(molname):
    '''For finding positions in one strand that bind to positions in a separate strand.'''
    
    if molname in ['RDN5-1', 'RDN25-1', 'RDN58-1']:
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/3U5H.ct' %data, 'r')
    elif molname == 'RDN18-1':
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/3U5F.ct' %data, 'r')
    elif molname == 'TETp4p6':
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/TETp4p6.ct' %data, 'r')
#        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/TETp4p6-PDB.ct' %data, 'r')
    else:
        ct_file = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.ct' %(data,molname), 'r')
        
    nucs = ['A', 'C', 'G', 'U', 'T', 'X', 'a', 'c', 'g', 'u', 't', 'x']
    cts = []
    append_cts = cts.append
    for i, line in enumerate(ct_file): #May be no good since it will also count any header lines for i
        line = [int(x) if x not in ['A','G','T','C','U'] else x for x in line.split()]
        if len(line) > 1 and line[1] in nucs:
            #Collecting all positions for given molecule from ct file containing multiple molecules
            #FIXME: not ideal efficiency since must manually look into ct files to know where particular molecules begin
            if molname == 'RDN25-1':
                if line[0]-line[5] == 0:
                    append_cts((line[0], line[4], line[5]))  #(pos in ct, pos of pair, pos in strand)
            elif molname == 'RDN5-1':
                if line[0]-line[5] == 3396:
                    append_cts((line[0], line[4], line[5]))
            elif molname == 'RDN58-1':
                if line[0]-line[5] == 3517:
                    append_cts((line[0], line[4], line[5]))
            elif molname == 'TETp4p6':
                if line[0] <= 158:
                    append_cts((line[0], line[4], line[5]))
            else:
                append_cts((line[0], line[4], line[5]))
    ct_file.close()
    
    externals = []
    append_ext = externals.append
    for c in cts:
        if c[1] != 0:   #Only looking for pairs
            if c[1]<cts[0][0] or c[1]>cts[len(cts)-1][0]:   #When position of pair is not within the given strand
                append_ext((c[0],c[1],c[2],c[2]-1))  #externals[3] is the index for the given strand
    return externals
####################################################################
def get_bad_indices(dots, exts):
    '''Determine positions in ext/int pseudoknots or missing info from PDB'''
    bad_indices = []
    append_bi = bad_indices.append
    remove_bi = bad_indices.remove
    i=0
    for d in list(dots):
        if d in ['-', '[', ']', '{', '}', '<', '>', 'A', 'a', 'B', 'b']:    #'-' are positions with no pairing info, others are psuedos
            if d == '-':
                append_bi(str(i)+'m')  #missing info
            else:
                append_bi(str(i)+'i')  #internal knot
        i+=1
    #Add external knots to bad_indices, replacing any that were incorrectly labeled internal due to '[','{', etc.
    for i in exts:
        if i[3] not in [int(b.replace('m','').replace('i','').replace('e','')) for b in bad_indices]:
            append_bi(str(i[3])+'e')
        else:
            remove_bi(str(i[3])+'i')   #For any positions that are actually between molecules, but appear as '[]' or '{}', etc.
            append_bi(str(i[3])+'e')
    return sorted(bad_indices, key=lambda x:int(x.replace('m','').replace('i','').replace('e','')))
#####################################################################################
def indices_to_ignore(bad_indices, knots):
    '''Read file containing indices for positions to be ignored because missing structural info or in pseudoknots'''
    
    if knots == 'all':
        ignore = [n.replace('m','').replace('i','').replace('e','') for n in bad_indices]
    elif knots == 'external':
        ignore = [n for n in bad_indices if n.endswith('m') or n.endswith('e')]
        ignore = [b.replace('m','').replace('e','') for b in ignore]
    elif knots == 'internal':
        ignore = [n for n in bad_indices if n.endswith('m') or n.endswith('i')]
        ignore = [b.replace('m','').replace('i','') for b in ignore]
    else:
        ignore = [n for n in bad_indices if n.endswith('m')]
        ignore = [b.replace('m','') for b in ignore]
    return [int(b) for b in ignore]
#####################################################################################

def align_calls_ignore(molname, bpnucs, verbose, last50, knots=None):
    '''Returns a list of indices only for those positions aligned in clw file and neither in pseudoknots nor missing info in PDB files,
        as well as the consensus sequence and a list of exp_indices that will not be considered with options.full'''

    try:
#        clustal_file = open('%s/tgabaldon/PROJECTS/NONCODEVOL/DATA/SEQS/CLW/%s.clw' %(data,molname), 'r')
        clustal_file = open('%s/SEQS/CLW/%s.clw' %(data,molname), 'r')
        refseq = ''
        expseq = ''
        for line in clustal_file:
            if line.startswith(molname+'-ref') or line.startswith(molname+'-PDB') or line.startswith(molname+'-3U5H') or line.startswith(molname+'-RNASTRAND'):
                refseq = refseq+line.split()[1]
            elif line.startswith(molname+'-1') or line.startswith(molname+'\t') or line.startswith(molname+' '):
                expseq = expseq+line.split()[1]
        clustal_file.close()
        
        if bpnucs and not verbose:
            print colored("refseq", 'green'), refseq, colored(len(refseq), 'cyan')
            print colored("expseq", 'green'), expseq, colored(len(expseq), 'cyan')
                        
        #collect list of nucleotides with position in reference, determine if it should be a knot to be ignored
        n = 0
        total = 0
        reflist = []
        append_ref = reflist.append
        actual_refseq = ''
        for a in refseq:
            if a != '-':
                actual_refseq += a
                if knots and n in knots:
                    append_ref((a,n,total))   #total gives position in the alignment sequences (includes gaps) of pseudos to be ignored
                else:
                    append_ref((a,n))   #n gives position of nuc within sequence
                n += 1
            else:
                append_ref(a)
            total+=1
    #    print reflist, len(reflist)
        
        #collect list of nucleotides with position in experimental sequence
        m = 0
        explist = []
        append_exp = explist.append
        actual_expseq = ''
        for b in expseq:
            if b != '-':
                actual_expseq += b
                append_exp((b,m))   #m gives position of nuc within sequence
                m += 1
            else:
                append_exp(b)
    #    print explist, len(explist)
    
        ref_indices = []
        exp_indices = []
        not_in_full = []
        append_refi = ref_indices.append
        append_expi = exp_indices.append
        append_noti = not_in_full.append
        for i in xrange(len(refseq)):
            if refseq[i] == expseq[i] or (refseq[i] == 'U' and expseq[i] == 'T') or (refseq[i] == 'T' and expseq[i] == 'U'):
                if len(reflist[i]) <= 2:    #Take only those that match in alignment and are not in list of ignored positions
                    append_refi(reflist[i][1])
                    append_expi(explist[i][1])
            elif refseq[i] == '-':
                append_noti(explist[i][1])
    
        if last50:
            ref_indices = [x for x in ref_indices if x < len(actual_refseq)-50] #Removes last 50 positions by default to avoid bias in calculations due to lack of reads at these positions
            exp_indices = [x for x in exp_indices if x < len(actual_expseq)-50] #Retain by using command line option '-l' or '--last50'
    
        if verbose:
            print colored("refseq", 'green'), refseq, colored(len(refseq), 'cyan')
            print colored("expseq", 'green'), expseq, colored(len(expseq), 'cyan')
            print colored("ref_indices", 'green'), ref_indices, colored(len(ref_indices), 'cyan')
            print colored("exp_indices", 'green'), exp_indices, colored(len(exp_indices), 'cyan')
        mol_seq = refseq.replace('-','')
    except IOError:
        try:
            fasta_file = open('%s/SEQS/PROBES/%s.fa' %(data,molname), 'r')
        except IOError:
            try:
                fasta_file = open('%s/SEQS/PROBES/%s.fasta' %(data,molname), 'r')
            except IOError:
                fasta_file = open('%s/SEQS/PROBES/%s' %(data,molname), 'r')
        mol_seq = [line.strip() for line in fasta_file.readlines() if not line.startswith('>')]
        fasta_file.close()
        mol_seq = ''.join(mol_seq)
        ref_indices = [x for x in xrange(len(mol_seq))]
        if last50:
            ref_indices = ref_indices[:-50]
        exp_indices = ref_indices
        not_in_full = []
    
        if verbose:
            print colored("sequence", 'green'), mol_seq, colored(len(mol_seq), 'cyan')
        
    return ref_indices, exp_indices, mol_seq, not_in_full
#####################################################################################
def align_calls(molname, verbose, last50):
    '''Returns a list of indices only for those positions aligned in clw file and neither in pseudoknots nor missing info in PDB files,
        as well as the consensus sequence and a list of exp_indices that will not be considered with options.full'''

    try:
        try:
            clustal_file = AlignIO.read('%s/SEQS/CLW/%s.clw' %(data,molname), "clustal")
        except AssertionError:
            raise IOError
        if clustal_file[0].id.upper().endswith(('REF','PDB', '3U5H', 'RNASTRAND')):
            refseq = str(clustal_file[0].seq)
            expseq = str(clustal_file[1].seq)
        else:
            expseq = str(clustal_file[0].seq)
            refseq = str(clustal_file[1].seq)
              
        match_in_ref, exp_indices = [], []
        exp_fill_in = []
        ref_nuc, exp_nuc = 0, 0
        full_ref_seq = ''
        full_exp_seq = ''
        for i in xrange(len(refseq)):
            if refseq[i] != '-':
                full_ref_seq += refseq[i]
                if refseq[i] == expseq[i] or (refseq[i] == 'U' and expseq[i] == 'T') or (refseq[i] == 'T' and expseq[i] == 'U'):
                    #indices that will be given scores from experiment
                    match_in_ref.append(ref_nuc)
                else:
                    #will appened scores of 0.0 later to these positions for comparison to the ref structure
                    exp_fill_in.append(ref_nuc)
                ref_nuc += 1
                
            if expseq[i] != '-':
                full_exp_seq += expseq[i]
                if expseq[i] == refseq[i] or (expseq[i] == 'U' and refseq[i] == 'T') or (expseq[i] == 'T' and refseq[i] == 'U'):
                    #those indices in the exp sequence that will be scored and compared to the reference
                    #if not in this list, will be ignored because cannot compare to 
                    exp_indices.append(exp_nuc)
                exp_nuc += 1

        #last 50 have no info in experiment so must ignore them when calculating average later
        #so must remove those sites in the last 50 that are actually scored and compared to reference
        last50_matches = len([x for x in exp_indices if x >= len(full_exp_seq)-50])
#        print 'last50 matches', last50_matches
        
        
#        newref_actual = [refseq.replace('-','')[x] for x in match_in_ref]
#        newexp_actual = [expseq.replace('-','')[x] for x in exp_indices]
#
#        print
#        print ''.join(newref_actual), len(newref_actual)
#        print ''.join(newexp_actual), len(newexp_actual)
#        print
#        print full_ref_seq, colored(len(full_ref_seq), 'cyan')
#        print full_exp_seq, colored(len(full_exp_seq), 'cyan')
#        print

        if verbose:
            print colored("refseq", 'green'), refseq, colored(len(refseq), 'cyan')
            print colored("expseq", 'green'), expseq, colored(len(expseq), 'cyan')
            print colored("match_in_ref", 'green'), match_in_ref, colored(len(match_in_ref), 'cyan')
            print colored("exp_indices", 'green'), exp_indices, colored(len(exp_indices), 'cyan')
            print

    except IOError:
        try:
            fasta_file = open('%s/SEQS/PROBES/%s.fa' %(data,molname), 'r')
        except IOError:
            try:
                fasta_file = open('%s/SEQS/PROBES/%s.fasta' %(data,molname), 'r')
            except IOError:
                try:
                    fasta_file = open('%s/SEQS/PROBES/%s' %(data,molname), 'r')
                except IOError: # create fasta file from sequence found in ct structure file if necessary
                    try:
                        nucs = ['A', 'C', 'G', 'U', 'T', 'X', 'a', 'c', 'g', 'u', 't', 'x']  #To ensure use of only those lines containing the position and pairing info
                        ct = open('%s/STRUCTURES/REFERENCE_STRUCTURES/%s.ct' %(data,molname), 'r')
                        ctl = ct.readlines()
                        ct.close()
                        seq = [line.split()[1].upper() for line in ctl if not line.startswith('#') and len(line.split())>1 and line.split()[1] in nucs]
                        fasta_write = open('%s/SEQS/PROBES/%s.fa' %(data,molname), 'w')
                        fasta_write.write('>%s\n%s' %(molname, ''.join(seq)))
                        fasta_write.close()
                        fasta_file = open('%s/SEQS/PROBES/%s.fa' %(data,molname), 'r')
                    except IOError: # when no fasta or ct files available, will simply include all indices
                        exp_indices = 'all_indices'
                        match_in_ref = 'all_indices'
                        exp_fill_in = []
                        full_ref_seq = 'all_indices'
                        full_exp_seq = 'all_indices'
                        last50_matches = 50
                        
                        return exp_indices, match_in_ref, exp_fill_in, full_ref_seq, full_exp_seq, last50_matches
                        
                                
        full_exp_seq = [line.strip() for line in fasta_file.readlines() if not line.startswith('>')]
        fasta_file.close()
        full_exp_seq = ''.join(full_exp_seq)
        full_ref_seq = full_exp_seq
        
        exp_indices = [x for x in xrange(len(full_exp_seq))]
        match_in_ref = exp_indices
        exp_fill_in = []
        last50_matches = 50
    
        if verbose:
            print colored("sequence", 'green'), full_exp_seq, colored(len(full_exp_seq), 'cyan')
        
    return exp_indices, match_in_ref, exp_fill_in, full_ref_seq, full_exp_seq, last50_matches

#####################################################################################
#####################################################################################
def normList_to_average(L, uncapped_L, filename, last50, last50_matches, verbose, molname):
    '''Normalize values of a list to average'''
    
    #verify that average read count per site is at least 5
    if last50:
        mean=np.mean(uncapped_L[:-last50_matches])
    else:
        mean=np.mean(uncapped_L)
        
    if verbose:
        if mean < 5:
            print "%s \t %f \t" %(filename, mean), colored('*Not included', 'red')
        else:
            print "%s \t %f" %(filename, mean)
#    return [x/mean for x in L if mean >= 1], mean   ##Since some tab files have mostly 0s and thus mean = 0, do not include those
    if mean >= 5:# or molname == '226_neg_ROX2':
        return [x/mean for x in L], mean   ##Mean >= 1 since some tab files have mostly 0s and thus mean = 0, include only files with average call per position of at least 1
    else:
        return [], mean
#####################################################################################
def normList_to1(L, normalizeTo=1.0):
    '''Normalize values of a list to make its max = normalizeTo'''
    vMax = max(L)
    return [x/(vMax*1.0)*normalizeTo for x in L]
#####################################################################################
def norm_pos_neg(L, normalizeTo=1):
    '''Normalize positive values and negative values separately to 1'''
    pos = []
    neg = []
    zer = []
    append_p = pos.append
    append_n = neg.append
    append_z = zer.append
    i = 0
    for a in L:
        if a > 0:
            append_p((a, i))
        elif a < 0:
            append_n((a, i))
        elif a == 0:
            append_z((a, i))
        i += 1
            
    pmax = max(pos)[0]
    nmin = min(neg)[0]
    
    pnorm = [(p/(pmax*1.0)*normalizeTo, j) for (p, j) in pos]
    nnorm = [(n/(nmin*-1.0)*normalizeTo, j) for (n, j) in neg]
    
    normlist = pnorm + nnorm + zer    
    normlist = sorted(normlist, key=lambda x: x[1])    
    return [x[0] for x in normlist]
    
#####################################################################################
def exclude_underexpressed(L, filename, verbose):
    '''Exclude digestion files with average calls under 1'''
    mean=np.mean(L)
    if verbose:
        if mean < 1:
            print "%s \t %f \t *Not included" %(filename, mean)
        else:
            print "%s \t %f" %(filename, mean)
    return [x for x in L if mean >= 1], mean
#####################################################################################
def normList_to_average_log(L, filename, verbose):
    '''Normalize values of a list to average and get the log2'''
    mean=np.mean(L)
    if verbose:
        print "%s \t" %(filename), mean
#    l=[]
#    print l
#    for x in L: 
#        if mean != 0:
#            if x != 0:
#                l.append(np.log2(x/mean))
#            else:
#                l.append(0)
#    print l
#    return l
    if mean != 0:
        return [np.log2(x/mean) if x != 0 else 0 for x in L]     #Force value of 0 for positions with 0 in tab file, because cannot take log of 0
    else:
        return []          ##Since some tab files have mostly 0s and thus mean = 0, do not include those

    
#####################################################################################
#####################################################################################

def cap_to_percentile(L, cap):
    '''Cap values to the given percentile'''
    cap=np.percentile(L,cap,axis=None, out=None, overwrite_input=False)
    capped=[]
    append_cap = capped.append
    for a in L:
        if a>=cap:
            append_cap(cap)
        else:
            append_cap(a)
    
    return capped

#####################################################################################
#####################################################################################

def windows(combined, window):
    '''Return list of average scores in given window for each position.'''
    win_avg = []
    for i in xrange(len(combined)):
        avg = []
        for j in xrange(-window, window+1):   #Uses j adjacent positions
            if len(combined) >= i+j+1:
#                print i+j
                avg.append(combined[i+j])       #At beginning of sequence, will wrap around to end of sequence for -j positions
            else:
#                print 'not'                
                avg.append(combined[window-j])   #When at end of list, wrap back around to beginning of sequence
#            print avg
        win_avg.append(np.mean(avg))
#    print win_avg
    return win_avg

#####################################################################################
#####################################################################################

def calculations(real, strict):
    '''Perform ACC, other percentage calculations for combined scores'''
    tp = tn = fp = fn = u =0
    correct_x = []
    incorrect_x = []
    undet_x = []
    append_c = correct_x.append
    append_i = incorrect_x.append
    append_u = undet_x.append
    for i in xrange(len(real) if len(real) < len(strict) else len(strict)):
        if (strict[i] == 1.0 and real[i] == 1):# or (strict[i] == -1.0 and real[i] == 0):
            tp += 1
            append_c(i)
        elif (strict[i] == -1.0 and real[i] == 0):
            tn += 1
            append_c(i)
        elif (strict[i] == 1.0 and real[i] == 0):# or (strict[i] == -1.0 and real[i] == 1):
            fp += 1
            append_i(i)
        elif (strict[i] == -1.0 and real[i] == 1):
            fn += 1
            append_i(i)
        else:
            u += 1
            append_u(i)

#    print tp, tn, fp, fn, u

    correct_y = []
    incorrect_y = []
    undet_y = []
    append_cy = correct_y.append
    append_iy = incorrect_y.append
    append_uy = undet_y.append

    ACC = (float(tp+tn)/float(tp+fp+tn+fn))*100 if tp+fp+tn+fn!=0 else 0.0 #if no ref structure available, these values will all be 0
    PPV = (float(tp)/float(tp+fp))*100 if tp+fp!=0 else 0.0
    sen = (float(tp)/float(tp+fn))*100 if tp+fn!=0 else 0.0
    fpr = (float(fp)/float(tn+fp))*100 if tn+fp!=0 else 0.0
    for i in correct_x:
        append_cy(real[i]/2) #Divide real by 2 because will be used in plot to show location of correct calls
    cor = 100*(float(tp+tn)/float(len(strict)))   ## % of positions in PARS prediction that are correct
    for i in incorrect_x:
        append_iy(real[i]/2) #Divide real by 2 because will be used in plot to show location of incorrect calls
    inc = 100*(float(fp+fn)/float(len(strict)))
    for i in undet_x:
        append_uy(real[i])
    und = 100*(float(u)/float(len(strict)))
    
    return ACC, PPV, sen, fpr, cor, inc, und, correct_x, correct_y, incorrect_x, incorrect_y, undet_x, undet_y
#####################################################################################

def only_reliable_calls (L,threshold):
    """Gives +1 for DS above threshold and -1 for SS below negative threshold"""
    
    negative_threshold=(threshold)*(-1)
    strict_calls=[]
    append_sc = strict_calls.append
    for a in L:
        if a>=threshold:           #When value is positive, looking for DS positions
            call = 1.0
        elif a<=negative_threshold:  #When value is negative, looking for SS positions
            call = -1.0
        else:
            call = 0.0              #For values that do not meet pos or neg threshold
        append_sc(call)
   
    return strict_calls
#####################################################################################
#####################################################################################

def sensitivity_falsepositiverate (L, real, threshold,positive='1.0'):
    """Looks for false positives of DS calls"""
    
    if positive=='1.0':
        negative='0.0'
    else:
        negative='1.0'

    TP=0.0
    FN=0.0
    TN=0.0
    FP=0.0
    if len(real) > 0: #only if ref structure was available to determine real calls
        index=0
        for a in L:
    #        print a,threshold
            if a>=threshold:
                call=positive
    #            print a,threshold,call,real_calls[index]
                if call==str(real[index]):    #So here (when positive = '1.0') a TP is a correctly called DS position 
                    TP=TP+1.0
                else:
                    FP=FP+1        #FP is a falsely called DS position
            elif a<=threshold*(-1.0):   #else:
                call=negative   #TODO: (^) was originally just else:, but should it be this elif(???)
    #            print a,threshold,call,real_calls[index]
                if call==str(real[index]):    #A TN would be a correctly called SS position
                    TN=TN+1.0
                else:
                    FN=FN+1        #FN is a falsely called SS position
            index=index+1
                
#    print TP,FN,TN,FP
    if FN+TP==0.0:
        sensitivity=0.0
    else:
        sensitivity=(TP/(TP+FN))*100
    if TN+FP==0.0:
        falsepositiverate=0.0
    else:
        falsepositiverate=((FP/(TN+FP))*100)
    #print threshold,sensitivity, falsepositiverate,TP,FP,TN,FN
    return sensitivity,falsepositiverate
    
#####################################################################################
#####################################################################################

def colormap(varnafile, scores, threshold=None):
    '''Create file with scores in format readable for Varna.'''
    varna_file = open(varnafile, 'w')
    if threshold:
        for i in scores:
            if i < threshold and i > threshold*(-1.0):
                i = 0.0
            varna_file.write(str(i)+'\n')
    else:
        for i in scores:
            if i <= 0.0001 and i >= -0.0001:
                i = 0.0
            varna_file.write(str(i)+'\n')
    varna_file.close()
    
#####################################################################################
#####################################################################################