
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  9 11:31:15 2017

@author: ahafez
"""

import numpy as np
import pandas as pd

import os
import sys
scriptpath = os.path.dirname(os.path.abspath(__file__))
#sys.path.append(os.path.abspath(scriptpath))
#sys.stderr = os.devnull

import argparse

#%%
parser = argparse.ArgumentParser(description='RNN Classifier for RNA sequences .')
parser.add_argument('-f' ,  dest='fastaFilename', metavar='fastaFilename', type=str, help="Input fasta file" , required=True)
parser.add_argument('-p' ,  dest='parScoreFileName',metavar='parScoreFileName', type=str, help="Input Pars Score tab file", required=True)
parser.add_argument('-w1' ,  dest='w1',metavar='w1', type=float, help="Weight for RNN score. Default 0.5", default=0.5)
parser.add_argument('-w2' ,  dest='w2',metavar='w2', type=float, help="Weight for Pars score. Default 0.5", default=0.5)

parser.add_argument('-o' , dest="outfileName",metavar="outputfileName",default="scores.tab", type=str, help="Final Score tab file.")

args = parser.parse_args()

outputPath = os.path.dirname(args.outfileName)

outputScoreTabFileName = args.outfileName

fastaFileName = args.fastaFilename

scoresTabFileName = args.parScoreFileName
#%%
import keras.models

__Key_SeqName_ID = "name"
# for getting seqeunce infomrtion if any is with the title
def seqName_cuts(line,clm_Names=None,sep=" "):
    dic = {}
    tokens = line.split(sep)
    for i in range(0,len(tokens)):
        keyName = str(i)
        if i == 0 :
            keyName = __Key_SeqName_ID
        if clm_Names is not None:
            if i < len(clm_Names):
                keyName = clm_Names[i]
            else:
                 keyName = str(i-len(clm_Names))
        if "=" in tokens[i]:
            kv = tokens[i].split("=")
            tokens[i] = kv[1]
            keyName = kv[0]
        dic[keyName] = tokens[i].rstrip()
    return dic
# read fasta file 
def readFasta_file(inputfilename, clm_Names=None,sep=" "):
    inputFile = open(inputfilename, "r")
    lines = inputFile.readlines()
    
    seqs  = {}
    seqs_data  = {}

    allSeq = {}
    allSeqData = {}
    seqLines = ""
    seqName = ""
    nameKeyIDStr = __Key_SeqName_ID
    if clm_Names is not None:
        nameKeyIDStr = clm_Names[0]
    for line in lines:
        if line.startswith(">"):
            if seqName is not "":
                seqs["id"] = seqName.rstrip()
                seqs_data["data"] = np.array(list(seqLines),dtype=np.str)
                seqs_data["id"] = seqName.rstrip()
                metaInfo = seqName_cuts(seqName[1:],clm_Names,sep)
                seqs.update(metaInfo)
                seqs_data.update(metaInfo)
                seqLines = ""
                
                allSeq[seqs[nameKeyIDStr]] = seqs
                allSeqData[seqs[nameKeyIDStr]] = seqs_data
            seqs = {}
            seqs_data= {}
            seqName = line
        else:
            seqLines += line.rstrip()
            #seqData.extend(list(line.rstrip()))
            
    seqs = {}
    seqs_data= {}

    seqs["id"] = seqName
    
    seqs["id"] = seqName.rstrip()
    seqs_data["data"] = np.array(list(seqLines),dtype=np.str)
    seqs_data["id"] = seqName.rstrip()
    metaInfo = seqName_cuts(seqName[1:],clm_Names,sep)
    seqs.update(metaInfo)
    seqs_data.update(metaInfo)    
    
   
    allSeq[seqs[nameKeyIDStr]] = seqs
    allSeqData[seqs[nameKeyIDStr]] = seqs_data

    return   allSeqData



def readScoreTab_file(inputfilename,sep="\t"):
    '''
    read score tab file
    each line contains one sequence 
    sequence name sep. by tab from the score
    the score are sep. by ;
    '''

    inputFile = open(inputfilename, "r")
    lines = inputFile.readlines()
    
    
    allSeqData = {}
    seqName = ""
    for line in lines:
        lineTokens = line.split(sep)
        # first token is the seq name 
        seqName = lineTokens[0].rstrip()
        # second token is the score separated by ; 
        parScores = lineTokens[1].rstrip()
        # if last char is ; remove it
        if parScores[len(parScores)-1] == ';':
            parScores = parScores[0:(len(parScores)-1)]
        scoreVector = np.array(parScores.split(';'),dtype=float)
        #allSeqData[seqName] = pd.DataFrame(scoreVector,columns=['score'])
        allSeqData[seqName] = scoreVector
    return  allSeqData


def writeScoreTab_file(finalScores,outputFilename):
    '''
    read score tab file
    each line contains one sequence 
    sequence name sep. by tab from the score
    the score are sep. by ;
    '''

    outputFileStream = open(outputFilename, "w+")
    for seqName in finalScores:
        outputFileStream.write("{0}\t".format(seqName))
        for i in range(0,len(finalScores[seqName])):
            outputFileStream.write("{0};".format(finalScores[seqName][i]))
        outputFileStream.write("\n")
    outputFileStream.close()
#%%
def consShifted(sample_in, f_list = ['C','INDEX']  ,  n = 3):


    sample = sample_in [ f_list  ] 
    # n is odd number
    m = int(n/2)
    
    l = len(sample_in)
    #print(m,l)
    pivotS = sample.iloc[m:(l-n+1+m),:]
    pivotS.index = range(0,pivotS.shape[0]) 
    # pivotS.columns = pivotS.columns + "_" + str(m)

    sample = sample_in [f_list]

    for i in range(0,m):
        subSample = sample.iloc[i:(l-n+i+1),:]
        subSample.columns = subSample.columns + "_" + str(i)
        #print(subSample.shape)
        subSample.index = range(0,subSample.shape[0]) 
        pivotS = pd.concat([pivotS,subSample], axis=1)
    for i in range(m+1,n):
        subSample = sample.iloc[i:(l-n+i+1),:]
        #print(subSample.shape)
        subSample.columns = subSample.columns + "_" + str(i)
        subSample.index = range(0,subSample.shape[0]) 
        pivotS = pd.concat([pivotS,subSample], axis=1)
        
        
    #cList = ['Class'+"_" + str(m)]
    cList = []
    for feature in f_list:
        for i in range(0,n):
            if i == m :
                cList.append(feature)
                continue
            cList.append(feature + '_'+str(i))
    
    

    
    return pivotS[ cList]
def getFeatureList(wSize):
    m = int(wSize/2)
    #fList = []
    fList = []
    for i in range(0,wSize):
        if i == m :
            fList.append('C')
            continue
        fList.append('C_'+str(i))
    return fList
def process(seqVector):
    # INDEX and char only
    seq_df = pd.DataFrame(seqVector,columns=['char'])
    seq_df['INDEX']= seq_df.index
    colName = "char"
    seq_df.loc[seq_df[colName] == "u" ,colName] = "U"
    seq_df.loc[seq_df[colName] == "g" ,colName] = "G"
    seq_df.loc[seq_df[colName] == "c" ,colName] = "C" 
    seq_df.loc[seq_df[colName] == "a"  ,colName] = "A" 
    seq_df.loc[seq_df[colName] == "t" ,colName] = "U" 
    seq_df.loc[seq_df[colName] == "T" ,colName] = "U"
    seq_df['C'] = 0
    seq_df.loc[seq_df[colName] == "U" ,'C']  = 1
    seq_df.loc[seq_df[colName] == "G" ,'C']  = 2
    seq_df.loc[seq_df[colName] == "C" ,'C']  = 3 
    seq_df.loc[seq_df[colName] == "A"  ,'C'] = 4 
    seq_df.loc[seq_df[colName] == "T" ,'C']  = 1
    return seq_df
def calc_p_RNN(sample,wSize,model):
    fList =  getFeatureList(wSize)
    sample_shift = consShifted(sample ,f_list=['C','INDEX'] ,  n = wSize)
    x_data = sample_shift[fList].values
    #y_data = sample['Class'].values
    res_p = nt_RNN_pred(model,wSize,x_data)
    sample_shift["RNN_"+str(wSize)] = res_p
    sampleF = sample_shift[["INDEX","RNN_"+str(wSize)]]
    sample = pd.merge(sample, sampleF, left_on=['INDEX'],right_on=['INDEX'],how='left')
    sample.fillna(value=0.5,inplace=True)
    sample.index = sample.index
    return sample
def nt_RNN_pred(nt_model,n,x_data):
    x_data = toDense(x_data,n)
    x_data = np.reshape(x_data, (len(x_data), n, 4))
    res = nt_model.predict(x_data,batch_size=16)
    return res

# convert matrix rows into one hot encoding
def toDense(t,n):
    dense = np.zeros((len(t),n*4))
    for di in range(0,len(t)):
        for i in range(0,n):
            nV = t[di,i]
            dense[di,((i)*4)+nV -1] = 1
    return dense


def calcFinalScore(rnn_score,parScore,w1=0.5,w2=0.5):
    # TODO :: make sure that both array hava same dim...
    return w1*rnn_score+w2*parScore

#%%
## calculate score for different k classifiers
def calcRNNScore(q_sample):
    ClsName = 'RNN'
    for wSize in range(7,16,2):
        model = RNN_models[wSize]
        sample = calc_p_RNN(q_sample,wSize,model)
        q_sample["RNN_"+str(wSize)] =  sample["RNN_"+str(wSize)]
    q_sample[ClsName] =  q_sample[ClsName+"_7"]
    for wSize in range(9,16,2):
        q_sample[ClsName] =  q_sample[ClsName] + q_sample[ClsName+"_"+str(wSize)]
    q_sample[ClsName]/=5
    q_sample[ClsName] = q_sample[ClsName]*2 - 1
    return q_sample


#%% Load models

databaseName = 'pdb'
modelFiles_Loc =  scriptpath + "/RNN_LSTM_" + databaseName + "_"

sys.stdout.write("Loading RNN Models :\n")
RNN_models = {}
for wSize in range(7,16,2):
    sys.stdout.write("\t\tLoading {0}-RNN model ....\n".format(wSize))
    model = keras.models.load_model(modelFiles_Loc + str(wSize) )
    RNN_models[wSize] = model

#%%

sys.stdout.write("Reading score tab file .... \n")
seqScore = readScoreTab_file(scoresTabFileName)


sys.stdout.write("Reading Fasta file .... \n")


rna_seqs = readFasta_file(fastaFileName)
finalScores = {}
for seqName in rna_seqs:
    sys.stdout.write(seqName+"\n")
    ## construct data matrix
    sys.stdout.write("Calculating RNN score for {0} sequence\n".format(seqName))
    q_sample = process(rna_seqs[seqName]['data'])
    q_sample = calcRNNScore(q_sample)
    if seqName in  seqScore:
        finalScore = calcFinalScore(q_sample['RNN'].values,seqScore[seqName],args.w1,args.w2)
    else:
        sys.stdout.write("\t\t No Pars Score provided for {0}. The final score will be only RNN.\n".format(seqName))
        finalScore = q_sample['RNN'].values
    finalScores[seqName] = finalScore

sys.stdout.write('Writing final scores file {0}\n'.format(outputScoreTabFileName))
writeScoreTab_file(finalScores,outputScoreTabFileName)
