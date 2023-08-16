#coding=utf-8



import math
import os
import numpy as np
import re
ISOTIMEFORMAT='%Y-%m-%d %X'
 

 
# our orthogonal encoding method...

orth_protein={"A":0,"C":1,"D":2,"E":3,"F":4,"G":5,"H":6,"I":7,"K":8,\
          "L":9,"M":10,"N":11,"P":12,"Q":13,"R":14,"S":15,"T":16,"V":17,"W":18,"Y":19, \
            'NoSeq':20,'X':21,'!':21,'U':21}
def load_one_2_vec(one_2_vec_file=
        "data/one-vec.npy"):
    one_hot_2_vec=[]
    data_set = np.load(one_2_vec_file ,allow_pickle=True) 
    one_hot_2_vec.extend(data_set )  
    return one_hot_2_vec
    

def prepare_input(fastaFile,pssmFile,hhmFile):
     
    one_hot_2_vec=load_one_2_vec()
    dim1=1;time_step=700
    input_dim_a=20; input_dim_b=22;input_dim_c=30

    hhmProfiles=read_hmm(hhmFile)
    if hhmProfiles==None:
        print ('no .....',hhmFile)
        return None, None, None, None
    pssmLines=extract_pssm(pssmFile)  
    if pssmLines==None:
        print ('no .....',pssmFile)
        return None, None, None, None        
    feat20=np.zeros(shape=( dim1, time_step,input_dim_a ), dtype=np.float32)
    feat22=np.zeros(shape=( dim1, time_step,input_dim_b ), dtype=np.float32)
    feat30=np.zeros(shape=( dim1, time_step,input_dim_c ), dtype=np.float32)

    totalResidues=len(pssmLines)
 
     
    i=0
    sequence=''
    while i<totalResidues:

        pssmLine = pssmLines[i]


        hhmResidue, hhmVector=hhmProfiles[0][i],hhmProfiles[1][i]
        pssmResidue  =   pssmLine[-1]
        if pssmResidue!=hhmResidue:
            return None, None, None, None
        [tmp_hhm]=np.asarray(hhmVector)
        feat22[0][i][0:22]  =one_hot_2_vec[0][ orth_protein[ pssmResidue] ]
        feat20[0][i][0:20]  = pssmLine[0:20]
        feat30[0][i][0:30] =tmp_hhm[:]

        sequence+=pssmResidue
        i+=1

    return feat22,feat20,feat30,sequence
 
    #1st for end
# end function    


 

 

  

def extract_lines(pssmFile): 
    try:  
        fr2 =   open(pssmFile)   # pssm file open
        pssmLines=[]            #read pssp lines to memory
        if  fr2==None:    #files null
            return 
        #read pssm data ........................................
        for i in range(3):
            fr2.readline()      #skip the 3rd line. 
        while True:
            psspLine=fr2.readline()
            if  psspLine.strip()=='' or psspLine.strip()==None:
                break   #end reading.
            pssmLines.append(psspLine)  #append line.
        fr2.close()
        return   pssmLines  
    except:
        return None


def extract_pssm(pssm_file):
    '''
    extract one sequence PSSM matrix from file
    input: filename
    output:2-dim matrix [residues][20dim PSSM + one residue]
    '''
    pssmLines = extract_lines(pssm_file)  # read pssp lines to memory
    if pssmLines == None:
        return None
    i = 0;totalResidues = len(pssmLines)
    featLines = []
    while i < totalResidues:
        tmp = []
        pssmLine = pssmLines[i]
        pssmResidue = pssmLine[5:7].strip()
        nums = [int(x) for x in re.findall(r"-\d+\.?\d*|\d+\.?\d*", pssmLine[7:])[:20]]
        for count in range(20):


            item = 1.0 / (1 + math.exp(-1 * nums[count]))  # normalize to (0,1)

            tmp.append(item)
        tmp.append(pssmResidue)

        featLines.append(tmp)
        i += 1
        # break
    print(pssm_file, len(featLines))
    return featLines


def extract_pssm2(pssm_file): #
    '''
    extract one sequence PSSM matrix from file
    input: filename
    output:2-dim matrix [residues][20dim PSSM + one residue]
    '''
    pssmLines=extract_lines(pssm_file)            #read pssp lines to memory
    if pssmLines==None:
        return None
    i=0; totalResidues=len(pssmLines)
    featLines=[]
    while i<totalResidues:
        tmp=[]
        pssmLine    =   pssmLines[i]
        pssmLine = pssmLine.split(" ")

        pssmLine = [i for i in pssmLine if i != '']
        pssmResidue =  pssmLine[0:2]
        pssmLine = pssmLine[2:22]

        for count in range(20):
                # pssm matrix value like:-9-10-10 ............20170904
                #tmp.append(string.atoi(spLine[count]))
            # item=int(pssmLine[9+count*3:9+(count+1)*3])
            #item=int( pssmLine[9+count*4:9+(count+1)*4] )
            # nums = [int(x) for x in re.findall(r"-\d+\.?\d*|\d+\.?\d*", pssmLine[7:])[:20]]
            item= 1.0/(1+math.exp(-1* float(pssmLine[count]) ))  # normalize to (0,1)

            tmp.append(item )
        tmp.append(pssmResidue)

        featLines.append(tmp)
        i+=1
        #break
    print (pssm_file,len(featLines) )


    return featLines

def read_hmm(hhm_file):
    if not os.path.exists(hhm_file):
        print ('----------------',hhm_file)
        return None #files null
    f = open(hhm_file)
    if f==None :    #files null
        return None
    
    line = f.readline()
    while line[0] != '#':
        line = f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    seq = []
    extras = np.zeros([0, 10])
    prob = np.zeros([0, 20])
    line = f.readline()
    while line[0:2] != '//':
        lineinfo = line.split()
        seq.append(lineinfo[0])
        probs_ = [2 ** (-float(lineinfo[i]) / 1000) if lineinfo[i] != '*' else 0. for i in range(2, 22)]
        prob = np.concatenate((prob, np.matrix(probs_)), axis=0)
        line = f.readline()
        lineinfo = line.split()
        extras_ = [2 ** (-float(lineinfo[i]) / 1000) if lineinfo[i] != '*' else 0. for i in range(0, 10)]
        extras = np.concatenate((extras, np.matrix(extras_)), axis=0)
        line = f.readline()
        assert len(line.strip()) == 0
        line = f.readline()
    # return (''.join(seq),prob,extras)
    return (seq, np.concatenate((prob, extras), axis=1))    


 
 
