# coding=utf-8

  

import pickle as cPickle #import cPickle
import os
import string
import time
import math
import numpy as np
import subprocess 
import shutil
 
 
    
def getPrefixLists(inDir):   
    if os.path.isdir( inDir)==False  :    #directory is not existed
        return None
              #create directory
    fileNames=os.listdir(inDir)
    res=[]
    for  file in fileNames:
        res.append(get_prefix(file))    
    return res
def get_prefix(file_path):
    try:
        lists=file_path.split(os.sep)
        names=lists[-1].split(".fasta")
        return names[0]
    except:
        return None

 
def hhblits_schedule_(inDir,outDir,
                      dbName='/pub/db/uniref30_hhsuite_db'):
                 
    '''
    When more than one fasta files are input.
    '''
    if os.path.isdir( inDir)==False :
        print("input file directory error")
        return None
    if os.path.exists(outDir)==False:    #mkdir output directory
        os.mkdir(outDir)                                                    
    fileNames=os.listdir(inDir)
                                                     
    for file in fileNames :
        outPrefix=get_prefix(file)            #get file prefix                               
        in_file=inDir+os.sep+file                #input fasta file
        outFile=outDir+os.sep+outPrefix+".hhm"    #output pssm file
        dbName=' -d ' +dbName
         
        hhblits_schedule(in_file,out_file,dbName)
    return 

                                          
def hhblits_schedule(fasta_file,out_file,
                     dbName=' -d /pub/uniref30_hhsuite_db/UniRef30_2020_06'):
                  
    txtOutFile=out_file+"1"    
    command="hhblits -i "+fasta_file+dbName+ " -ohhm " +out_file+' -cpu 4 -maxres 40000'
    try:
        #subprocess.check_call(command,shell=True)
        os.system(command)
        
    except:
        print( "execute hhblits error!!!")
        return false
    return True

def gen_temp_file(inputFasta):
    '''
    input: fasta file.
    output: fasta file,pssm file, hmm file, pssm text file, temp work directory.
    '''
    tmp_dir="../work_tmp"
    print(tmp_dir)
    if os.path.exists(tmp_dir)==False:
        os.mkdir(tmp_dir) 
    outPrefix=get_prefix(inputFasta)            #get file prefix
                  #input fasta file
    
    pssmFile=tmp_dir+os.sep+outPrefix+".pssm"    #output pssm file
    txtPssmFile=tmp_dir+os.sep+outPrefix+".txt"
    hhmFile=tmp_dir+os.sep+outPrefix+".hhm"    #output hhm file
    fastaFile=tmp_dir+os.sep+outPrefix+".fasta"   
    
    shutil.copyfile(inputFasta, fastaFile)

    return fastaFile,pssmFile,txtPssmFile,hhmFile,tmp_dir
        
    
        
def blast_schedule_(inDir,outDir,tmp_dir,
                    dbName='/pub/unirefdb/uniref90'):
                    
    '''
    When more than one fasta files are input.
    '''
     
     
    if os.path.isdir( inDir)==False :
        print("input file directory error")
        return False
    if os.path.exists(outDir)==False:    #mkdir output directory
        os.mkdir(outDir)    
    if os.path.exists(tmp_dir)==False:    #mkdir temp directory
        os.mkdir(tmp_dir)    
    fileNames=os.listdir(inDir)
    for file in fileNames :
        outPrefix=get_prefix(file)            #get file prefix
        inFile=inDir+os.sep+file                #input fasta file
        outFile=outDir+os.sep+outPrefix+".pssm"    #output pssm file
        txtOutFile=tmp_dir+os.sep+outPrefix+".pssm"+"1"
        blast_schedule(inFile,outFile,txtOutFile,dbName)
      

def blast_schedule(inFile,outFile,txtOutFile,
                   dbName='/pub/unirefdb/uniref90'):
                  

    command="psiblast -query "+inFile+" -db "+dbName+ \
            " -num_threads 4   -out "+txtOutFile+"-save_pssm_after_last_round  -out_ascii_pssm "+ +outFile+" -num_iterations 3"
    try:
        p=subprocess.Popen(command,shell=True)
        p.wait()
    except:
        print( "execute psiblast error!!!")
        return False
    return True 


  
        
def main():
     
    file='test_fasta.txt'   
    gen_fasta_files(file,'web_site')
    
   
if __name__=="__main__": 
    main()
