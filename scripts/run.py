# coding=utf-8


import pickle as cPickle #import cPickle
import os
import string
import time
import math
import numpy as np
import argparse
import shutil
import tensorflow as tf
from gen_profile import blast_schedule,hhblits_schedule,gen_temp_file
from prepare_for_input import prepare_input
from exec_model import run_model, write_to_file

 
        
    

def execute(inputFasta,outputFile):
    
        
    fastaFile,pssmFile,txtPssmFile,hhmFile,tmp_dir=gen_temp_file(inputFasta)
    '''
    if not blast_schedule(fastaFile,pssmFile,txtPssmFile) :
        return False
    '''
    if not hhblits_schedule(fastaFile,hhmFile):
        return False

    feat22,feat20,feat30,sequence=prepare_input(fastaFile,pssmFile,hhmFile)
    
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' # shut debug info... 
    pred_output=run_model(feat22,feat20,feat30,sequence)
    write_to_file(pred_output,sequence,outputFile)
    

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputFasta', help='input sequence fasta file')
    parser.add_argument('-o', '--outputFile', help='output prediction file')
    args = parser.parse_args()
    inputFasta =args.inputFasta    # input
    outputFile=args.outputFile  #output
    print ("outputFile,",outputFile)
    execute(inputFasta,outputFile)



if __name__ == "__main__":
    main()
    
    
