
#coding=utf-8



import pickle as cPickle #import cPickle
import os
import string
import numpy as np
import tensorflow as tf
from keras import backend as K
from keras.models import save_model,load_model
     


 
 
def get_model(saved_model_file):
    
    if os.path.exists( saved_model_file)==False  :   
        return None     
    model=load_model(saved_model_file)
    
    return model

def convert_to_angle (seq_len,outs_sin,outs_cos):
    outs =[]
    for i in range(seq_len):
        predict_angle=np.rad2deg(np.arctan2(outs_sin[0][i][0], outs_cos[0][i][0]))  
        outs.append(round(predict_angle,1))  
    return outs


            

 
def   write_to_file(pred_output,sequence,output_file):
    '''
    write prediction to output file.
    '''   
    with open(output_file,'w') as f:
        f.write('    #\tAA\t'  +'PHI\tPSI\t'+'\n')
        for j,residue in enumerate(sequence):
            f.write('%5d\t%5s'%(j,residue)+'\t%3.1f\t%3.1f\t\n'%(pred_output[j][:]))





def  run_model(feat22,feat20,feat30,sequence,
				model_file='../model/dcma_model.h5'):
            
    time_step=700
    batch_size=100
 
    model = load_model(model_file, custom_objects={'tf': tf})
    out_phi_sin,out_phi_cos,out_psi_sin,out_psi_cos \
                =model.predict([feat20,feat22,feat30],batch_size=batch_size)
    pred_phi=convert_to_angle (len(sequence),out_phi_sin,out_phi_cos)
    pred_psi=convert_to_angle (len(sequence),out_psi_sin,out_psi_cos)        
    return list(zip(pred_phi,pred_psi))
    

   
