Readme.md 

This repository provides the companying data sets and scripts for the paper:

Buzhong Zhang, Meili Zheng, Yuzhou Zhang, Lijun Quan. DCMA: a Faster Protein Backbone Dihedral Angles Prediction by Dilated Convolutional Attention-Based Neural Network.  

which presents DCMA for predicting protein backbone angles .

The running platform:

​    Linux, Python 3.7 , Tensorflow 1.13, Keras 2.2.5 and other utilities for preparing data. 

Download the package DCMA.zip and unzip it. The folder structure and its explanation:

**/scripts**: All runnable files. 

- exec_model.py:	execute model and format output data.
- gen_profile.py:	generate PSSM and HHM profile.
- prepare_for_input.py: obtain input data feeding to DCMA.  Features of PSSM  are normalized to (0, 1) by logistic function. Features of sequence coding(one-hot vector) are mapped to dense vector.
- run.py:	run DCMA model.
 
 

**/work_tmp**: temp directory for storing raw PSSM and HHM profile.

**/model**: The trained model file of DCMA.

**/data**:  
    one-vec.dat  Sequence coding vectors where are used to map a one-hot vector to  a Densed vector.

Thank you!
If you have any suggestions or questions, Please email to:
zhbzhong@aqnu.edu.cn
