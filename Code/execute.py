import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
#Loading the Coding RNA dataset 
with open('/path_to/coding.fasta') as fasta_file:
    identifiers=[]
    lengths=[]
    name=[]
    for seq_record in SeqIO.parse(fasta_file,'fasta'):
        identifiers.append(str(seq_record.seq))
        lengths.append(len(seq_record.seq))
        name.append(seq_record.id)
        
d={'sequence':identifiers,'Len':lengths,'Name':name}
coding=pd.DataFrame(d)

coding['class']=coding.apply(lambda x:1, axis=1)
coding_df=coding[coding['Len'].between(30,60)]
coding_df.drop('Name',axis='columns',inplace=True)
#Loading the non-coding RNA dataset
with open('/path_to/non-coding.fasta') as fasta_file:
    identifiers=[]
    lengths=[]
    name=[]
    for seq_record in SeqIO.parse(fasta_file,'fasta'):
        identifiers.append(str(seq_record.seq))
        lengths.append(len(seq_record.seq))
        name.append(seq_record.id)
 
 d={'sequence':identifiers,'Len':lengths,'Name':name}
nocoding=pd.DataFrame(d)
nocoding['class']=coding.apply(lambda x:0, axis=1) 
nocoding_df=nocoding[nocoding['Len'].between(28, 60)] 
nocoding_df.drop('Name',axis='columns',inplace=True)

#Combining the coding and non-coding RNA data into a single dataframe
import re
df_All = [coding_df, nocoding_df]
df_sequence = pd.concat(df_All)
p = re.compile(r'[^ACGT]')
df_sequence['sequence'] = df_sequence['sequence'].map(lambda x: x.replace('N', ''))

# This FCGR is a variant version of the R package kaos by Dominic Eger and Hannah Franziska Löchel, implemented in Python by Sandra Clemens
from math import sin, cos, pi, floor, ceil
BASES = {
"digits": [0,1,2,3,4,5,6,7,8,9],
"AMINO": ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"],
"amino": ["a","c","d","e","f","g","h","i","k","l","m","n","p","q","r","s","t","v","w","y"],
"DNA": ["A","G","T","C"],
"dna": ["a","g","t","c"]
}

class CGR:
    r = 1
    
    #preciion 8 = f8 = 64 bit float
    def __init__(self, data, seq_base=None, sf=None, precision='8'):
        self.data = data
        self.dtype = f"f{precision}"
        
        if seq_base is not None:
            if type(seq_base) is list:
                self.seq_base = seq_base
            else: 
                self.seq_base = BASES[seq_base]
        else:
            self.seq_base = list(sorted(set(data)))

        base_num = len(self.seq_base)
       
        if sf is not None:
            #todo:sf <= 1 sf >=1
            self.sf = sf
        else:
            self.sf =  1 - (sin(pi * (1/base_num)) / (sin(pi * (1/base_num)) + sin(pi * (1/base_num + 2*(floor(base_num/4) / base_num)))))
            
        if base_num == 4:
            corners = {
                "x":np.array([1,-1,-1,1], dtype=self.dtype),
                "y":np.array([-1,-1,1,1], dtype=self.dtype)
            }
        else:
            corners = self.__calc_corners()
        
        self.corners = pd.DataFrame.from_dict(corners, orient="index", columns=self.seq_base)
        
        self.coords = self.__calc_coords()
        
    def calc_fcgr(self, res=100):
        A = np.zeros((res,res), dtype=int)
        for i in range(len(self.data)):
            x = ceil((self.coords["x"][i]+self.r) * res/(2*self.r)) -1
            y = ceil((self.coords["y"][i]+self.r) * res/(2*self.r)) -1
            A[x,y] += 1
            
        return np.rot90(A)
    
    def __calc_coords(self):
        coords = {
            "x" : np.zeros(len(self.data), dtype=self.dtype),
            "y" : np.zeros(len(self.data), dtype=self.dtype)
        }
        last_x = 0
        last_y = 0
        for i, base in enumerate(self.data):
            corner = self.corners[base]
            last_x = last_x + (corner.x - last_x) * self.sf
            last_y = last_y + (corner.y - last_y) * self.sf
            coords["x"][i] = last_x
            coords["y"][i] = last_y
         
        return coords 
         
    def __calc_corners(self):
        base_num = len(self.seq_base)
        corners = {
            "x" : np.zeros(base_num, dtype=self.dtype),
            "y" : np.zeros(base_num, dtype=self.dtype)
        }
        for i in range(base_num):
            tmp = 2*i+1
            corners["x"][i] = self.r*sin(pi * (tmp/base_num))
            corners["y"][i] = self.r*cos(pi * (tmp/base_num))
        
        return corners
        
eval_sequence = df_sequence.iloc[:, 0].values
J = []
for x in eval_sequence:
    cgr = CGR(data=x, seq_base=["A","G","T","C"])
    fcgr_seq = cgr.calc_fcgr(res=8)
    seq = fcgr_seq
    J.append(seq)
    eval_seq = np.array(J)
    
eval_labels = df_sequence.iloc[:, 2].values
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(eval_seq, eval_labels,test_size = 0.25, random_state=45)

from tensorflow.keras.models import load_model
cnnModel = load_model ('RNAIntel.h5')
cnnModel.summary()

y_pred_proba = cnnModel.predict(X_test)
fpr, tpr, thresholds = roc_curve(y_test, y_pred_proba)
y_pred = (y_pred_proba > 0.5).astype(int)

from sklearn.datasets import make_circles
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import matthews_corrcoef

accuracy = accuracy_score(y_test, y_pred)
print('Accuracy: %f' % accuracy)
precision = precision_score(y_test, y_pred)
print('Precision: %f' % precision)
f1 = f1_score(y_test, y_pred)
print('F1 score: %f' % f1)
auc = roc_auc_score(y_test, y_pred)
print('ROC AUC: %f' % auc)
matrix = confusion_matrix(y_test, y_pred)
print(matrix)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
print(tn, fp, fn, tp)
specificity = tn / (tn+fp)
print('Specificity: %f' % specificity)
sensitivity = tp/ (tp+fn)
print('Sensitivity: %f' % sensitivity)
matthews_corrcoef =  matthews_corrcoef(y_test, y_pred)
print(' MCC: %f' % matthews_corrcoef)




   
