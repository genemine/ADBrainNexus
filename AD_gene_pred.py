#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 17:21:00 2021

@author: lincxcsu
"""


import numpy as np
import pandas as pd
import random,sys,os
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score,roc_curve,average_precision_score,precision_recall_curve
from sklearn.linear_model import RidgeCV

#### input and output settings
script,out,featurematrix,totaltimes=sys.argv

model=RidgeCV(alphas=[0.01, 0.1, 1, 10, 20, 30])
os.system('mkdir '+out+'/')
############# below DO NOT change  ##################
def to_list(filex):
    f=open(filex)
    result=[]
    for line in f:
        t=line.strip()
        if len(t)>0:
            result.append(t)
    return result


### 5-fold CV
print('5-fold CV')
def kfoldcv(data,CV_file,gene,negs):
    random.seed(0)
    score={}
    count={}
    AUROC=[]
    AUPRC=[]
    for i in range(totaltimes):
        print('step:'+str(i))
        nposigenes= len(gene) 
        negative=random.sample(list(negs),len(negs))
        index=gene+negative
        y=[1]*nposigenes+[0]*len(negative)
        y=pd.Series(y)
        y.index=index
        X=data.loc[index,:]
        skf = StratifiedKFold(n_splits=5)
        y_pred=pd.Series([],dtype='float')
        y_new=pd.Series([],dtype='float')
        for train,test in skf.split(X,y):
            X_train, X_test = X.iloc[train,:], X.iloc[test,:]
            y_train, y_test = y.iloc[train], y.iloc[test]
            model.fit(X_train,y_train)
            pred=model.predict(X_test)
            y_predi=pd.Series(pred)
            y_predi.index=y_test.index
            y_new=y_new.append(y_test)
            y_pred=y_pred.append(y_predi)
            left=data.index.difference(X_train.index)
            X_left=data.loc[left,:]
            Pred_left=model.predict(X_left)
            Pred_left=pd.Series(Pred_left)
            Pred_left.index=left
            for genei in Pred_left.index:
                if genei not in score:
                    score[genei]=Pred_left[genei]
                    count[genei]=1
                else:
                    score[genei]=score[genei]+Pred_left[genei]
                    count[genei]=count[genei]+1
        AUROCi=round(roc_auc_score(y_new, y_pred),4)
        AUPRCi=round(average_precision_score(y_new, y_pred),4)
        print('AUROC:'+str(AUROCi))
        print('AUPRC:'+str(AUPRCi))
        AUROC.append(AUROCi)
        AUPRC.append(AUPRCi)
    fnew=open(CV_file,'w')
    result={}
    result['AUPRC']=AUPRC
    result['AUROC']=AUROC
    fnew.write(str(result))
    fnew.close()
    result={}
    for genei in data.index:
        if genei in score:
            result[genei]=score[genei]/count[genei]
    x=pd.Series(result)
    x=x.sort_values(ascending=False)
    x.to_csv(out+'/prediction.txt',sep='\t')


###main####
        
totaltimes=int(totaltimes)
posigenes=to_list('data/positive.txt')
negs=to_list('data/negative.txt')
CV_file=out+'/5folds_CV'+'.tsv'


#preprocessing
data=pd.read_csv(featurematrix,index_col=0,header=0,sep='\t')
data.index=[str(indexi) for indexi in data.index]
data=data.dropna()
data = data[~data.index.duplicated(keep = "first")]
gene=list(set(data.index.intersection(data.columns.intersection(posigenes))))
negs=list(set(data.index.intersection(negs)).difference(gene))


#model performance
kfoldcv(data,CV_file,gene,negs)
# end of 5-fold CV
