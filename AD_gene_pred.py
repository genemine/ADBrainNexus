#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 11:45:09 2021

@author: lincxcsu
"""


import numpy as np
import pandas as pd
import random,sys
from sklearn.ensemble import  ExtraTreesClassifier,RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score,roc_curve,average_precision_score,precision_recall_curve
from sklearn import svm 
from sklearn.linear_model import LogisticRegression 

#### input and output settings
script,featurematrix,outputfileprefix,totaltimes=sys.argv

totaltimes=int(totaltimes)
allnegagenes ='./data/negative.txt'
posigenes='./data/gene.txt'
CV_file=outputfileprefix+'_5folds_CV'+'.tsv'

#preprocessing
data=pd.read_csv(featurematrix,index_col=0,header=0,sep='\t')
data.index=[str(indexi) for indexi in data.index]
data=data.abs()
columns=data.loc[(data.std(axis=1)<0.01),:].index.intersection(data.columns)
data=data.drop(columns,axis=1)
indexs=data.loc[(data.std(axis=1)>0.01),:].index
data=data.loc[indexs,:]
data=data.dropna()
data = data[~data.index.duplicated(keep = "first")]
models={}
models['ExtraTrees']=ExtraTreesClassifier(n_estimators = 500,class_weight = "balanced_subsample", n_jobs = -1, bootstrap = True) 
#models['RandomForest']=RandomForestClassifier(n_estimators = 500,class_weight = "balanced_subsample", n_jobs = -1, bootstrap = True) 
#models['SVM']=svm.SVC(kernel = "linear",probability = True,class_weight = "balanced",C=0.01,gamma=100)
#models['LogisticRegression']=LogisticRegression(solver='lbfgs',multi_class='multinomial',class_weight={0:0.5, 1:0.5})

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
def kfoldcv(data,CV_file,gene):
    AUROC={}
    AUPRC={}
    score={}
    for modeli in models:
        AUROC[modeli]=[]
        AUPRC[modeli]=[]
    for i in range(totaltimes):
        print('step:'+str(i))
        nposigenes= len(gene)
        negs=to_list(allnegagenes)
        negs=list(data.index.intersection(negs))
        negative=random.sample(list(negs),len(negs))
        index=gene+negative
        y=[1]*nposigenes+[0]*len(negative)
        y=pd.Series(y)
        y.index=index
        X=data.loc[index,:]
        skf = StratifiedKFold(n_splits=5)
        for modeli in models:
            y_pred=pd.Series([])
            y_new=pd.Series([])
            for train,test in skf.split(X,y):
                X_train, X_test = X.iloc[train,:], X.iloc[test,:]
                y_train, y_test = y.iloc[train], y.iloc[test]
                model=models[modeli]
                model.fit(X_train,y_train)
                pred=model.predict_proba(X_test)[:,1]
                y_predi=pd.Series(pred)
                y_predi.index=y_test.index
                y_new=y_new.append(y_test)
                y_pred=y_pred.append(y_predi)
                left=data.index.difference(X_train.index)
                X_left=data.loc[left,:]
                P_left=model.predict_proba(X_left)[:,1]
                P_left=pd.Series(P_left)
                P_left.index=left
                for genei in P_left.index:
                    if genei not in score:
                        score[genei]=[P_left[genei]]
                    else:
                        score[genei].append(P_left[genei])
            AUROC[modeli].append(round(roc_auc_score(y_new, y_pred),4))
            AUPRC[modeli].append(round(average_precision_score(y_new, y_pred),4))
            result=pd.DataFrame(y_pred)
            result['label']=y_new
            result.to_csv(outputfileprefix+'pred_aurc.txt',sep='\t')
            print(modeli+':'+str(AUROC[modeli]))
            print(modeli+':'+str(AUPRC[modeli]))
    fnew=open(CV_file,'w')
    result={}
    result['AUPRC']=AUPRC
    result['AUROC']=AUROC
    fnew.write(str(result))
    fnew.close()
    ave={}
    ave=pd.Series(ave)
    ave=ave.sort_values(ascending=False).round(4)
    ave.to_csv(outputfileprefix+'_pred.txt',sep='\t')


gene=list(set(data.index.intersection(to_list(posigenes))))
print(len(gene))
#model performance

kfoldcv(data,CV_file,gene)

# end of 5-fold CV






