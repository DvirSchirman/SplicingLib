# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:07:01 2019

@author: dvirs
"""

import lightgbm as lgb
import pandas as pd
import numpy as np
import scipy as sp
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import shap
from random import sample
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold
from sklearn import metrics

df_classes = pd.read_csv('alternative_isoforms_multiclass.csv',index_col=0)
df_features = pd.read_csv('alternative_sp_model_features_all.csv',index_col=0)
X_features = df_features.drop(['splicing_eff_median'], axis=1)
X_features = np.concatenate((X_features.iloc[0::3,:].values,X_features.iloc[1::3,:].values,X_features.iloc[2::3,:].values),axis=1)

y=df_classes['class']

seed = 20191024
test_size = 0.25


X=np.concatenate((X_features[:,:76],X_features[:,79:81]),axis=1)
                
X_train = X
y_train=y

params = {
    'boosting_type': 'gbdt',
    'objective': 'multiclass',
    'num_class':8,
    'metric': 'multi_logloss',
    'num_leaves': 50,
    'learning_rate': 0.1,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'importance_type': 'gain',
    'categorical': [0,1,2,38,39,40],
    'verbose': 0
    
}

kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)

shap_values = None

i=0
class_pred_mat=np.empty((X.shape[0]))
for train_index, val_index in kf.split(X_train, y_train):
    X_train_i, X_val = X_train[train_index,:], X_train[val_index,:]
    y_train_i, y_val = y[train_index], y[val_index]
    
    lgb_train = lgb.Dataset(X_train_i, y_train_i)
    lgb_eval = lgb.Dataset(X_val, y_val, reference=lgb_train)
    

    gbm = lgb.train(params,
                    lgb_train,
                    num_boost_round=500,
                    valid_sets=lgb_eval,
                    early_stopping_rounds=5,
                    )
    
    y_pred=gbm.predict(X_val)
    class_pred_mat[val_index]=np.argmax(y_pred,axis=1)
 
    explainer = shap.TreeExplainer(gbm)
    tmp = np.array(explainer.shap_values(X_train))
    
    if shap_values is None:
        shap_values = tmp
    else:
        shap_values += tmp
        
    i+=1

shap_values /= 5

suffixes = ['_1','_2']
feature_names=[]
for suffix in suffixes:
    for name in df_features.columns[1:]:
        name_tmp = name + suffix
        feature_names.append(name_tmp)

feature_names.append('intron_GC_J')
feature_names.append('intron_len_J')

if X_train.shape[1]>79:
    feature_names.append('gbm_1')
    feature_names.append('gbm_2')
    feature_names.append('gbm_J')


##
cat_features=[0,1,2,38,39,40]
for feature in cat_features:
    X_train[:,feature]=X_train[:,feature].astype('float')
    

df_ss5_map=pd.read_csv('ss5_map.csv')
for f in df_ss5_map.values:
    X_train[np.where(X_train[:,0]==f[1]),0]=f[2]
    X_train[np.where(X_train[:,38]==f[1]),38]=f[2]
    
df_branch_map=pd.read_csv('branch_map.csv')
for f in df_branch_map.values:
    X_train[np.where(X_train[:,1]==f[1]),1]=f[2]
    X_train[np.where(X_train[:,39]==f[1]),39]=f[2]
    
df_ss3_map=pd.read_csv('ss3_map.csv')
for f in df_ss3_map.values:
    X_train[np.where(X_train[:,2]==f[1]),2]=f[2]
    X_train[np.where(X_train[:,40]==f[1]),40]=f[2]   
    
##
    
for i in range(8):
    shap.summary_plot(shap_values[i,:,:], X_train,max_display=4,feature_names=feature_names,show=False)
    filename = 'feature_importance_multiclass_' + str(i) + '.png'
    # plt.savefig(filename)
    plt.show()
    filename = 'alternative_isoform_shap_values_multiclass_' + str(i) + '.csv'
    np.savetxt(filename,shap_values[i,:,:])
    
inds = [0,1,5,6,39,43]
n=0
X_inds=np.empty((X_train.shape[0],len(inds)))
shap_inds=np.empty((8,X_train.shape[0],len(inds)))
features_inds=[]
for i in inds:
    X_inds[:,n]=X_train[:,i]
    shap_inds[:,:,n]=shap_values[:,:,i]
    features_inds.append(feature_names[i])
    n+=1
for i in range(8):
    shap.summary_plot(shap_inds[i,:,:], X_inds,feature_names=features_inds,sort=False)
    plt.show()
    filename = 'shap_values_isoforms/isoform_shap_values_multiclass_' + str(i) + '.csv'
    np.savetxt(filename,shap_inds[i,:,:])
    filename = 'shap_values_isoforms/isoform_feature_values_multiclass_' + str(i) + '.csv'
    np.savetxt(filename,X_inds)
    
print(metrics.confusion_matrix(y, class_pred_mat))
print(metrics.classification_report(y, class_pred_mat, digits=3))
    

#shap.summary_plot(shapJ_values[:,-3:], X_train[:,-3:],max_display=8,feature_names=feature_names[-3:])

#np.savetxt('alternative_isoform_gbm_combined_predictions_binary.csv',class_pred_mat)

#state=1
#shap.summary_plot(shap1_values, X_train,max_display=8,feature_names=feature_names)

#shap.summary_plot(shap1_values[:,37:74], X_train[:,37:74],max_display=8,feature_names=feature_names[37:74])
#y1_pred = gbm1.predict(X_test, num_iteration=gbm1.best_iteration)
#class1_pred=y1_pred
#class1_pred[class1_pred<0.5]=0
#class1_pred[class1_pred>=0.5]=1
#
#y2_pred = gbm2.predict(X_test, num_iteration=gbm2.best_iteration)
#class2_pred=y2_pred
#class2_pred[class2_pred<0.5]=0
#class2_pred[class2_pred>=0.5]=1
#
#yJ_pred = gbmJ.predict(X_test, num_iteration=gbmJ.best_iteration)
#classJ_pred=yJ_pred
#classJ_pred[classJ_pred<0.5]=0
#classJ_pred[classJ_pred>=0.5]=1
#
#y_test=yJ_test
#class_pred=classJ_pred
#tp=(y_test.values)*(class_pred)
#fp=(1-y_test.values)*(class_pred)
#tn=(1-y_test.values)*(1-class_pred)
#fn=(y_test.values)*(1-class_pred)
#precision=sum(tp)/(sum(tp)+sum(fp))
#recall=sum(tp)/(sum(tp)+sum(fn))
#f1=2*precision*recall/(precision+recall)
#print('[Precision, Recall, f1] = [',precision,recall,f1,']')
#
#explainer = shap.TreeExplainer(gbmJ)
#shap_values = explainer.shap_values(X_train)
#shap.summary_plot(np.array(shap_values)[1,:,:], X_train,max_display=8)
#
#
#
