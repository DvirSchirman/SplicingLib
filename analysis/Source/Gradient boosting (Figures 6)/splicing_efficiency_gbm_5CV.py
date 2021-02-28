# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 17:40:14 2019

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


df=pd.read_csv('model_features_all.csv',index_col=0)
y=df['splicing_eff_median']
X=df.drop(['splicing_eff_median'], axis=1)



seed = 20191007
test_size = 0.25
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=seed)


cat_features=['ss5_seq','branch_seq','ss3_seq']
for feature in cat_features:
    X_train[feature]=X_train[feature].astype('category')
    X_test[feature]=X_test[feature].astype('category')
    
shap.initjs()

params = {
    'boosting_type': 'gbdt',
    'objective': 'regression',
    'metric': {'l2', 'l1'},
    'num_leaves': 50,
    'learning_rate': 0.1,
    'feature_fraction': 0.8,
    'bagging_fraction': 0.8,
    'bagging_freq': 5,
    'importance_type': 'gain',
    'verbose': 0
}

gbm = lgb.LGBMRegressor(params,
                num_boost_round=500,
                )

kf = KFold(n_splits=5, shuffle=True, random_state=seed)
predicts = []


shap_values = None
for train_index, val_index in kf.split(X_train, y_train):
    print("###")
    X_train_i, X_val = X_train.iloc[train_index], X_train.iloc[val_index]
    y_train_i, y_val = y_train.iloc[train_index], y_train.iloc[val_index]
    for feature in cat_features:
        X_train_i[feature]=X_train_i[feature].astype('category')
        X_val[feature]=X_val[feature].astype('category')
        
    lgb_train = lgb.Dataset(X_train_i, y_train_i)
    #                        categorical_feature=['ss5_seq','branch_seq','ss3_seq','ss5_stem_dir','branch_stem_dir','ss3_stem_dir'])
    lgb_eval = lgb.Dataset(X_val, y_val, reference=lgb_train)
    gbm = lgb.train(params,
                lgb_train,
                num_boost_round=500,
                valid_sets=lgb_eval,
                early_stopping_rounds=5,
                )
    predicts.append(gbm.predict(X_test))
    
    explainer = shap.TreeExplainer(gbm)
    if shap_values is None:
        shap_values = explainer.shap_values(X_train)
    else:
        shap_values += explainer.shap_values(X_train)      
    
y_pred=np.array(predicts).mean(axis=0)

shap_values /= 5

print('The rmse of prediction is:', mean_squared_error(y_test, y_pred) ** 0.5)

plt.scatter(y_pred,y_test,alpha=0.5)
plt.ylim([0,1])
plt.xlim([0,1])
fig=plt.gcf()
plt.show()
fig.savefig('../../Figures/Figure5/A - corr.png')
r=np.corrcoef(y_pred,y_test)
r=r[0,1]
print('Pearson corr:', r)


##
for feature in cat_features:
    X_train[feature]=X_train[feature].astype('float')
    X_test[feature]=X_test[feature].astype('float')

df_ss5_map=pd.read_csv('ss5_map.csv')
for f in df_ss5_map.values:
    X_train['ss5_seq'].iloc[np.where(X_train['ss5_seq']==f[1])]=f[2]
    
df_branch_map=pd.read_csv('branch_map.csv')
for f in df_branch_map.values:
    X_train['branch_seq'].iloc[np.where(X_train['branch_seq']==f[1])]=f[2]
    
df_ss3_map=pd.read_csv('ss3_map.csv')
for f in df_ss3_map.values:
    X_train['ss3_seq'].iloc[np.where(X_train['ss3_seq']==f[1])]=f[2]    
    
##
    
shap.summary_plot(shap_values, X_train ,max_display=8,show=False)
plt.savefig('../../Figures/Figure5/B - feature_importance_individual.png')
plt.show()

shap_reference = shap_values[:,[0,1,5,6]]
features_reference = X_train.values
features_reference = features_reference[:,[0,1,5,6]]
filename = 'shap_values_isoforms/shap_reference.csv'
np.savetxt(filename,shap_reference)
filename = 'shap_values_isoforms/features_reference.csv'
np.savetxt(filename,features_reference)
