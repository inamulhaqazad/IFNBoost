# -*- coding: utf-8 -*-
"""

@author: iuaa
"""

import pandas as pd
import numpy as np
import random
from functions import *

import sys

SEED = 42
np.random.seed(SEED)
random.seed(SEED)

import os
os.environ['PYTHONHASHSEED'] =str(SEED)
#%% read and preprocess

df = pd.read_excel('allhost.xlsx', sheet_name='Sheet1', 
                    usecols=lambda column: column.startswith('Assay') or column.startswith('Epitope')
                    or column.startswith('Host'))


dfing = df[df['Assay - Response measured'] == 'IFNg release'].reset_index(drop=True)
df_30 = clean(dfing)
#%% Process

df_asli = df_30.drop_duplicates(ignore_index=True)

df_tel = sort_check(df_asli) # host and epitope specie only


#%% train-test split

from sklearn.model_selection import train_test_split

df_enc, df_dec = train_test_split(df_tel, test_size = 0.20, stratify=df_tel['Label'], random_state=SEED)

df_encode = df_enc.reset_index(drop=True)
df_decode = df_dec.reset_index(drop=True)


#%% 2 input featutres : Metadata

hosts = df_encode['Host - Name'].value_counts() 

epitope_species = df_encode['Epitope - Species'].value_counts() 


ratio = epitope_species[:120].sum()/epitope_species.sum() # 120 99% afteer removing assay
# print('Epitope species', ratio)

ratio_hosts = hosts[:42].sum()/hosts.sum() # 42 99% after removing assay
# print('Hosts ratio', ratio_hosts)


#%% Mappings

hosts = df_encode['Host - Name'].value_counts()
top_hosts = hosts.head(42).index
other_hosts = hosts.index.difference(top_hosts)

host_mapping = {host: i + 1 for i, host in enumerate(top_hosts)}
for host in other_hosts:
    host_mapping[host] = 43   


# Encode 'Epitope - Species' column
epitope_species = df_encode['Epitope - Species'].value_counts()
top_specie = epitope_species.head(120).index
other_specie = epitope_species.index.difference(top_specie)

epitope_mapping = {epitope: i + 1 for i, epitope in enumerate(top_specie)}
for epitope in other_specie:
    epitope_mapping[epitope] = 121


#%% Missing entries go to others' category

df_encode['Host index'] = df_encode['Host - Name'].map(host_mapping).fillna(43)

df_encode['Specie index'] = df_encode['Epitope - Species'].map(epitope_mapping).fillna(121)


#%% Encode the training data

epit = df_encode['Epitope - Name']

enc_epitope = enc(epit) # ORDINAL ENCODING

enc_df = pd.DataFrame(enc_epitope, columns = [f'{i+1}' for i in range(enc_epitope.shape[1])])

df_train = pd.concat([df_encode, enc_df], axis=1)

# df_train['Length'] = df_train['Epitope - Name'].str.len()

# print(df_train)

#%% encoding of test set

df_decode['Host index'] = df_decode['Host - Name'].map(host_mapping).fillna(43)

df_decode['Specie index'] = df_decode['Epitope - Species'].map(epitope_mapping).fillna(121)

epit_test = df_decode['Epitope - Name']

dec_epitope = enc(epit_test)  #ORDINAL


dec_df = pd.DataFrame(dec_epitope, columns = [f'{i+1}' for i in range(dec_epitope.shape[1])])

df_test = pd.concat([df_decode, dec_df], axis=1)

# df_test['Length'] = df_test['Epitope - Name'].str.len()

# print(df_test)


#%% train test inputs for model


target = 'Label'

features = df_test.columns[5:].tolist()

train_set = df_train[features].values
label = df_train[target].values

test_set = df_test[features].values

label_test = df_test[target].values

split = len(label_test)/(len(label)+len(label_test))
# print('Test train split is', split)


#%% xgboost

import xgboost as xgb

model = xgb.XGBClassifier( # 
    n_estimators=1170,  # Number of boosting rounds (trees)
    learning_rate=0.05,  # Step size shrinkage
    max_depth=36,  # Maximum depth of each tree
     objective='binary:logistic',  # For binary classification
     seed = SEED # Evaluation metric for training
    ,reg_lambda=0.9, reg_alpha=0.5, gamma=2, min_child_weight = 2
    ,tree_method='auto', subsample=0.9, colsample_bytree=0.6
)


#%% evaluate model

model.fit(train_set, label)

from sklearn.metrics import accuracy_score, matthews_corrcoef, f1_score, confusion_matrix

y_pred = model.predict(test_set)

evaluate_model(label_test, y_pred) # from functions


#%% bootstrapping
import numpy as np
from sklearn.utils import resample
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef, confusion_matrix


#random seed for reproducibility
rng = np.random.RandomState(SEED)

# Number of bootstrap samples
n_bootstraps = 1000

acc_scores = []
sens_scores = []
spec_scores = []
f1_scores = []
mcc_scores = []

for _ in range(n_bootstraps):
    # Generate bootstrap sample indices
    indices = rng.choice(len(y_pred), len(y_pred), replace=True)

    # Get bootstrap samples
    y_true_bootstrap = np.array(label_test)[indices]
    y_pred_bootstrap = np.array(y_pred)[indices]

    # Compute metrics
    acc_scores.append(accuracy_score(y_true_bootstrap, y_pred_bootstrap))
    sens_scores.append(recall_score(y_true_bootstrap, y_pred_bootstrap))  # Sensitivity (Recall)
    
    # Compute specificity
    tn, fp, fn, tp = confusion_matrix(y_true_bootstrap, y_pred_bootstrap).ravel()
    spec_scores.append(tn / (tn + fp))  # Specificity = TN / (TN + FP)

    f1_scores.append(f1_score(y_true_bootstrap, y_pred_bootstrap))
    mcc_scores.append(matthews_corrcoef(y_true_bootstrap, y_pred_bootstrap))

# Compute 95% confidence intervals
def bootstrap_ci(scores):
    return np.percentile(scores, [2.5, 97.5])

acc_lower, acc_upper = bootstrap_ci(acc_scores)
sens_lower, sens_upper = bootstrap_ci(sens_scores)
spec_lower, spec_upper = bootstrap_ci(spec_scores)
f1_lower, f1_upper = bootstrap_ci(f1_scores)
mcc_lower, mcc_upper = bootstrap_ci(mcc_scores)

# Print results

print(f"Sensitivity (Recall): {np.mean(sens_scores):.4f} (95% CI: {sens_lower:.4f} - {sens_upper:.4f})")
print(f"Specificity: {np.mean(spec_scores):.4f} (95% CI: {spec_lower:.4f} - {spec_upper:.4f})")
print(f"Accuracy: {np.mean(acc_scores):.4f} (95% CI: {acc_lower:.4f} - {acc_upper:.4f})")
print(f"F1-score: {np.mean(f1_scores):.4f} (95% CI: {f1_lower:.4f} - {f1_upper:.4f})")
print(f"MCC: {np.mean(mcc_scores):.4f} (95% CI: {mcc_lower:.4f} - {mcc_upper:.4f})")


