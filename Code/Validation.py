# -*- coding: utf-8 -*-
"""

@author: iuaa
"""

import numpy as np
import pandas as pd
from functions import *

#%% LOAD model and mappings
import joblib

loaded_model = joblib.load('ifnboost_trained.joblib')
model = loaded_model

# Load the mappings from the file
# import joblib
loaded_mappings = joblib.load('all_mappings.joblib')

host_mapping = loaded_mappings['host_mapping']
epitope_mapping = loaded_mappings['species_mapping']
# method_mapping = loaded_mappings['method_mapping']

#%% Open the file

path24 = "Tcell2024.xlsx" # path to 2024 data file

df24 = pd.read_excel(path24, sheet_name='Sheet1', 
                    usecols=lambda column: column.startswith('Assay') or column.startswith('Epitope')
                    or column.startswith('Host'))


#%% read and clean file

dfing = df24[df24['Assay - Response measured'] == 'IFNg release'].reset_index(drop=True)
df_30 = clean(dfing)

#%% process data

df_asli = df_30.drop_duplicates(ignore_index=True)

df_2024 = sort_check(df_asli) # species and host

df_decode = df_2024

#%% common epitopes in 2024 and training data #df_encode from MODEL.py

common_epitopes = pd.merge(df_encode[['Epitope - Name']], df_decode[['Epitope - Name']], on='Epitope - Name')

common_epitope_names = common_epitopes['Epitope - Name'].tolist()

# Filter out common rows from df_decode
df_filter = df_decode[~df_decode['Epitope - Name'].isin(common_epitope_names)].reset_index(drop=True)

# Display the filtered DataFrame
# print(df_filter)

dfcom = df_decode[df_decode['Epitope - Name'].isin(common_epitope_names)].reset_index(drop=True)
#173 common epitopes

#%% test set...choose 1 of these and test

df_decode = df_2024 # 2024 data -> 970 samples

# df_decode = df_filter # 2024 data not in training data -> 822 samples

#%% encoding the test set

df_decode['Host index'] = df_decode['Host - Name'].map(host_mapping).fillna(43)

df_decode['Specie index'] = df_decode['Epitope - Species'].map(epitope_mapping).fillna(121)

epit_test = df_decode['Epitope - Name']

dec_epitope = enc(epit_test)  #ORDINAL

dec_df = pd.DataFrame(dec_epitope, columns = [f'{i+1}' for i in range(dec_epitope.shape[1])])

df_test = pd.concat([df_decode, dec_df], axis=1)


#%%
target = 'Label'
features = df_test.columns[5:].tolist()

test_set = df_test[features].values

label_test = df_test[target].values

# split = len(label_test)/(len(label)+len(label_test))
# print('Test train split is', split)

 #%%make predictions using model

y_pred = model.predict(test_set)

evaluate_model(y_pred, label_test)



#%% Hamming distance

# Example column names
encode_seqs = df_encode['Epitope - Name']
decode_seqs = df_2024['Epitope - Name']

# Function to calculate Hamming distance
def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

# Store results
results = []

# Iterate through each sequence in the training set
for i, seq1 in enumerate(encode_seqs):
    for j, seq2 in enumerate(decode_seqs):
        if len(seq1) == len(seq2):
            dist = hamming_distance(seq1, seq2)
            results.append({
                'encode_index': i,
                'decode_index': j,
                'seq1': seq1,
                'seq2': seq2,
                'hamming_distance': dist
            })

# Convert to a DataFrame for easy exploration
df_hamming = pd.DataFrame(results)

# Optional: find closest matches
df_closest = df_hamming.sort_values(['encode_index', 'hamming_distance']).groupby('encode_index').first().reset_index()

# Get minimum Hamming distance for each encode sequence
min_dists = df_hamming.groupby(['decode_index', 'seq2'])['hamming_distance'].min().reset_index()

#Unseen epitopes
diff = min_dists[min_dists['hamming_distance']!=0]

# min_dists.to_excel('hamming.xlsx')

#%% visualise the hamming distances

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(5, 3), dpi=600)
sns.histplot(diff['hamming_distance'], bins=range(diff['hamming_distance'].min(), diff['hamming_distance'].max()+1), kde=True)
# plt.title("Distribution of Minimum Hamming Distances to Validation Set")
plt.xlabel("Minimum Hamming Distance")
plt.ylabel("Number of Sequences")
plt.tight_layout()
plt.show()


#%% Creating data for ifnepitope 

a = df_2024
# a = df_decode_filtered
web = a.copy(deep=True)
webs = web.iloc[:,:3]
webs = webs.drop_duplicates(ignore_index=True)

df_web = sort_amino(webs) #filter based on amino acids #944 epitopes

#%% function to convert to fasta

from Bio import SeqIO

from io import StringIO

# df_web.to_excel('df_web.xlsx', index=False)
def dataframe_to_fasta(dataframe, id_col='Epitope ID', seq_col='Epitope - Name'):
    fasta_lines = []
    for _, row in dataframe.iterrows():
        fasta_lines.append(f'>{row[id_col]}\n{row[seq_col]}')
        # fasta_lines.append(f'>Epitope_{row[id_col]}\n{row[seq_col]}')
    return '\n'.join(fasta_lines)


#%%
# Convert DataFrame to FASTA format
fasta_content = dataframe_to_fasta(df_web)

#save fasta
# with StringIO(fasta_content) as fasta_file:
#     SeqIO.write(SeqIO.parse(fasta_file, "fasta"), "df24.fasta", "fasta")

# df24.fasta was tested in IFNepitope webserved

#%% Test IFNepitope: df_web.xlsx contains results of IFNepitope

dfc= pd.read_excel('df_web.xlsx', sheet_name = 'Sheet3')       
dfc = dfc.sort_values(by='id') # sort by ID to compare            
dfc['Label'] = dfc['result'].map({'NEGATIVE':0, 'POSITIVE':1})
y_pred = dfc['Label'].values
y_label = df_web['Label'].values

evaluate_model(y_pred, y_label)    

#%% IFNepitope2 - Mouse and human host and 12<length<20

df_epi2_mouse = df_2024[df_2024['Host - Name'].str.startswith('Mus musculus')].reset_index(drop=True) #Mouse

df_epi2 = df_2024[df_2024['Host - Name'] == 'Homo sapiens'].reset_index(drop=True) #Human

#HUMAN
dft = df_epi2

dft['length'] = dft['Epitope - Name'].str.len()

dft20 = dft[(12 <= dft['length']) & (dft['length'] <= 20)].reset_index(drop=True)

hum = dft20.iloc[:,:3]
hum = hum.drop_duplicates(ignore_index=True)

df_hum = sort_amino(hum) # human

#MOUSE
dftm = df_epi2_mouse

dftm['length'] = dftm['Epitope - Name'].str.len()

dft20m = dftm[(12 <= dftm['length']) & (dftm['length'] <= 20)].reset_index(drop=True)

mouse = dft20m.iloc[:,:3]
mouse = mouse.drop_duplicates(ignore_index=True)

df_mouse = sort_amino(mouse) # mouse

#%% convert to fasta

# fasta_content = dataframe_to_fasta(df_hum)

# with StringIO(fasta_content) as fasta_file:
#     SeqIO.write(SeqIO.parse(fasta_file, "fasta"), "dfhuman.fasta", "fasta") #dfmouse
    
#%% ifnepitope2 predictions on human/mouse

dfc = pd.read_excel('data2024org.xlsx', sheet_name='Sheet3') # sheet2-mouse, sheet3-human # sheet4 2324human prediction

dfc['predict'] = dfc['Prediction'].map({'Non-inducer':0, 'IFN-γ inducer':1})
y_pred = dfc['predict'].values
y_label = df_hum['Label'].values 

# evaluate_model(y_pred, y_label)   
conf_matrix_hum = confusion_matrix(y_label, y_pred)
print("Confusion matrix for IFNepitope2 human prediction: \n", conf_matrix_hum)

#MOUSE results
dfc = pd.read_excel('data2024org.xlsx', sheet_name='Sheet2') # sheet2-mouse, sheet3-human # sheet4 2324human prediction
dfc['predict'] = dfc['Prediction'].map({'Non-inducer':0, 'IFN-γ inducer':1})
y_pred = dfc['predict'].values
y_label = df_mouse['Label'].values 

# evaluate_model(y_pred, y_label)   
conf_matrix_mouse = confusion_matrix(y_label, y_pred)
print("Confusion matrix for IFNepitope2 mouse prediction: \n", conf_matrix_mouse)

#%% Aggregate confusion matrix

# 709 mouse samples
# [[231 256]
# [ 85 137]]

# 145 human samples
# [[  0   2]
#  [ 14 129]]

# Aggregate
# [[ 231 258]
#  [ 99  266]]

# 854 samples total
# Sensitivity (Recall): 0.729
# Specificity: 0.472
# Accuracy is:  0.582
# F1 score is: 0.598
# MCC score is: 0.204


#%% Generalizability and Comparison

import matplotlib.pyplot as plt

default_weight = plt.rcParams['font.weight']

plt.rcParams.update({
    'font.size': 8,
    'font.family': 'Arial',
    'font.weight': 'normal',
    # Add other font properties as needed
})

#%% generated figures
import seaborn as sns
import numpy as np

categories = ['IFNBoost \n(All)', 'IFNBoost \n(Unseen)', 'IFNepitope', 'IFNepitope2'] # all 1012 samples and unseen 839 samples
#2024 independent data  1012 samples IFNBoost, 944 IFNepitope, 854 IFNepitope2
# remove assay - 970 samples and 822 (unseen)

sens = [0.778, 0.776,0.546, 0.729]
spf = [0.914, 0.916, 0.578, 0.472]
accuracy = [0.840, 0.8418, 0.566, 0.582]
f1 = [0.841, 0.8395, 0.480, 0.598]
mcc = [0.691, 0.6945, 0.120, 0.204]


# Create a matrix with MCC and Accuracy
data = np.array([sens, spf,accuracy, f1, mcc])

# Create a DataFrame
import pandas as pd
dfo = pd.DataFrame(data, columns=categories, index=['Sensitivity', 'Specificity','Accuracy', 'F1 score', 'MCC'])

# Create a heatmap using seaborn
plt.figure(figsize=(4.5,3.5), dpi=600)
# colormap = sns.color_palette("Greens")
colormap = sns.light_palette('#115316')
heatmap = sns.heatmap(dfo, annot=True, cmap=colormap, fmt=".3f", linewidths=1,annot_kws={"size": 8},square=True, vmin=0,vmax=1,cbar=True)   

heatmap.set_xticklabels(heatmap.get_xticklabels(), fontsize=8 , rotation=30, ha='center')
heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=8,rotation=0)

plt.tight_layout()
plt.show()

#%% ifnepitope and ifnepitope2 figures

# # # 2024 independent data  1012 samples IFNBoost, 944 IFNepitope, 854 IFNepitope2

# categories = ['IFNepitope', 'IFNepitope2']

# sens = [ 0.546, 0.729]
# spf = [ 0.578, 0.472]
# accuracy = [ 0.566, 0.582]
# f1 = [ 0.480, 0.598]
# mcc = [ 0.120, 0.204]


# # Create a matrix with MCC and Accuracy
# data = np.array([sens, spf,accuracy, f1, mcc])

# # Create a DataFrame
# import pandas as pd
# dfo = pd.DataFrame(data, columns=categories, index=['Sensitivity', 'Specificity','Accuracy', 'F1 score', 'MCC'])

# # Create a heatmap using seaborn
# plt.figure(figsize=(3,3), dpi=300)
# # colormap = sns.color_palette("Greens")
# colormap = sns.light_palette('#115316')
# heatmap = sns.heatmap(dfo, annot=True, cmap=colormap, fmt=".3f", linewidths=1,annot_kws={"size": 8},square=True, vmin=0,vmax=1)   

# heatmap.set_xticklabels(heatmap.get_xticklabels(), fontsize=8 , rotation=30, ha='center')
# # heatmap.set_yticklabels([])
# # heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=8,rotation=0)

# plt.tight_layout()
# plt.show()