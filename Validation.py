# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 21:59:46 2024

@author: iuaa
"""

import numpy as np
import pandas as pd
# #%%
# data = np.load('tape_total_layer0.npz', allow_pickle=True)

# # Access the array stored in the file
# emb = data['emb']

# # Perform operations on the loaded array
# # print(emb.shape)


# df_embed = pd.DataFrame(emb)

# # dfm = pd.concat([df_tel, df_embed], axis=1)

# print(df_embed)

# #%%

# df_new = pd.concat([df_tel, df_embed], axis=1)

#%% LOAD model and mappings
import joblib

# Save the model to a file
# joblib.dump(model, 'bestfit_model.joblib')
loaded_model = joblib.load('bestfit_model2025.joblib')
model = loaded_model

# Load the mappings from the file
# import joblib
loaded_mappings = joblib.load('all_mappings.joblib')

host_mapping = loaded_mappings['host_mapping']
epitope_mapping = loaded_mappings['species_mapping']
# method_mapping = loaded_mappings['method_mapping']

#%% Open the file

path24 = r"D:\Project_IEDB\newdata\FINAL\test\Tcell2024.xlsx"
# path = r"D:\Project_IEDB\newdata\FINAL\test\data2324.xlsx"

df24 = pd.read_excel(path24, sheet_name='Sheet1', 
                    usecols=lambda column: column.startswith('Assay') or column.startswith('Epitope')
                    or column.startswith('Host'))

# df = pd.read_excel('test_2024.xlsx', sheet_name='Sheet1', 
#                     usecols=lambda column: column.startswith('Assay') or column.startswith('Epitope')
#                     or column.startswith('Host'))

#%% read and clean file

dfing = df24[df24['Assay - Response measured'] == 'IFNg release'].reset_index(drop=True)
df_30 = clean(dfing)

#%% process data

df_asli = df_30.drop_duplicates(ignore_index=True)

df_2024 = sort_check(df_asli) # species and host

df_decode = df_2024

#%% Hamming distance

# Example column names
encode_seqs = df_encode['Epitope - Name']
decode_seqs = df_decode['Epitope - Name']

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

#%% visualise

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(5, 3), dpi=600)
sns.histplot(diff['hamming_distance'], bins=range(diff['hamming_distance'].min(), diff['hamming_distance'].max()+1), kde=True)
# plt.title("Distribution of Minimum Hamming Distances to Validation Set")
plt.xlabel("Minimum Hamming Distance")
plt.ylabel("Number of Sequences")
plt.tight_layout()
plt.show()


#%%

common_epitopes = pd.merge(df_encode[['Epitope - Name']], df_decode[['Epitope - Name']], on='Epitope - Name')

common_epitope_names = common_epitopes['Epitope - Name'].tolist()

# Filter out common rows from df_decode
df_filter = df_decode[~df_decode['Epitope - Name'].isin(common_epitope_names)].reset_index(drop=True)

# Display the filtered DataFrame
print(df_filter)

dfcom = df_decode[df_decode['Epitope - Name'].isin(common_epitope_names)].reset_index(drop=True)
#173 common epitopes

#%% CORRECT
# cols = [ 'Epitope ID', 'Specie index', 'Method index', 'Host index']
cols = ['Epitope ID', 'Epitope - Name', 'Epitope - Species', 'Host - Name']#,'Assay - Method']
unique_combinations = set(
    tuple(row[cols]) 
    for index, row in df_encode.iterrows())

# unique_combinations = [tuple(row[cols]) for index, row in df_train.iterrows()]


df_test2 = df_2024[~df_2024[cols]
                    .apply(tuple, axis=1).isin(unique_combinations)].reset_index(drop=True)

print(df_test2)

#%%
df_test3 = df_2024[df_2024[cols]
                    .apply(tuple, axis=1).isin(unique_combinations)].reset_index(drop=True)

print(df_test3)

df_decode = df_test3

#%% test set

df_decode = df_test2  ## apply the encodings after this

df_decode = df_filter ## apply the encodings after this

df_decode = df_2024

df_decode = dfcom # common

df_decode = df_epi2.reset_index(drop=True) #mouse/human host of 2024

# df_decode = dft
# df_decode = df_decode.drop('length', axis=1)
################  apply the encodings after this

#%%
df_decode['Host index'] = df_decode['Host - Name'].map(host_mapping).fillna(43)

df_decode['Specie index'] = df_decode['Epitope - Species'].map(epitope_mapping).fillna(121)

epit_test = df_decode['Epitope - Name']

dec_epitope = enc(epit_test)  #ORDINAL
# dec_epitope = dpc(epit_test)  #DPC
# dec_epitope = enc_phys(epit_test) # PHYSIOCHEMICAL

# decoding = np.concatenate((dec_epitope, dec_phys), axis=1)

dec_df = pd.DataFrame(dec_epitope, columns = [f'{i+1}' for i in range(dec_epitope.shape[1])])
# f'enc_{i}'
# dec_df = pd.DataFrame(dec_epitope, columns=column_enc)

# df_test = pd.concat([encoded_dftest, dec_df], axis=1)
df_test = pd.concat([df_decode, dec_df], axis=1)

# df_test['Length'] = df_test['Epitope - Name'].str.len()

print(df_test)

#%%
target = 'Label'
features = df_test.columns[5:].tolist()

test_set = df_test[features].values
# test_set = dft20[features].values

label_test = df_test[target].values
# label_test = dft20[target].values

split = len(label_test)/(len(label)+len(label_test))
print('Test train split is', split)

 #%%make predictions using model

y_pred = model.predict(test_set)
accuracy = accuracy_score(label_test, y_pred)

print( "Accuracy of extremeBoost: ", accuracy)

from sklearn.metrics import classification_report
classification_report = classification_report(label_test, y_pred)

print(classification_report)

# mcc_score = matthews_corrcoef(label_test, y_pred)
# print('MCC score is:',mcc_score)

evaluate_model(y_pred, label_test)

conf_matrix = confusion_matrix(label_test, y_pred)
print(conf_matrix)

# 839 samples df_decode_filtered
# Sensitivity (Recall): 0.8525
# Specificity: 0.9225512528473804
# Accuracy is:  0.8891537544696066
# F1 score is: 0.88
# MCC score is: 0.778606552483698

#%%

#CORRECT PREDICTIONS
dfcheck = df_test[df_test['Label']==df_test['prediction']]

dfinc = df_test[df_test['Label']!=df_test['prediction']] # incorrect predictions



#%% INCORRECT
# col = [ 'Epitope - Name', 'Epitope - Species', 'Host - Name','Assay - Method']

# df_test3 = df_decode[~df_decode[col].apply(lambda row: row.equals(df_encode[col]), axis=1)].reset_index(drop=True)

# df_test3

#%%

xx = remove_duplicate(df_2024, df_train)

#%% CORRECT

cols = ['Epitope ID',  'Epitope - Species', 'Host - Name','Assay - Method']

merged_df = pd.merge(df_2024, df_encode, on=cols, how='left', indicator=True)
df_test4 = merged_df[merged_df['_merge'] == 'left_only'].drop(columns='_merge').reset_index(drop=True)

#%%

filtered_df = df_tel[df_tel['Epitope - Species'].notna() & df_tel['Epitope - Species'].str.contains('ebola', na=False)]

# filtered_df['Epitope - Species'].value_counts()
# Out[454]: 
# Epitope - Species
# Zaire ebolavirus    25
# Sudan ebolavirus    24


# epitope_species[:10]
# Out[457]: 
# Epitope - Species
# Roseolovirus humanbeta6b                           5150
# Dengue virus                                       2405
# Vaccinia virus                                     2260
# Severe acute respiratory syndrome coronavirus 2    1752
# Orthoflavivirus zikaense                           1722
# Mycobacterium tuberculosis                         1689
# Homo sapiens                                       1538
# Alphainfluenzavirus influenzae                     1433
# Orthoflavivirus nilense                            1229
# Phleum pratense                                    1211

#%% ifnepitope 

a = df_2024
# a = df_decode_filtered
web = a.copy(deep=True)
webs = web.iloc[:,:3]
webs = webs.drop_duplicates(ignore_index=True)

df_web = sort_amino(webs)

#%% human host for webserver

from Bio import SeqIO

from io import StringIO

# df_web.to_excel('df_web.xlsx', index=False)
def dataframe_to_fasta(dataframe, id_col='Epitope ID', seq_col='Epitope - Name'):
    fasta_lines = []
    for _, row in dataframe.iterrows():
        fasta_lines.append(f'>{row[id_col]}\n{row[seq_col]}')
        # fasta_lines.append(f'>Epitope_{row[id_col]}\n{row[seq_col]}')
    return '\n'.join(fasta_lines)


# def dataframe_to_fasta(dataframe, id_col='Epitope ID', seq_col='Epitope - Name'):
#     fasta_lines = []
#     for index, row in dataframe.iterrows():
#         epitope_number = index + 1
#         # fasta_lines.append(f'>{row[id_col]}\n{row[seq_col]}')
#         fasta_lines.append(f'>Epitope_{epitope_number}\n{row[seq_col]}')
#     return '\n'.join(fasta_lines)


#%%
# Convert DataFrame to FASTA format
fasta_content = dataframe_to_fasta(df_web)

# Print or save the FASTA content
# print(fasta_content)

with StringIO(fasta_content) as fasta_file:
    SeqIO.write(SeqIO.parse(fasta_file, "fasta"), "df24.fasta", "fasta")

#%%

dfe = pd.read_excel('df_web.xlsx', sheet_name='Sheet4')

dfe = dfe.sort_values(by = 'id')
# dfc['length'] = dfc['Seq'].str.len()
# dfc.set_index('ID', inplace=True)

dfe['Label'] = dfe['result'].map({'NEGATIVE':0, 'POSITIVE':1})
y_pred = dfe['Label'].values
# y_label = df_web['Label'].values

y_label = df_hum['Label'].values #mouse host in df_hum
## OR
# dfc= pd.read_excel('df_web.xlsx', sheet_name = 'Sheet3')                   
# dfc['Label'] = dfc['result'].map({'NEGATIVE':0, 'POSITIVE':1})
# y_pred = dfc['Label'].values
# y_label = df_web['Label'].values

accuracy = accuracy_score(y_label, y_pred)

print( "Accuracy of webserver: ", accuracy)

from sklearn.metrics import classification_report
classification_report = classification_report(y_label, y_pred)

print(classification_report)

evaluate_model(y_pred, y_label)    

#%% ifnepitope
# 792 epitopes
# Sensitivity (Recall): 0.43869209809264303
# Specificity: 0.6847058823529412
# Accuracy is:  0.5707070707070707
# F1 score is: 0.48640483383685795
# MCC score is: 0.12727594984428733

#%% IFNepitope2

df_epi2 = df_2024[df_2024['Host - Name'].str.startswith('Mus musculus')].reset_index(drop=True)

# df_epi2 = df_2024[df_2024['Host - Name'] == 'Homo sapiens'].reset_index(drop=True)

dft = df_epi2
# dft = df_decode_filtered[df_decode_filtered['Host index'] == 2].reset_index(drop=True)

# dft = df_decode[df_decode['Host - Name'] == 'Homo sapiens'].reset_index(drop=True)

dft['length'] = dft['Epitope - Name'].str.len()

dft20 = dft[(12 <= dft['length']) & (dft['length'] <= 20)].reset_index(drop=True)

hum = dft20.iloc[:,:3]
hum = hum.drop_duplicates(ignore_index=True)

df_hum = sort_amino(hum)

#%%

fasta_content = dataframe_to_fasta(df_hum)

# Print or save the FASTA content
# print(fasta_content)

with StringIO(fasta_content) as fasta_file:
    SeqIO.write(SeqIO.parse(fasta_file, "fasta"), "dfhuman2324.fasta", "fasta")
    
#%% ifnepitope2 predictions on mouse/human

dfc = pd.read_excel('data2024org.xlsx', sheet_name='Sheet3') # sheet2-mouse, sheet3-human # sheet4 2324human prediction

dfc['predict'] = dfc['Prediction'].map({'Non-inducer':0, 'IFN-γ inducer':1})
y_pred = dfc['predict'].values
y_label = df_hum['Label'].values

accuracy = accuracy_score(y_label, y_pred)

print( "Accuracy of webserver: ", accuracy)

from sklearn.metrics import classification_report
classification_report = classification_report(y_label, y_pred)

print(classification_report)

evaluate_model(y_pred, y_label)   

conf_matrix = confusion_matrix(y_label, y_pred)
print(conf_matrix)

#%%


dfc = pd.read_csv('predict_dft20.csv')
# dfc['length'] = dfc['Seq'].str.len()
# dfc.set_index('ID', inplace=True)

dfc['Label'] = dfc['Hybrid_Score'].map({'Non-inducer':0, 'IFN-γ inducer':1})
y_pred = dfc['Label'].values
y_label = df_hum['Label'].values
## OR
# dfc= pd.read_excel('df_web.xlsx', sheet_name = 'Sheet3')                   
# dfc['Label'] = dfc['result'].map({'NEGATIVE':0, 'POSITIVE':1})
# y_pred = dfc['Label'].values
# y_label = df_web['Label'].values

accuracy = accuracy_score(y_label, y_pred)

print( "Accuracy of webserver: ", accuracy)

from sklearn.metrics import classification_report
classification_report = classification_report(y_label, y_pred)

print(classification_report)

evaluate_model(y_pred, y_label)        

#%% Generalizability and Comparison

import matplotlib.pyplot as plt

default_weight = plt.rcParams['font.weight']

plt.rcParams.update({
    'font.size': 8,
    'font.family': 'Arial',
    'font.weight': 'normal',
    # Add other font properties as needed
})

#%% generated same figure and separated in powerpoint
import seaborn as sns
import numpy as np

categories = ['IFNBoost \n(All)', 'IFNBoost \n(Unseen)', 'IFNepitope', 'IFNepitope2'] # all 1012 samples and unseen 839 samples
#2024 independent data  1012 samples IFNBoost, 944 IFNepitope, 854 IFNepitope2
# remove assay - 970 samples and 822 (unseen)

sens = [0.778, 0.776,0.428, 0.729]
spf = [0.914, 0.916, 0.687, 0.472]
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

#%% ifnepitope and 2

# # # 2024 independent data  1012 samples IFNBoost, 944 IFNepitope, 854 IFNepitope2

# categories = ['IFNepitope', 'IFNepitope2']

# sens = [ 0.428, 0.729]
# spf = [ 0.687, 0.472]
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