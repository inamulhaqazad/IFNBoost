# -*- coding: utf-8 -*-
"""
@author: iuaa
"""

import pandas as pd
import numpy as np
import random
SEED = 42
np.random.seed(SEED)
random.seed(SEED)

import os
os.environ['PYTHONHASHSEED']= str(SEED)

def clean(df):
    selected_columns = [
        'Epitope - Name',
        'Epitope - IEDB IRI',
        'Assay - Qualitative Measurement',
        'Epitope - Species',
        'Host - Name'
          
    ]

    # Create a new DataFrame with selected columns
    new_df = df.loc[:, selected_columns] 
    
    new_df['Label'] = new_df['Assay - Qualitative Measurement'].map\
                               ({'Negative': 0, 'Positive': 1, 'Positive-Low': 1, 
                                 'Positive-Intermediate': 1,'Positive-High': 1})

    new_df['Label'].value_counts()     

    new_df = new_df.drop('Assay - Qualitative Measurement', axis=1)


    # remove inavlid epitopes from data

    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    # Define the condition
    condition = lambda row: any(letter not in amino_acids for letter in row['Epitope - Name'])

    filtered_df = new_df[new_df.apply(condition, axis=1)]


    df_minus_filtered = new_df.drop(filtered_df.index).reset_index(drop=True)

    #print(df_minus_filtered)
    df_minus_filtered['Epitope ID'] = df_minus_filtered['Epitope - IEDB IRI'].str.split('/').str[-1]

    df_minus_filtered["Epitope ID"] = pd.to_numeric(df_minus_filtered["Epitope ID"])

    # Sort the DataFrame by the "Epitope ID" column
    dclean = df_minus_filtered.sort_values(by="Epitope ID", ignore_index=True)
    dclean.drop('Epitope - IEDB IRI', axis=1, inplace=True)

    del new_df, filtered_df, df_minus_filtered

    # keep epitopes of under 30 length

    leng = dclean['Epitope - Name'].str.len()

    df_30 = dclean[(12<=leng) & (leng<= 30)].reset_index(drop=True)
    
    columns_ordered = ['Epitope ID', 'Epitope - Name', 'Label',
                       'Epitope - Species', 'Host - Name']#, 'Assay - Method']
    df_30 = df_30[columns_ordered]
    
    df_30['Host - Name'] = df_30['Host - Name'].apply(combine_homo_sapiens)
    return df_30
#%%
"sort_initial checks for samples with positive and negative IFNg response and if an epitope has \
 both positive and negative response, it is put in positive labels"
 
def sort_initial(df):
    cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Assay - Method', 'Host - Name']
    # cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Assay - Method', 'Host Species']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_amino(df):
   
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]
     
    df_n = df_n[~df_n['Epitope - Name'].isin(df_p['Epitope - Name'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    return df_

#%%

def sort_host(df):
    cols = ['Epitope - Name', 'Host - Name']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_specie(df):
    cols = ['Epitope - Name', 'Epitope - Species']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_assay(df):
    cols = ['Epitope - Name', 'Assay - Method']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_check(df):
    cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Host - Name']
    # cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Host Species']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_assayhost(df):
    cols = ['Epitope - Name', 'Epitope ID', 'Assay - Method', 'Host - Name']
    # cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Host Species']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

def sort_assayspecie(df):
    cols = ['Epitope - Name', 'Epitope ID', 'Assay - Method', 'Epitope - Species']
    # cols = ['Epitope - Name', 'Epitope ID', 'Epitope - Species', 'Host Species']
    df_p = df[df['Label'] == 1]
    df_n = df[df['Label'] == 0]

    df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
    df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
    df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
    df_ = df_.drop(columns = ['nc'])
    return df_

#%% evaluation metrics

from sklearn.metrics import classification_report, matthews_corrcoef, f1_score, confusion_matrix, accuracy_score

def evaluate_model(la_test, y_pred):
    # y_pred = model.predict(test_set)
    
    # # Confusion Matrix
    # conf_matrix = confusion_matrix(la_test, y_pred)
    
    # # Sensitivity (Recall) and Specificity
    # sensitivity = conf_matrix[1, 1] / (conf_matrix[1, 0] + conf_matrix[1, 1])
    # specificity = conf_matrix[0, 0] / (conf_matrix[0, 0] + conf_matrix[0, 1])
    
    # Confusion Matrix
    tn, fp, fn, tp = confusion_matrix(la_test, y_pred).ravel()
    
    # Sensitivity (Recall) and Specificity
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    
    print('Sensitivity (Recall):', sensitivity)
    print('Specificity:', specificity)
          
    # Accuracy
    accuracy = accuracy_score(la_test, y_pred)
    print("Accuracy is: ", accuracy)
    
    # F1 Score
    f1 = f1_score(la_test, y_pred)
    print('F1 score is:', f1)
    
    # Classification Report
    # class_report = classification_report(la_test, y_pred)
    # print(class_report)
    
    # Matthews Correlation Coefficient (MCC)
    mcc_score = matthews_corrcoef(la_test, y_pred)
    print('MCC score is:', mcc_score)

#%%

def combine_homo_sapiens(host):
    if 'Homo sapiens' in host:
        return 'Homo sapiens'
    else:
        return host

def mapping(df):
    hosts = df['Host - Name'].value_counts()
    top_hosts = hosts.head(43).index
    other_hosts = hosts.index.difference(top_hosts)
    
    host_mapping = {host: i + 1 for i, host in enumerate(top_hosts)}
    for host in other_hosts:
        host_mapping[host] = 44    
    
    # hosts = df['Host Species'].value_counts()
    # host_mapping = {host:i+1 for i, host in enumerate(hosts)}
    
    # Encode 'Epitope - Species' column
    epitope_species = df['Epitope - Species'].value_counts()
    top_specie = epitope_species.head(122).index
    other_specie = epitope_species.index.difference(top_specie)
    
    epitope_mapping = {epitope: i + 1 for i, epitope in enumerate(top_specie)}
    for epitope in other_specie:
        epitope_mapping[epitope] = 123
        
    # Encode 'Assay - Method' column
    # method = df['Assay - Method'].value_counts()
    
    # method_mapping = {method: i + 1 for i, method in enumerate(method.index)}
    
    return host_mapping, epitope_mapping #, method_mapping
    
#%% encoding

def encode_categorical_columns(df, host_map, epitope_map, method_map):
    
    df['Host index'] = df['Host - Name'].map(host_map).fillna(49)
        
    df['Specie index'] = df['Epitope - Species'].map(epitope_map).fillna(128)
        
    df['Method index'] = df['Assay - Method'].map(method_map).fillna(0)
    
    return df

#%% return features and labels from df

def encoding(df):

    epitopes = df['Epitope - Name']

    enc_epitope = enc(epitopes)
    enc_df = pd.DataFrame(enc_epitope, columns=[f'enc_{i}' for i in range(enc_epitope.shape[1])])

    target = 'Label'
    enc_df = pd.concat([df, enc_df], axis=1)
    
    features = enc_df.columns[6:].tolist()

    train_set = enc_df[features].values
    labels = enc_df[target].values

    return train_set, labels



#%%
def remove_duplicate(dftest, dftrain):
    a = dftest.copy(deep=True)
    b = dftrain.copy(deep=True)
    cols = ['Epitope - Name', 'Specie index', 'Method index', 'Host index']
    
    # Create a new column 'nc' in both dataframes
    b['nc'] = b.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
    a['nc'] = a.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)

    # Merge dataframes based on the 'nc' column
    merged_df = pd.merge(a, b[['nc']], how='left', on='nc', indicator=True)

    # Keep only the rows that are unique to the left dataframe
    a_unique = merged_df[merged_df['_merge'] == 'left_only'].drop(columns=['_merge', 'nc'])

    return a_unique

# #%%


# def sort(df):
#     cols = ['Epitope - Name', 'Epitope ID', 'Specie index', 'method index', 'Host index']
#     df_p = df[df['Label'] == 1]
#     df_n = df[df['Label'] == 0]

#     df_p.loc[:, 'nc'] = df_p.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
#     df_n.loc[:, 'nc'] = df_n.apply(lambda x: ''.join(str(i) for i in x[cols]), 1)
     
#     df_n = df_n[~df_n['nc'].isin(df_p['nc'])]
     
#     df_ = pd.concat([df_p, df_n]).sort_values(by = 'Epitope ID').reset_index(drop = True)
     
#     df_ = df_.drop(columns = ['nc'])
#     return df_

# def remove_dupli(dftest, dftrain):
#     a = dftest.copy(deep=True)
#     b = dftrain.copy(deep=True)
#     cols = ['Epitope - Name', 'Specie index', 'method index', 'Host index']
#     b.loc[:,'nc'] = b.apply(lambda x: ''.join(str(i) for i in x[cols]),1)
#     a.loc[:,'nc'] = a.apply(lambda x: ''.join(str(i) for i in x[cols]),1)
    
#     a = a[~a['nc'].isin(b['nc'])].reset_index(drop=True)
    
#     a = a.drop(columns = ['nc'])
    
#     return a


#%%    
aminos = ['X', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

amino_acid_encoding = {amino_acid: index  for index, amino_acid in enumerate(aminos)}

def pad(seq, length):
    n = len(seq)
    if n>length:
        seq = seq[:length]
    elif n<length:
        seq = seq + ('X' * (length-n))
    return seq    

def enc(sequence):
    emb_list = [pad(seq, 30) for seq in sequence]
    #print(emb_list)
    a=[]
    for seq in emb_list:
        encode = [amino_acid_encoding[a] for a in seq]
        a.append(encode)
    
    return np.asarray(a)    

#%%
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

dipeptide_alphabet = [aa1 + aa2 for aa1 in amino_acids for aa2 in amino_acids]

# epit = dfinal['Epitope Name']

def calculate_dpc(epitope):
    # Initialize an array to store the DPC values
    dpc = np.zeros(400)

    # Count the occurrences of each dipeptide in the epitope
    total_dipeptides = 0
    for i in range(len(epitope) - 1):
        dipeptide = epitope[i:i+2]
        if dipeptide in dipeptide_alphabet:
            index = dipeptide_alphabet.index(dipeptide)
            dpc[index] += 1
            total_dipeptides += 1

    # Calculate the DPC values as percentages
    if total_dipeptides > 0:
        dpc = (dpc / total_dipeptides) * 100
        
    #print(total_dipeptides)
    return dpc

def dpc(epitopes):
    x=[]
    for i in range(len(epitopes)):
        seq = epitopes[i]
        emb = calculate_dpc(seq)
        x.append(emb)
    return np.asarray(x)

# dp_epit = dpc(epit)

# dpc_df = pd.DataFrame(dp_epit, columns = [f'f_{i}' for i in range(dp_epit.shape[1])])

#%%

column_enc = columns=['Encoding\u2081',
                  'Encoding\u2082',
                  'Encoding\u2083',
                  'Encoding\u2084',
                  'Encoding\u2085',
                  'Encoding\u2086',
                  'Encoding\u2087',
                  'Encoding\u2088',
                  'Encoding\u2089',
                  'Encoding\u2081\u2080',
                  'Encoding\u2081\u2081',
                  'Encoding\u2081\u2082',
                  'Encoding\u2081\u2083',
                  'Encoding\u2081\u2084',
                  'Encoding\u2081\u2085',
                  'Encoding\u2081\u2086',
                  'Encoding\u2081\u2087',
                  'Encoding\u2081\u2088',
                  'Encoding\u2081\u2089',
                  'Encoding\u2082\u2080',
                  'Encoding\u2082\u2081',
                  'Encoding\u2082\u2082',
                  'Encoding\u2082\u2083',
                  'Encoding\u2082\u2084',
                  'Encoding\u2082\u2085',
                  'Encoding\u2082\u2086',
                  'Encoding\u2082\u2087',
                  'Encoding\u2082\u2088',
                  'Encoding\u2082\u2089',
                  'Encoding\u2083\u2080']
                  
                  
            

#%%

def sort_old(df):
    
    pos_records = []
    neg_records = []

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
       # epitope_id = row['Epitope - IEDB IRI']
        epitope_id = row['Epitope ID']
        epitope_name = row['Epitope - Name']
        
        result = row['Label']
        
        host_name = row['Host index']
        
        epitope_species = row['Specie index']
        method = row['method index']
        
        

        # Check if assay response is 'IFNg' and quality is 'positive'
        if result == 1:
            #epitope_id = epitope_id.split('/')[-1]

            # Add the epitope ID to the list
            pos_records.append((epitope_id, epitope_name, host_name, epitope_species, method))
                               
        
        elif result == 0:
            #epitope_id = epitope_id.split('/')[-1]
            neg_records.append((epitope_id, epitope_name, host_name, epitope_species, method))
                               

    # POSITIVE IF ONLY 1 IS POSITIVE

    pos_records = list(set(pos_records))
    neg_records = list(set(neg_records))

    # Remove entries from negative_epitopes that are present in positive_epitopes
    neg_records = [epitope_id for epitope_id in neg_records if epitope_id not in pos_records]
    
    print(f'Positive epitopes {len(pos_records)}')
    print(f'Negative EPitopes {len(neg_records)}')

    #
    pos_df = pd.DataFrame(pos_records, columns=
                          ('Epitope ID', 'Epitope Name', 'Host', 'Epitope specie', 'assay method'))
                           
                           

    neg_df = pd.DataFrame(neg_records, columns=
                          ('Epitope ID', 'Epitope Name', 'Host', 'Epitope specie', 'assay method'))
                           
                          

    #
    pos_df['Label'] = 1
    neg_df['Label'] = 0 

    dfinal = pd.concat([pos_df, neg_df]).sort_values(by = ['Epitope ID']).reset_index(drop=True)

    dfinal = dfinal.sample(frac=1.0, random_state=SEED).reset_index(drop=True)


    #
    columns_ordered = ['Epitope ID', 'Epitope Name', 'Label',
                       'Epitope specie', 'Host', 'assay method']
                       

    # Reorder the columns based on the 'columns_ordered' list
    dfinal = dfinal[columns_ordered]
    
    return dfinal
#%% function to chose labels based on majority rule 
    # this has not been used in final model
def sort_majority(df):
    
    pos_records = []
    neg_records = []

    # Iterate over each row in the DataFrame
    for index, row in df.iterrows():
       # epitope_id = row['Epitope - IEDB IRI']
        epitope_id = row['Epitope ID']
        epitope_name = row['Epitope - Name']
        
        result = row['Label']
        
        host_name = row['Host index']
        
        epitope_species = row['Specie index']
        method = row['method index']
        
        

        # Check if assay response is 'IFNg' and quality is 'positive'
        if result == 1:

            # Add the epitope ID to the list
            pos_records.append((epitope_id, epitope_name, host_name, epitope_species, method))
                                
        
        elif result == 0:
            #epitope_id = epitope_id.split('/')[-1]
            neg_records.append((epitope_id, epitope_name, host_name, epitope_species, method))
                                
                                
    pos_counts = {}
    neg_counts = {}

    # Count positive records
    for record in pos_records:
        epitope_id = record[0]
        if epitope_id in pos_counts:
            pos_counts[epitope_id] += 1
        else:
            pos_counts[epitope_id] = 1

    # Count negative records
    for record in neg_records:
        epitope_id = record[0]
        if epitope_id in neg_counts:
            neg_counts[epitope_id] += 1
        else:
            neg_counts[epitope_id] = 1


    common_epitopes = (pos_counts.keys()) & (neg_counts.keys())


    positive_final = [(epitope_id, epitope_name, host_name, epitope_species, method) 
                      for epitope_id, epitope_name, host_name, epitope_species, method in pos_records if 
                      pos_counts.get(epitope_id,0) >= neg_counts.get(epitope_id, 0)]

    negative_final = [(epitope_id, epitope_name, host_name, epitope_species, method) 
                      for epitope_id, epitope_name, host_name, epitope_species, method in neg_records 
                      if neg_counts.get(epitope_id,0) > pos_counts.get(epitope_id, 0)]


    positive_final = list(set(positive_final))

    negative_final = list(set(negative_final))

    # Print the final positive and negative lists
    print(len(positive_final))
    print(len(negative_final))

    common_epitopes = []
    for epitope_id in negative_final:
        if epitope_id in positive_final:
            common_epitopes.append(epitope_id)

    print(len(common_epitopes))
        
    pos_df = pd.DataFrame(positive_final, columns=
                          ('Epitope ID', 'Epitope Name', 'Host', 'Epitope specie', 'assay method'))
                           
                           

    neg_df = pd.DataFrame(negative_final, columns=
                          ('Epitope ID', 'Epitope Name', 'Host', 'Epitope specie', 'assay method'))
                           

    #
    pos_df['Label'] = 1
    neg_df['Label'] = 0 

    dfinal = pd.concat([pos_df, neg_df]).reset_index(drop=True)

    dfinal = dfinal.sample(frac=1.0, random_state=SEED).reset_index(drop=True)


    #
    columns_ordered = ['Epitope ID', 'Epitope Name', 'Label',
                       'Epitope specie', 'Host', 'assay method']
                       

    # Reorder the columns based on the 'columns_ordered' list
    dfinal = dfinal[columns_ordered]
    
    return dfinal

#%%
import matplotlib.pyplot as plt
import numpy as np

def analyse(df_encode, df_decode, train, test):
    

    # Get the unique species names in both dataframes
    decode_species = df_decode['Epitope - Species'].unique()
    encode_species = df_encode['Epitope - Species'].unique()

    specie = df_encode['Epitope - Species'].value_counts().head(20).index.tolist()
    
    count = []
    for specs in specie:
        count.append(df_decode[df_decode['Epitope - Species']==specs].shape[0])

    count_encode = df_encode['Epitope - Species'].value_counts().head(20)
    
    plt.figure(figsize=(12, 8))

    # Create the bar plot for df_decode in blue
    plt.bar(np.arange(len(specie)), count, width=0.4, label= test, color='blue')

    # Create the bar plot for df_encode in orange
    plt.bar(np.arange(len(specie)) + 0.4, count_encode, width=0.4, label=train, color='orange')

    # Set the x-tick labels to be species names
    plt.xticks(np.arange(len(specie)) + 0.2, specie, rotation=90)

    # Set the title and labels
    plt.title('Distribution of Epitope - Species (top 20 species from training)')
    plt.xlabel('Epitope - Species')
    plt.ylabel('Count')

    # Show the legend
    plt.legend()

    # Show the plot
    plt.tight_layout()
    plt.show()
    
#%%
    
