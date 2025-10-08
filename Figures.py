"""
Model performance, AUROC AND PRC Curves with 95% CI

"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score
from sklearn.utils import resample
import seaborn as sns

#%% Predicted probabilities from model

y_pred_proba = model.predict_proba(test_set)[:,1]


#%%

# Assuming you have the following variables defined:
# label_test: true labels from MODEL.py file
# y_pred_proba: predicted probabilities

# Number of bootstrap samples
n_bootstraps = 1000
rng = np.random.RandomState(SEED)

# Arrays to store the results for ROC
tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)

# Arrays to store the results for PRC
precisions = []
recalls = []
aucs_prc = []
mean_recall = np.linspace(0, 1, 100)

for i in range(n_bootstraps):
    indices = rng.randint(0, len(y_pred_proba), len(y_pred_proba))
    if len(np.unique(label_test[indices])) < 2:
        continue

    # ROC Curve
    fpr, tpr, _ = roc_curve(label_test[indices], y_pred_proba[indices])
    roc_auc = roc_auc_score(label_test[indices], y_pred_proba[indices])
    aucs.append(roc_auc)
    tprs.append(np.interp(mean_fpr, fpr, tpr))
    tprs[-1][0] = 0.0

    # PRC Curve
    precision, recall, _ = precision_recall_curve(label_test[indices], y_pred_proba[indices])
    auc_prc = average_precision_score(label_test[indices], y_pred_proba[indices])
    aucs_prc.append(auc_prc)
    precisions.append(np.interp(mean_recall, recall[::-1], precision[::-1]))

# Calculate mean and std for ROC
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = np.mean(aucs)
std_auc = np.std(aucs)
tprs_upper = np.minimum(mean_tpr + np.std(tprs, axis=0), 1)
tprs_lower = np.maximum(mean_tpr - np.std(tprs, axis=0), 0)

# Calculate mean and std for PRC
mean_precision = np.mean(precisions, axis=0)
mean_auc_prc = np.mean(aucs_prc)
std_auc_prc = np.std(aucs_prc)
precisions_upper = np.minimum(mean_precision + np.std(precisions, axis=0), 1)
precisions_lower = np.maximum(mean_precision - np.std(precisions, axis=0), 0)


#%% 95% CI
# Sort the bootstrap AUCs
sorted_aucs = np.sort(aucs)

# Get the 2.5th and 97.5th percentiles
lower = np.percentile(sorted_aucs, 2.5)
upper = np.percentile(sorted_aucs, 97.5)

print(f"95% Confidence Interval for ROC AUC: [{lower:.3f}, {upper:.3f}]")

sorted_aucs_prc = np.sort(aucs_prc)
lower_prc = np.percentile(sorted_aucs_prc, 2.5)
upper_prc = np.percentile(sorted_aucs_prc, 97.5)

print(f"95% Confidence Interval for PRC AUC: [{lower_prc:.3f}, {upper_prc:.3f}]")

#%% plot using 95% CI AUROC
import matplotlib.pyplot as plt

# Convert tprs to numpy array for easy slicing
tprs_np = np.array(tprs)

# Compute 95% confidence interval bounds at each mean_fpr point
tprs_lower = np.percentile(tprs_np, 2.5, axis=0)
tprs_upper = np.percentile(tprs_np, 97.5, axis=0)

# Compute 95% CI for AUC
sorted_aucs = np.sort(aucs)
ci_lower = np.percentile(sorted_aucs, 2.5)
ci_upper = np.percentile(sorted_aucs, 97.5)

# Plot the ROC curve
plt.figure(figsize=(3,2),dpi=600) # these 2 separate of 3,2 size
plt.plot(mean_fpr, mean_tpr, color='#377eb8', lw=1, label=f'AUROC = {mean_auc:.3f}')
plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='red', alpha=0.2, label="95% CI")
plt.plot([0, 1], [0, 1], color='black', lw=1, linestyle='--', label='Random classifier')
plt.xlabel('False Positive Rate (FPR)', fontsize=8)
plt.ylabel('True Positive Rate (TPR)', fontsize=8)
plt.legend(loc='lower right', fontsize=7)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


# Precision-Recall curve

# Compute 95% CI for PR AUC
sorted_aucs_prc = np.sort(aucs_prc)
ci_lower_prc = np.percentile(sorted_aucs_prc, 2.5)
ci_upper_prc = np.percentile(sorted_aucs_prc, 97.5)

precisions_np = np.array(precisions)
precisions_lower = np.percentile(precisions_np, 2.5, axis=0)
precisions_upper = np.percentile(precisions_np, 97.5, axis=0)

# Plot PRC curve
plt.figure(figsize=(3,2),dpi=600)
plt.plot(mean_recall, mean_precision, color='#377eb8', lw=1, label=f'AUPRC = {mean_auc_prc:.3f}')
plt.fill_between(mean_recall, precisions_lower, precisions_upper, color='red', alpha=0.2, label='95% CI')
positive_ratio = sum(label_test) / len(label_test)
plt.plot([0, 1], [positive_ratio, positive_ratio], color='black', lw=1, linestyle='--', label='Random classifier')
plt.xlabel('Recall', fontsize=8)
plt.ylabel('Precision', fontsize=8)
plt.xlim([-0.05, 1.05])
plt.ylim([-0.05, 1.05])
plt.legend(loc='lower right', fontsize=7)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()


#%% plot iFNBoost performance error bars

plt.figure(figsize=(3, 2), dpi=300)
metrics = ['Sensitivity', 'Specificity', 'Accuracy', 'F1 score', 'MCC']
means = [0.8152, 0.8113, 0.8130, 0.7933, 0.6235]
ci_lower = [0.8032, 0.7989, 0.8041, 0.7829, 0.6061]
ci_upper = [0.8278, 0.8236, 0.8220, 0.8038, 0.6421]

# Compute error bars
errors = [[m - l for m, l in zip(means, ci_lower)], [u - m for u, m in zip(ci_upper, means)]]

bars = plt.barh(metrics, means, xerr=errors, capsize=5, color='#115316', alpha=0.7)  #377eb8 blue

# Add values inside bars with 95% CI
for bar, mean, lower, upper in zip(bars, means, ci_lower, ci_upper):
    text = f"{mean:.3f} ({lower:.3f} - {upper:.3f})"
    plt.text(bar.get_width() / 2,  # Position at center of bar
             bar.get_y() + bar.get_height() / 2,  
             text, ha='center', va='center', fontsize=8, color='white')



plt.xlabel("Score")
plt.xlim(0, 1)
# plt.title("Model Performance with 95% CI")
plt.gca().invert_yaxis()
# plt.grid(axis='x', linestyle='--', alpha=0.6)

plt.show()

#%% plot species and host wise

categories = ['Human', 'Mouse']# 'M.m. C57BL/6', 'Roseolovirus', 'Dengue'] #840 trees
#this mouse is M.m C57BL/6

sens = [0.8465, 0.7181]#, 0.9576, 0.6911,0.7714, 0.81]
spf = [0.666, 0.9755]#, 0.7197, 0.7339, 0.724, 0.834]
accuracy = [0.766, 0.9192]#, 0.6826, 0.7151, 0.7469, 0.8232]
f1 = [0.801, 0.7955]#, 0.6595, 0.679, 0.7466, 0.8023]
mcc = [0.524, 0.7527]#, 0.3632, 0.42366, 0.4954, 0.6425]

# Create a matrix with MCC and Accuracy
data = np.array([sens, spf,accuracy, f1, mcc])

# Create a DataFrame
import pandas as pd
dfm = pd.DataFrame(data, columns=categories, index=['Sensitivity', 'Specificity','Accuracy', 'F1 score', 'MCC'])

dfmt = dfm.transpose()
# Create a heatmap using seaborn
plt.figure(figsize=(3, 3), dpi=600)
# colormap = sns.color_palette('#377eb8')
colormap = sns.light_palette('#115316')

heatmap = sns.heatmap(dfm, annot=True, cmap=colormap, fmt=".3f", linewidths=1,annot_kws={"size": 8}, square=True,vmin=0, vmax=1)

# heatmap.xaxis.set_visible(False) # removes the x-ticks completely

heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=0, fontsize=7)
# heatmap.set_xticklabels([])
heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0, fontsize=7)

plt.xlabel('Host', fontsize=8)
# Display the plot
plt.tight_layout()
plt.show()




