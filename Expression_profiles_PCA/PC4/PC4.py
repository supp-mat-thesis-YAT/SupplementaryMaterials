# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 22:10:53 2019

@author: Tamal
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from sklearn.decomposition import PCA
plt.rc('font', size = 8)
sns.set_palette(sns.color_palette('Paired'))
from sklearn.cluster import KMeans
sns.set()
sns.set_style("whitegrid", {'axes.grid' : False})

data = pd.read_csv('chRNAseq_data.csv', index_col = 0)
data = data.T
data.drop('L3T2_1', axis = 0, inplace = True)


def cutoff_filter(mat):
    f = mat.loc[:,((mat == 0).any(axis = 0) >= 1)]
    
    return f.max(axis = 0) >= 20

data = data.loc[:,(data != 0).all(axis = 0) | cutoff_filter(data)]

dataset = data.copy()

sample_list = ['L1T1_1', 'L1T1_2', 'L1T1_3', 'L1T2_1', 'L1T2_2', 'L1T2_3', 'L1T3_1', 'L1T3_2', 'L1T3_3', 'L1T4_1', 'L1T4_2', 'L1T4_3', 'L3T2_1', 'L3T2_2', 'L3T3_1', 'L3T3_2', 'L3T3_3', 'L3T4_1', 'L3T4_2', 'L3T4_3', 'L5T3_1', 'L5T3_2', 'L5T3_3', 'L5T4_1', 'L5T4_2', 'L5T4_3', 'L7T4_1', 'L7T4_2', 'L7T4_3']

y_true = sorted(list(range(10)) *3)
del y_true[12]

X = dataset.values 

sc = StandardScaler()
sc_X = sc.fit_transform(X)

pca = PCA(n_components = 4)
X_train = pca.fit_transform(sc_X)

pcs = ['PC1', 'PC2', 'PC3', 'PC4']
loadings = pd.DataFrame(pca.components_,columns = data.columns, index = pcs)
axes_pos = np.arange(len(pcs))

def variable_contrib(index_n):
    pc_loadings = loadings.loc[loadings.index[index_n]].to_frame().abs()
    sorted_loadings = pc_loadings.sort_values(pcs[index_n], ascending = False)
    gene_contrib = sorted_loadings.loc[(sorted_loadings[pcs[index_n]] >= 0.01)].T
    genes = gene_contrib.columns.values.tolist()
    
    return gene_contrib, genes


genes_pc1 = variable_contrib(0)[1][:75]
genes_pc2 = variable_contrib(1)[1][:75]
genes_pc3 = variable_contrib(2)[1][:75]
genes_pc4 = variable_contrib(3)[1][:75]



sample_group = ['L1T1', 'L1T2','L1T3', 'L1T4', 'L3T2', 'L3T3', 'L3T4', 'L5T3', 'L5T4', 'L7T4']
group_position = np.arange(len(sample_group))

#by genes
fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[0]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])
plt.tight_layout()

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[1]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])
plt.tight_layout()

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[2]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])
plt.tight_layout()

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[3]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])
plt.tight_layout()

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[4]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])
plt.tight_layout()

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[5]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])
plt.tight_layout()

plt.savefig('PC4_expression_1_6_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[6]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[7]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[8]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[9]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[10]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[11]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_7_12_PCA.png')
plt.show()


fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[12]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[13]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[14]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[15]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[16]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[17]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_13_18_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[18]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[19]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[20]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[21]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[22]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[23]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_19_24_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[24]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[25]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[26]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[27]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[28]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[29]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_25_30_PCA.png')
plt.show()


fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[30]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[31]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[32]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[33]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[34]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[35]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_31_36_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[36]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[37]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[38]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[39]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[40]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[41]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_37_42_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[42]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[43]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[44]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[45]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[46]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[47]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_43_48_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[48]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[49]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[50]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[51]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[52]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[53]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_49_54_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[54]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[55]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[56]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[57]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[58]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[59]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_55_60_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[60]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[61]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[62]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[63]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[64]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[65]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_61_66_PCA.png')
plt.show()

fig = plt.figure(figsize = (16, 12))

plt.subplot(2, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[66]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(2, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[67]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(2, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[68]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.subplot(2, 3, 4)
plt.scatter(y_true,data.loc[:,genes_pc4[69]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[3])

plt.subplot(2, 3, 5)
plt.scatter(y_true,data.loc[:,genes_pc4[70]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[4])

plt.subplot(2, 3, 6)
plt.scatter(y_true,data.loc[:,genes_pc4[71]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[5])

plt.savefig('PC4_expression_67_72_PCA.png')
plt.show()


fig = plt.figure(figsize = (16, 6))

plt.subplot(1, 3, 1)
plt.scatter(y_true,data.loc[:,genes_pc4[72]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[0])

plt.subplot(1, 3, 2)
plt.scatter(y_true,data.loc[:,genes_pc4[73]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[1])

plt.subplot(1, 3, 3)
plt.scatter(y_true,data.loc[:,genes_pc4[74]],  edgecolors = 'y', s = 50, alpha = 0.7)
plt.xticks(group_position, sample_group, rotation = 'vertical', fontsize = 10)
plt.ylabel(genes_pc4[2])

plt.savefig('PC4_expression_73_75_PCA.png')
plt.show()
