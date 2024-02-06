#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pyteomics
import lxml
from pyteomics import mzml, mass
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib import pyplot
import matplotlib
import seaborn as sns
import re


# In[2]:


sample=[0,2,3,4]
def remove_dul(tar_list=sample):
    res=[]
    [res.append(x) for x in tar_list if x not in res]
    return(res)


# In[3]:


#  import the quantification files of all NPs and DD

dir = "F:\\Seer_project"

plasma_pro = []
plasma_pep = []

for file in os.listdir(dir):
    if 'common_proteins' in file:
        plasma_pro.append(os.path.join(dir, file))
    if 'common_peptides' in file:
        plasma_pep.append(os.path.join(dir, file))

plasma_pro


# In[4]:


filename=plasma_pro

pro_dd=pd.read_csv(filename[0])
pro_dd.set_index('Unnamed: 0',inplace=True)
pro_np1=pd.read_csv(filename[1])
pro_np1.set_index('Unnamed: 0',inplace=True)
pro_np2=pd.read_csv(filename[2])
pro_np2.set_index('Unnamed: 0',inplace=True)
pro_np3=pd.read_csv(filename[3])
pro_np3.set_index('Unnamed: 0',inplace=True)
pro_np4=pd.read_csv(filename[4])
pro_np4.set_index('Unnamed: 0',inplace=True)
pro_np5=pd.read_csv(filename[5])
pro_np5.set_index('Unnamed: 0',inplace=True)


# In[ ]:





# In[5]:


all_proteins=pro_dd.index.tolist()+pro_np1.index.tolist()+pro_np2.index.tolist()+pro_np3.index.tolist()+pro_np4.index.tolist()+pro_np5.index.tolist()
np_proteins=pro_np1.index.tolist()+pro_np2.index.tolist()+pro_np3.index.tolist()+pro_np4.index.tolist()+pro_np5.index.tolist()


# In[6]:


list_all=remove_dul(tar_list=all_proteins)
list_dd=pro_dd.index.tolist()
list_np1=pro_np1.index.tolist()
list_np2=pro_np2.index.tolist()
list_np3=pro_np3.index.tolist()
list_np4=pro_np4.index.tolist()
list_np5=pro_np5.index.tolist()


# In[11]:


len(list_all)


# In[13]:


len(remove_dul(tar_list=np_proteins))


# In[37]:


toge_all=[list_dd,list_np1,list_np2,list_np3,list_np4,list_np5]


# In[38]:


wo=[]
for item in toge_all:
    wo.append(len(item)/len(list_all))
wo


# In[39]:


df_bio=pd.read_csv(r'F:\\Seer_project\\R_data\\chord_diagram_biomarkers2.csv')


# In[40]:


df_bio


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[12]:


def get_pro_names(protein = pro_np1):
    sample  = []
    for i in protein.index.tolist():
        sample.append(i.split('|')[1])
    return(sample)


# In[13]:


len(get_pro_names(protein = pro_np2))


# In[ ]:





# In[7]:


def get_zip(df=pro_dd):
    '''get the index ("protein name") and ion counts value into a new list from csodiaq outputs' dataframe'''
    newlist = []
    # Using zip to iterate through two lists simultaneously
    for protein, ionCount in zip(df.index.tolist(), df['mean'].tolist()):
        newlist.append([protein, ionCount])
    return(newlist)


# In[9]:


list1 = get_zip(df=pro_dd)
list2 = get_zip(df=pro_np1)
list3 = get_zip(df=pro_np2)
list4 = get_zip(df=pro_np3)
list5 = get_zip(df=pro_np4)
list6 = get_zip(df=pro_np5)

# Convert lists to dictionaries
dict1 = dict(list1)
dict2 = dict(list2)
dict3 = dict(list3)
dict4 = dict(list4)
dict5 = dict(list5)
dict6 = dict(list6)

# Fetch unique keys (indices) from both dictionaries
all_keys = set(dict1.keys()) | set(dict2.keys()) | set(dict3.keys()) | set(dict4.keys()) | set(dict5.keys()) | set(dict6.keys())

# Merge lists based on keys
merged_list = []
for key in sorted(all_keys):
    value1 = dict1.get(key, float('nan'))
    value2 = dict2.get(key, float('nan'))
    value3 = dict3.get(key, float('nan'))
    value4 = dict4.get(key, float('nan'))
    value5 = dict5.get(key, float('nan'))
    value6 = dict6.get(key, float('nan'))
    merged_list.append([key, value1, value2, value3, value4, value5, value6])

print(merged_list) 


# In[10]:


data_list=merged_list
# Splitting indices and data
indices = [item[0] for item in data_list]
data = [item[1:] for item in data_list]

# Creating DataFrame
newdf = pd.DataFrame(data, index=indices, columns=['DD', 'NP1','NP2','NP3','NP4','NP5'])

newdf


# In[68]:


uk_log2=np.log2(uk)


# In[69]:


uk_log2.mean()


# In[63]:


uk=newdf
#um=uk.transpose()
uk_log10=np.log10(uk)
uk_log10


# In[ ]:





# In[12]:


data=uk_log10.fillna(0)


# In[13]:


newindex=[item.split('|')[-1].split('_')[0] for item in uk_log10.index]


# In[14]:


pro_en=[item.split('|')[1] for item in uk_log10.index]


# In[15]:



uk_log10['pro_names']=pro_en
uk_log10['new_index']=newindex


# In[16]:


uk_log10


# In[17]:


df_biomarker=pd.read_excel(r'F:\\Seer_project\\FDA_biomarkers.xlsx',sheet_name='Sheet3')


# In[18]:


biomarker=df_biomarker['biomarkers'].tolist()


# In[20]:


uk_log10['Biomarkers'] = np.where(uk_log10['pro_names'].isin(biomarker), 'yes', 'no')


# In[21]:


uk_log10


# In[22]:


#newdf.to_excel(r'F:\\Seer_project\\NPs_DD_Highcon_only.xlsx')


# In[23]:


data.index=newindex


# In[24]:


sav_df=uk_log10.fillna(1)


# In[25]:


sav_df


# In[26]:


#uk_log10.set_index('new_index',inplace=True)


# In[27]:


#uk_log10.to_excel(r'F:\\Seer_project\\NPs_DD_Highcon_only.xlsx')


# In[28]:


#data=uk_log10.drop(columns=['pro_names','new_index'])


# In[29]:


uk_log10


# In[30]:


data['Biomarkers']=uk_log10['Biomarkers'].tolist()
data


# In[31]:


dfdf=data.iloc[:,:6]
dfdf


# In[32]:


data_transposed=dfdf.transpose()
data_transposed


# In[52]:



matplotlib.rcParams['font.family'] = "Arial"
# Map the species to colors
color_dict = {"yes": 'red', "no": '#DFA828'}

col_colors = data['Biomarkers'].map(color_dict).tolist()

# Transpose the  data to have species as columns (or x-axis observations)
data_transposed = data.drop('Biomarkers', axis=1).transpose()

# we can map the corresponding species colors to these columns
colors_transposed = dict(zip(data.index, col_colors))

# Display the transposed clustermap with colored labels on the x-axis
g=sns.clustermap(data_transposed, col_colors=data.index.map(colors_transposed), cmap='viridis',figsize=(18.7, 7), fmt=".3f",
                cbar_pos=(0.16, 0.76, 0.02, 0.20) )   #(colorbar size and position)
#g=sns.clustermap(data_transposed, cmap='viridis',figsize=(20, 7))

# 调整行树状图的线条粗细
for l in g.ax_row_dendrogram.collections:
    l.set_linewidth(1)  # 设置线宽为2.5

# 调整列树状图的线条粗细
for l in g.ax_col_dendrogram.collections:
    l.set_linewidth(1)  # 设置线宽为2.5
# Set tick label size
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=12)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=17)
#g.ax_heatmap.set_yticklabels([])
g.ax_row_dendrogram.set_position([0.126, 0.12, 0.20, 0.65])  # [x, y, width, height]
cbar = g.cax
for label in cbar.get_yticklabels():
    label.set_size(19)  # 设置字体大小为15

# 调整列树形图的大小
#g.ax_col_dendrogram.set_position([0.3, 0.9, 0.6, 0.3])  # [x, y, width, height]
fig_path= 'F:\\Seer_project\\figures\\' 
figname='heatmap_targeted_highcon'
plt.savefig(fig_path + "%s.svg" % figname,dpi=800, bbox_inches='tight')  
plt.show()


# # plot heatmap of all proteins identified by scouting experiments

# In[96]:


df_scout=pd.read_excel(r'F:\\Seer_project\\NPs_scout.xlsx',sheet_name="Sheet1")


# In[97]:


#df_scout.set_index('Unnamed: 0',inplace=True)


# In[98]:


df_scout.index=df_scout['new_index'].tolist()
df_scout


# In[99]:


data2=df_scout.drop(columns=['pro_names','new_index','Unnamed: 0'])
data2=data2.fillna(0)
data2


# In[100]:


dfdf2=data2.iloc[:,:6]


# In[126]:


# Map the species to colors
color_dict2 = {"yes": 'red', "no": '#3695C7'}

col_colors2 = data2['Biomarkers'].map(color_dict2).tolist()

# Transpose the  data to have species as columns (or x-axis observations)
data_transposed2 = data2.drop('Biomarkers', axis=1).transpose()

# we can map the corresponding species colors to these columns
colors_transposed2 = dict(zip(data2.index, col_colors2))

# Display the transposed clustermap with colored labels on the x-axis
g=sns.clustermap(data_transposed2, col_colors=data2.index.map(colors_transposed2), cmap='viridis',figsize=(19, 7), fmt=".3f",
                cbar_pos=(0.156, 0.76, 0.02, 0.20) )   #(colorbar size and position)
#g=sns.clustermap(data_transposed, cmap='viridis',figsize=(20, 7))

# 调整行树状图的线条粗细
for l in g.ax_row_dendrogram.collections:
    l.set_linewidth(1)  # 设置线宽为2.5

# 调整列树状图的线条粗细
for l in g.ax_col_dendrogram.collections:
    l.set_linewidth(1)  # 设置线宽为2.5
# Set tick label size
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=13)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=18)
#g.ax_heatmap.set_yticklabels([])
g.ax_row_dendrogram.set_position([0.14, 0.12, 0.20, 0.65])  # [x, y, width, height]
cbar = g.cax
for label in cbar.get_yticklabels():
    label.set_size(18)  # 设置字体大小为15


# 调整列树形图的大小
#g.ax_col_dendrogram.set_position([0.3, 0.71, 0.6, 0.3])  # [x, y, width, height]
fig_path= 'F:\\Seer_project\\figures\\' 
figname='heatmap_nontargeted'
plt.savefig(fig_path + "%s.svg" % figname,dpi=800, bbox_inches='tight')  
plt.show()


# In[ ]:




