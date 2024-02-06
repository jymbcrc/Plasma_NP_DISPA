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

plasma_pep


# In[4]:


def make_lists_pro_pep(filenames=plasma_pro):
    '''put all proteins or peptides from each NP to a lists in list'''    
    dataframes = []

    for name in filenames:
        df = pd.read_csv(name)
        df.set_index('Unnamed: 0', inplace=True)
        dataframes.append(df)

    pro_dd, pro_np1, pro_np2, pro_np3, pro_np4, pro_np5 = dataframes
    list_dd=pro_dd.index.tolist()
    list_np1=pro_np1.index.tolist()
    list_np2=pro_np2.index.tolist()
    list_np3=pro_np3.index.tolist()
    list_np4=pro_np4.index.tolist()
    list_np5=pro_np5.index.tolist()
    toge_all=[list_dd,list_np1,list_np2,list_np3,list_np4,list_np5]
    return(toge_all)


# In[5]:


toge_pro = make_lists_pro_pep(filenames=plasma_pro)
toge_pep = make_lists_pro_pep(filenames=plasma_pep)


# In[6]:


def cal_any_ratios(values=toge_pro,n=2):
    '''calculate the length ratio of any n lists in a list to the original list'''
    from itertools import combinations    
    # get all values from the input list and remove duplicates
    total = remove_dul(tar_list=[item for sublist in values for item in sublist])
    results = []

    # Iterate over each combination of n values in the list
    for group in combinations(values, n):
        unique_group_values = list(set(item for sublist in group for item in sublist))
        #unique_group_values =remove_dul(tar_list=list(item for sublist in group for item in sublist))
        ratio = len(unique_group_values) / len(total)
        indices = [values.index(sublist) for sublist in group]
        results.append((indices, ratio))
    # Convert the results to a DataFrame
    df = pd.DataFrame(results, columns=['Group', 'Ratio'])
    df.set_index('Group',inplace=True)
    return(df)  


 


# In[7]:


df_bio=pd.read_csv(r'F:\\Seer_project\\R_data\\chord_diagram_biomarkers2.csv')
df_bio.set_index('Item',inplace=True)


# In[8]:



bio_dd=df_bio[df_bio['DD']==1].index.tolist()
bio_np1=df_bio[df_bio['NP1']==1].index.tolist()
bio_np2=df_bio[df_bio['NP2']==1].index.tolist()
bio_np3=df_bio[df_bio['NP3']==1].index.tolist()
bio_np4=df_bio[df_bio['NP4']==1].index.tolist()
bio_np5=df_bio[df_bio['NP5']==1].index.tolist()

toge_bio=[bio_dd,bio_np1,bio_np2,bio_np3,bio_np4,bio_np5]


# In[9]:


Intensity=[17.442438, 17.112451, 17.010393, 17.677678, 17.361636, 17.457459]
Quan=[0.4, 0.62, 0.64, 0.44, 0.85, 0.71]
CV=[0.06729475100942127,0.05672149744753262,0.06053268765133172,0.07251631617113852,0.04681647940074907,0.04819277108433735]


# In[11]:


def cal_list_means(values=Intensity,n=2):
    '''calculate the mean value of any n values in a list'''   
    from itertools import combinations
    results = []
    # Iterate over each combination of n values in the list
    for group in combinations(values, n):
        sumvalue = sum(item for item in group)
        means = sumvalue / len(group)
        indices = [values.index(item) for item in group]
        results.append((indices, means))
    # Convert the results to a DataFrame
    df = pd.DataFrame(results, columns=['Group', 'Ratio'])
    df.set_index('Group',inplace=True)
    return(df) 


# In[70]:


i=2

dfbio_2=cal_any_ratios(values=toge_bio,n=i).rename({'Ratio': 'Biomarker'}, axis=1)
dfpep_2=cal_any_ratios(values=toge_pep,n=i).rename({'Ratio': 'Peptide'}, axis=1)
dfpro_2=cal_any_ratios(values=toge_pro,n=i).rename({'Ratio': 'Protein'}, axis=1)
dfinten_2=cal_list_means(values=Intensity,n=i).rename({'Ratio': 'Intensity'}, axis=1)
dfquan_2=cal_list_means(values=Quan,n=i).rename({'Ratio': 'Quan'}, axis=1)
dfcv_2=cal_list_means(values=CV,n=i).rename({'Ratio': 'CV'}, axis=1)


# In[71]:


# merge dataframes

all_toge2 = pd.concat([dfpro_2,dfbio_2,dfpep_2, dfinten_2,dfquan_2,dfcv_2], axis=1)


# In[72]:


all_toge2


# In[73]:


# apply standard scaler to the ratar dataframe 
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

#scaler = StandardScaler()
scaler = MinMaxScaler()

scaled_data = scaler.fit_transform(all_toge2)

# Convert scaled data back to DataFrame for better visualization and to retain column names
normalized = pd.DataFrame(scaled_data, columns=all_toge2.columns,index=all_toge2.index)
normalized


# In[74]:


# calculate the final score of any combination

A,B,C,D,E,F=0.15, 0.2, 0.15, 0.2, 0.1, 0.2    # define the weight
df = normalized
scores=[]
for i in range(len(df)):
    score=df.CV[i]*A+df.Biomarker[i]*B+df.Peptide[i]*C+df.Protein[i]*D+df.Intensity[i]*E+df.Quan[i]*F
    scores.append(score)
scores


# In[ ]:





# In[75]:



import matplotlib as mpl

# Set global font
mpl.rcParams['font.family'] = 'Arial'

plt.rcParams['figure.figsize']=(12,7)

# Define the data
x = [str(item) for item in df.index]
y = scores

# Create the plot
plt.plot(x, y, marker='o')
# Rotate x-axis labels
plt.xticks(fontsize=22, rotation=90)
plt.yticks(fontsize=22)
#plt.xlabel('X-axis', fontsize=26)
plt.ylabel('Score',fontsize=24)
#plt.title('Simple Line Plot')



#plt.savefig('F:\\Seer_project\\figures\\combination analysis\\4nps.svg', dpi=800,bbox_inches='tight')

plt.show()


# In[77]:


# Create bubble plot
z=np.array(y)*600

plt.figure(figsize=(12, 7))
bubble = plt.scatter(x, y, s=z,color='#701EB3')

# Adding titles and labels
plt.xticks(fontsize=22, rotation=90)
plt.yticks(fontsize=22)
#plt.xlabel('X-axis', fontsize=26)
plt.ylabel('Score',fontsize=24)
plt.ylim(0.1,0.9)
plt.savefig('F:\\Seer_project\\figures\\combination analysis\\2nps_bubble.svg', dpi=800,bbox_inches='tight')
# Show plot
plt.grid(True, alpha=0.6)
plt.show()


# In[66]:


df.index=['DD','NP1','NP2','NP3','NP4','NP5']


# In[262]:


df


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[263]:


df_radar=pd.read_excel(r'F:\\Seer_project\\test.xlsx',sheet_name='Sheet2')
df_radar.set_index('Unnamed: 0',inplace=True)


# In[264]:


df_radar.iloc[:,0].tolist()


# In[265]:


# apply standard scaler to the ratar dataframe 
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler

#scaler = StandardScaler()
scaler = MinMaxScaler()

scaled_data = scaler.fit_transform(df_radar)

# Convert scaled data back to DataFrame for better visualization and to retain column names
normalized_df_radar = pd.DataFrame(scaled_data, columns=df_radar.columns,index=df_radar.index)
normalized_df_radar


# In[267]:


normalized_df_radar.to_csv(r'F:\\Seer_project\\R_data\\radar3.csv')


# In[168]:


df=normalized_df_radar
df = df.rename({'1/CV': 'CV'}, axis=1)


# In[159]:


df.Intensity


# In[160]:


A,B,C,D,E,F=0.15, 0.2, 0.15, 0.2, 0.1, 0.2

scores=[]
for i in range(6):
    score=df.CV[i]*A+df.Biomarker[i]*B+df.Peptide[i]*C+df.Protein[i]*D+df.Intensity[i]*E+df.Quan[i]*F
    scores.append(score)
scores


# In[ ]:





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





# In[ ]:





# In[ ]:





# In[ ]:





# In[60]:


def get_zip(df=pro_dd):
    '''get the index ("protein name") and ion counts value into a new list from csodiaq outputs' dataframe'''
    newlist = []
    # Using zip to iterate through two lists simultaneously
    for protein, ionCount in zip(df.index.tolist(), df['mean'].tolist()):
        newlist.append([protein, ionCount])
    return(newlist)


# In[61]:


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


# In[62]:


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


# In[67]:





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


# In[19]:


#len([item for item in biomarker if item in allpro])


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




