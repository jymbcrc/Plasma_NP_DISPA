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
import math


# In[2]:


# plot a violin figure for all NPs and neat plasma


# In[3]:


sample=[2,4,6,8,3]

def remove_dul(tar_list=sample):
    res=[]
    [res.append(x) for x in tar_list if x not in res]
    return(res)

def average(sample):
    return sum(sample) / len(sample)

def procount(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    dic={}
    for i in range(len(plasma_pro)):
        df=pd.read_csv(plasma_pro[i]) 
        yummy=len(remove_dul(tar_list=df['protein'].tolist()))
        dic[str(plasma_pro[i].split('\\')[-1])]=yummy
        dfnew=pd.DataFrame.from_dict(dic, orient='index')
    return(dfnew)
def pepcount(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    dic={}
    for i in range(len(plasma_pro)):
        df=pd.read_csv(plasma_pro[i]) 
        yummy=len(remove_dul(tar_list=df['peptide'].tolist()))
        dic[str(plasma_pro[i].split('\\')[-1])]=yummy
        dfnew=pd.DataFrame.from_dict(dic, orient='index')
    return(dfnew)

def plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq"):
    ['''take the common_peptides and common_protein files out''']
    plasma_pro = []
    plasma_pep = []
    for file in os.listdir(dir):
        if 'common_proteins' in file:
            plasma_pro.append(os.path.join(dir, file))
        if 'common_peptides' in file:
            plasma_pep.append(os.path.join(dir, file))
    dfpep=pd.read_csv(plasma_pep[0])
    dfpro=pd.read_csv(plasma_pro[0])
    dfpep.set_index('Unnamed: 0',inplace=True)
    dfpro.set_index('Unnamed: 0',inplace=True)
    toge=[dfpep,dfpro]
    return(toge)


# In[4]:


plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]    


# In[5]:


i = 0    # i=1 extracts protein file, i=0 extracts peptide file

dfpro_np5=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[i]
dfpro_np4=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[i]
dfpro_np3=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2")[i]
dfpro_np2=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[i]
dfpro_np1=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[i]
dfpro_dd=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[i]


# In[6]:


dfpro_np5


# In[7]:


start_index = 0
end_index = 3
pronp5 = dfpro_np5.iloc[:, start_index:end_index].T
pronp4 = dfpro_np4.iloc[:, start_index:end_index].T
pronp3 = dfpro_np3.iloc[:, start_index:end_index].T
pronp2 = dfpro_np2.iloc[:, start_index:end_index].T
pronp1 = dfpro_np1.iloc[:, start_index:end_index].T
prodd = dfpro_dd.iloc[:, start_index:end_index].T
start_index = 3
end_index = 6
pronp5_2 = dfpro_np5.iloc[:, start_index:end_index].T
pronp4_2 = dfpro_np4.iloc[:, start_index:end_index].T
pronp3_2 = dfpro_np3.iloc[:, start_index:end_index].T
pronp2_2 = dfpro_np2.iloc[:, start_index:end_index].T
pronp1_2 = dfpro_np1.iloc[:, start_index:end_index].T
prodd_2 = dfpro_dd.iloc[:, start_index:end_index].T
start_index = 6
end_index = 9
pronp5_3 = dfpro_np5.iloc[:, start_index:end_index].T
pronp4_3 = dfpro_np4.iloc[:, start_index:end_index].T
pronp3_3 = dfpro_np3.iloc[:, start_index:end_index].T
pronp2_3 = dfpro_np2.iloc[:, start_index:end_index].T
pronp1_3 = dfpro_np1.iloc[:, start_index:end_index].T
prodd_3 = dfpro_dd.iloc[:, start_index:end_index].T
start_index = 9
end_index = 12
pronp5_4 = dfpro_np5.iloc[:, start_index:end_index].T
pronp4_4 = dfpro_np4.iloc[:, start_index:end_index].T
pronp3_4 = dfpro_np3.iloc[:, start_index:end_index].T
pronp2_4 = dfpro_np2.iloc[:, start_index:end_index].T
pronp1_4 = dfpro_np1.iloc[:, start_index:end_index].T
prodd_4 = dfpro_dd.iloc[:, start_index:end_index].T


# In[8]:


cv = lambda x: np.std(x, ddof=1) / np.mean(x) * 100 


# In[9]:


newdf={'DD':prodd.apply(cv).tolist(),'DD_2':prodd_2 .apply(cv).tolist(),'DD_3':prodd_3 .apply(cv).tolist(),'DD_4':prodd_4.apply(cv).tolist(),
    'NP1':pronp1.apply(cv).tolist(),'NP1_2':pronp1_2 .apply(cv).tolist(),'NP1_3':pronp1_3 .apply(cv).tolist(),'NP1_4':pronp1_4.apply(cv).tolist(),
    'NP2':pronp2.apply(cv).tolist(),'NP2_2':pronp2_2 .apply(cv).tolist(),'NP2_3':pronp2_3.apply(cv).tolist(),'NP2_4':pronp2_4.apply(cv).tolist(),
    'NP3':pronp3.apply(cv).tolist(),'NP3_2':pronp3_2 .apply(cv).tolist(),'NP3_3':pronp3_3.apply(cv).tolist(),'NP3_4':pronp3_4.apply(cv).tolist(),
    'NP4':pronp4.apply(cv).tolist(),'NP4_2':pronp4_2 .apply(cv).tolist(),'NP4_3':pronp4_3.apply(cv).tolist(),'NP4_4':pronp4_4.apply(cv).tolist(),
    'NP5': pronp5.apply(cv).tolist(),'NP5_2': pronp5_2 .apply(cv).tolist(),'NP5_3': pronp5_3.apply(cv).tolist(),'NP5_4': pronp5_4.apply(cv).tolist()}


# In[10]:


uk=pd.DataFrame.from_dict(newdf, orient='index')


# In[11]:


um=uk.transpose()


# In[12]:


um


# In[13]:


xlabels=['DD(0.4ug/ul)','DD(0.2ug/ul)','DD(0.1ug/ul)','DD(0.05ug/ul)',
         'NP1(0.4ug/ul)','NP1(0.2ug/ul)','NP1(0.1ug/ul)','NP1(0.05ug/ul)',
         'NP2(0.4ug/ul)','NP2(0.2ug/ul)','NP2(0.1ug/ul)','NP2(0.05ug/ul)',
         'NP3(0.4ug/ul)','NP3(0.2ug/ul)','NP3(0.1ug/ul)','NP3(0.05ug/ul)',
         'NP4(0.4ug/ul)','NP4(0.2ug/ul)','NP4(0.1ug/ul)','NP4(0.05ug/ul)',
         'NP5(0.4ug/ul)','NP5(0.2ug/ul)','NP5(0.1ug/ul)','NP5(0.05ug/ul)']


# In[14]:


xlabels2=['DD(0.4)','DD(0.2)','DD(0.1)','DD(0.05)',
         'NP1(0.4)','NP1(0.2)','NP1(0.1)','NP1(0.05)',
         'NP2(0.4)','NP2(0.2)','NP2(0.1)','NP2(0.05)',
         'NP3(0.4)','NP3(0.2)','NP3(0.1)','NP3(0.05)',
         'NP4(0.4)','NP4(0.2)','NP4(0.1)','NP4(0.05)',
         'NP5(0.4)','NP5(0.2)','NP5(0.1)','NP5(0.05)']


# In[15]:


um.columns=xlabels2


# In[16]:


um


# In[23]:


sum(list(um.mean()))/len(list(um.mean()))


# In[31]:


i=20

np.mean(list(um.mean())[i:i+4])


# In[ ]:





# In[143]:


from matplotlib import rcParams
custom_palette = ['#8B5C96','#8B5C96','#8B5C96','#8B5C96',
                  '#69A959','#69A959','#69A959','#69A959',
                  '#4B7FAB','#4B7FAB','#4B7FAB','#4B7FAB',
                  '#E4B58D','#E4B58D','#E4B58D','#E4B58D',
                  '#ef6f6c','#ef6f6c','#ef6f6c','#ef6f6c',
                  '#465775','#465775','#465775','#465775']


# figure size in inches
plt.rcParams['figure.figsize']=(18,7)

sns.set(font_scale=2)
sns.set_style("whitegrid")
 
sns.violinplot(data=um,palette=custom_palette)
#sns.stripplot(data=um, jitter=True,alpha=0.6,size=4, color="black")
#plt.ylabel(' LFQ CV(protein, %)', fontsize=26)
plt.ylabel(' LFQ CV(peptide, %)', fontsize=26)

plt.xticks(fontsize=20,rotation=-60)
#plt.xticks([])
plt.yticks(fontsize=24)

plt.savefig('F:\\Seer_project\\figures\\CV_analysis\\LFQ_CV_pep.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[13]:


ter_np1=um['NP1'].tolist()
ter_np2=um['NP2'].tolist()
ter_np3=um['NP3'].tolist()
ter_np4=um['NP4'].tolist()
ter_np5=um['NP5'].tolist()
ter_dd=um['DD'].tolist()


# In[14]:


data_sorted_np1 = np.sort(ter_np1)
data_sorted_np2 = np.sort(ter_np2)
data_sorted_np3 = np.sort(ter_np3)
data_sorted_np4 = np.sort(ter_np4)
data_sorted_np5 = np.sort(ter_np5)
data_sorted_dd = np.sort(ter_dd)


# In[22]:


len(data_sorted_np1)


# In[26]:


def remove_nan(tarlist):
    return(np.array([item for item in tarlist if not math.isnan(item)]))
    


# In[28]:


cumulative_fraction_np1 = np.arange(1, len(remove_nan(data_sorted_np1)) + 1) / len(remove_nan(data_sorted_np1))
cumulative_fraction_np2 = np.arange(1, len(remove_nan(data_sorted_np2)) + 1) / len(remove_nan(data_sorted_np2))
cumulative_fraction_np3 = np.arange(1, len(remove_nan(data_sorted_np3)) + 1) / len(remove_nan(data_sorted_np3))
cumulative_fraction_np4 = np.arange(1, len(remove_nan(data_sorted_np4)) + 1) / len(remove_nan(data_sorted_np4))
cumulative_fraction_np5 = np.arange(1, len(remove_nan(data_sorted_np5)) + 1) / len(remove_nan(data_sorted_np5))
cumulative_fraction_dd  = np.arange(1, len(remove_nan(data_sorted_dd)) + 1) / len(remove_nan(data_sorted_dd))


# In[32]:


cumulative_fraction_np4


# In[40]:


plt.rcParams['figure.figsize']=(8,6)
sns.set_style("whitegrid")


plt.plot(remove_nan(data_sorted_np1), cumulative_fraction_np1, linewidth=3,label="NP1")
plt.plot(remove_nan(data_sorted_np2), cumulative_fraction_np2, linewidth=3,label="NP2")
plt.plot(remove_nan(data_sorted_np3), cumulative_fraction_np3, linewidth=3,label="NP3")
plt.plot(remove_nan(data_sorted_np4), cumulative_fraction_np4, linewidth=3,label="NP4")
plt.plot(remove_nan(data_sorted_np5), cumulative_fraction_np5, linewidth=3,label="NP5")
plt.plot(remove_nan(data_sorted_dd), cumulative_fraction_dd, linewidth=3,label="DD")

plt.xlabel('CV(Coefficient of Variation), %')
plt.ylabel('Cumulative Fraction')
#plt.title('Cumulative Fraction Plot')
plt.grid(True)
plt.legend()
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\CV_pro_lineplot.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[ ]:




