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
import venn
from matplotlib_venn import venn3
import seaborn as sns


# In[2]:


sample=[0,2,3,4]

def remove_dul(tar_list=sample):
    res=[]
    [res.append(x) for x in tar_list if x not in res]
    return(res)
def get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq"):
    '''find the csodiaq protein and peptide files to a list'''
    plasma_pro,plasma_pep = [],[]
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
        if file.endswith(".csv") and 'peptideFDR' in file:
            plasma_pep.append(os.path.join(dir, file))
    return([plasma_pro,plasma_pep])



# In[3]:


pro_list1=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[0]
pro_list2=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[0]
pro_list3=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[0]
pro_list4=get_pro_pep_files(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[0]
pro_list5=get_pro_pep_files(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[0]
pro_listdd=get_pro_pep_files(dir="E:\\yuming\\2023\\20230118\\csodiaq\DD")[0]


# In[4]:


pro_list5


# In[5]:


pep_list1=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[1]
pep_list2=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[1]
pep_list3=get_pro_pep_files(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[1]
pep_list4=get_pro_pep_files(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]
pep_list5=get_pro_pep_files(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]
pep_listdd=get_pro_pep_files(dir="E:\\yuming\\2023\\20230118\\csodiaq\DD")[1]


# In[6]:


pep_list2


# In[11]:


pro_rep1=pd.read_csv(pro_list2[0])['protein'].tolist()
pro_rep2=pd.read_csv(pro_list2[1])['protein'].tolist()
pro_rep3=pd.read_csv(pro_list2[2])['protein'].tolist()


# In[16]:


len(remove_dul(tar_list=pro_rep3))


# In[21]:


pep_rep1=pd.read_csv(pep_list3[3])['peptide'].tolist()
pep_rep2=pd.read_csv(pep_list3[4])['peptide'].tolist()
pep_rep3=pd.read_csv(pep_list3[5])['peptide'].tolist()


# In[7]:


def plot_venn_diagram(pro_list=pro_list5,pep_list=pep_list5, i=0,
                      savename1='venn5_protein_triplicates',
                     savename2='venn5_peptide_triplicates'):
    '''plot venn diagram from the triplicates of both proteins and peptides in all NPs and Neat plasma'''

    pro_rep1=pd.read_csv(pro_list[i])['protein'].tolist()
    pro_rep2=pd.read_csv(pro_list[i+1])['protein'].tolist()
    pro_rep3=pd.read_csv(pro_list[i+2])['protein'].tolist()

    pep_rep1=pd.read_csv(pep_list[i])['peptide'].tolist()
    pep_rep2=pd.read_csv(pep_list[i+1])['peptide'].tolist()
    pep_rep3=pd.read_csv(pep_list[i+2])['peptide'].tolist()

    #plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示
    matplotlib.rcParams['font.family'] = "Arial"
    plt.rcParams['figure.figsize']=(6,6)

    out1=venn3(subsets =[set(remove_dul(tar_list=pro_rep1)), set(remove_dul(tar_list=pro_rep2)),set(remove_dul(tar_list=pro_rep3))], 
          set_labels =('Rep1', 'Rep2','Rep3'),
          alpha=0.5,
          normalize_to=1)
    for text in out1.set_labels:
        text.set_fontsize(22)
    for text in out1.subset_labels:
        text.set_fontsize(22)
        text.set_color('black')
    #plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename1), dpi=800, bbox_inches='tight') 
    plt.show()

    out2=venn3(subsets =[set(remove_dul(tar_list=pep_rep1)), set(remove_dul(tar_list=pep_rep2)),set(remove_dul(tar_list=pep_rep3))], 
          set_labels =('Rep1', 'Rep2','Rep3'),
          alpha=0.5,
          normalize_to=1)
    for text in out2.set_labels:
        text.set_fontsize(22)
    for text in out2.subset_labels:
        text.set_fontsize(22)
        text.set_color('black')
    #plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename2), dpi=800)
    plt.show()


# In[10]:


plot_venn_diagram(pro_list=pro_list2,pep_list=pep_list2,
                      savename1='venn5_protein_triplicates',
                     savename2='venn5_peptide_triplicates')


# In[44]:


plot_venn_diagram(pro_list=pro_listdd,pep_list=pep_listdd,
                      savename1='venndd_protein_triplicates',
                     savename2='venndd_peptide_triplicates')


# In[ ]:





# In[ ]:





# In[23]:


from matplotlib_venn import venn3
matplotlib.rcParams['font.family'] = "Arial"
plt.rcParams['figure.figsize']=(6,6)

out=venn3(subsets =[set(remove_dul(tar_list=pep_rep1)), set(remove_dul(tar_list=pep_rep2)),set(remove_dul(tar_list=pep_rep3))], 
      set_labels =('Rep1', 'Rep2','Rep3'),
      alpha=0.5,
      normalize_to=1)
for text in out.set_labels:
    text.set_fontsize(22)
for text in out.subset_labels:
    text.set_fontsize(22)
    text.set_color('black')
#plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\venn_peptide_triplicates.svg', dpi=800)
plt.show()


# In[86]:


def get_quan_files_(path1=r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_proteins.csv',
                    path2=r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_peptides.csv'):
    pro_quan = pd.read_csv(path1)
    pep_quan = pd.read_csv(path2)
    pro_quan.set_index('Unnamed: 0',inplace=True)
    pep_quan.set_index('Unnamed: 0',inplace=True)
    pro_quan.columns=[item.split("/")[-1].split(".")[0] for item in pro_quan.columns]
    pep_quan.columns=[item.split("/")[-1].split(".")[0] for item in pep_quan.columns]
    pro_quan_log2=np.log2(pro_quan)
    pep_quan_log2=np.log2(pep_quan)
    return([pro_quan_log2,pep_quan_log2])

from scipy import stats
def calculate_pearson_r_squared(list1, list2):
    """Compute Pearson R squared from two lists."""
    # Compute Pearson correlation coefficient
    correlation, _ = stats.pearsonr(list1, list2)
    
    # Compute R squared (the square of the correlation)
    r_squared = correlation**2

    return r_squared


# In[77]:


get_quan_files_(path1=r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_proteins.csv',
                path2=r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_peptides.csv')[1]


# In[46]:


# plot Seaborn jointplot of the replicates 

pro_quan=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP3\\csodiaq\\common_proteins.csv')
pep_quan=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP3\\csodiaq\\common_peptides.csv')

pro_quan.set_index('Unnamed: 0',inplace=True)
pep_quan.set_index('Unnamed: 0',inplace=True)


# In[47]:


pro_quan.columns=[item.split("/")[-1].split(".")[0] for item in pro_quan.columns]
pep_quan.columns=[item.split("/")[-1].split(".")[0] for item in pep_quan.columns]


# In[48]:


pro_quan_log2=np.log2(pro_quan)
pep_quan_log2=np.log2(pep_quan)


# In[111]:


dfDD0=get_quan_files_(path1=r'E:\\yuming\\2023\\20230118\\csodiaq\\DD\\common_proteins.csv',
                path2=r'E:\\yuming\\2023\\20230118\\csodiaq\\DD\\common_peptides.csv')[0]
dfDD1=get_quan_files_(path1=r'E:\\yuming\\2023\\20230118\\csodiaq\\DD\\common_proteins.csv',
                path2=r'E:\\yuming\\2023\\20230118\\csodiaq\\DD\\common_peptides.csv')[1]
dfDD0


# In[112]:



plt.rcParams['font.family'] = 'Arial'
# Plot the jointplot
g=sns.jointplot(data=dfDD0, x='DD2X', y='DD2X2',kind="reg", truncate=False,
                      color="purple", height=6)
R_square0=round(calculate_pearson_r_squared(dfDD0['DD2X'].tolist(), dfDD0['DD2X2'].tolist()),3)
savename1='dd_pro_jointplot' +str(R_square0)
g.ax_joint.tick_params(labelsize=20)
g.set_axis_labels('Log2(intensity), replicate1', 'Log2(intensity), replicate2', fontsize=22)
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename1), dpi=800,bbox_inches='tight') 
plt.show()
print(R_square0)

# Plot the jointplot
g1=sns.jointplot(data=dfDD1,x='DD2X', y='DD2X2',kind="reg", truncate=False,
                      color="purple", height=6)
R_square1=round(calculate_pearson_r_squared(dfDD1['DD2X'].tolist(), dfDD1['DD2X2'].tolist()),3)
savename2='dd_pep_jointplot'+str(R_square1)
g1.ax_joint.tick_params(labelsize=20)
g1.set_axis_labels('Log2(intensity), replicate1', 'Log2(intensity), replicate2', fontsize=22)
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename2), dpi=800,bbox_inches='tight') 
plt.show()
print(R_square1)


# In[ ]:





# In[ ]:





# In[104]:



df0=get_quan_files_(path1=r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_proteins.csv',
                path2=r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_peptides.csv')[0]
df1=get_quan_files_(path1=r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_proteins.csv',
                path2=r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_peptides.csv')[1]
df0


# In[105]:



plt.rcParams['font.family'] = 'Arial'
# Plot the jointplot
g=sns.jointplot(data=df0, x='NP5_2', y='NP5_2_2',kind="reg", truncate=False,
                      color="purple", height=6)
R_square0=round(calculate_pearson_r_squared(df0['NP5_2'].tolist(), df0['NP5_2_2'].tolist()),3)
savename1='NP5_pro_jointplot' +str(R_square0)
g.ax_joint.tick_params(labelsize=20)
g.set_axis_labels('Log2(intensity), replicate1', 'Log2(intensity), replicate2', fontsize=22)
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename1), dpi=800,bbox_inches='tight') 
plt.show()
print(R_square0)

# Plot the jointplot
g1=sns.jointplot(data=df1,x='NP5_2', y='NP5_2_2',kind="reg", truncate=False,
                      color="purple", height=6)
R_square1=round(calculate_pearson_r_squared(df1['NP5_2'].tolist(), df1['NP5_2_2'].tolist()),3)
savename2='NP5_pep_jointplot'+str(R_square1)
g1.ax_joint.tick_params(labelsize=20)
g1.set_axis_labels('Log2(intensity), replicate1', 'Log2(intensity), replicate2', fontsize=22)
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\{}.svg'.format(savename2), dpi=800,bbox_inches='tight') 
plt.show()
print(R_square1)


# In[ ]:





# In[ ]:





# In[ ]:





# In[10]:


np1_pro=pd.read_csv(plasma_pro[1])['protein'].tolist()
np1_pro2=pd.read_csv(plasma_pro[2])['protein'].tolist()
np2_pro=pd.read_csv(plasma_pro[3])['protein'].tolist()
np3_pro=pd.read_csv(plasma_pro[4])['protein'].tolist()
np3_pro2=pd.read_csv(plasma_pro[5])['protein'].tolist()
np4_pro=pd.read_csv(plasma_pro[6])['protein'].tolist()
np5_pro=pd.read_csv(plasma_pro[7])['protein'].tolist()
DD_pro=pd.read_csv(plasma_pro[0])['protein'].tolist()


# In[19]:


np1_pep=pd.read_csv(plasma_pep[1])['peptide'].tolist()
np1_pep2=pd.read_csv(plasma_pep[2])['peptide'].tolist()
np2_pep=pd.read_csv(plasma_pep[3])['peptide'].tolist()
np3_pep=pd.read_csv(plasma_pep[4])['peptide'].tolist()
np3_pep2=pd.read_csv(plasma_pep[5])['peptide'].tolist()
np4_pep=pd.read_csv(plasma_pep[6])['peptide'].tolist()
np5_pep=pd.read_csv(plasma_pep[7])['peptide'].tolist()
DD_pep=pd.read_csv(plasma_pep[0])['peptide'].tolist()


# In[ ]:





# In[20]:


pro_number=[]
for i in range(len(plasma_pro)):
    df=pd.read_csv(plasma_pro[i]) 
    pro_number.append(len(remove_dul(tar_list=df['leadingProtein'].tolist())))
pro_number  


# In[24]:


pep_number=[]
for i in range(len(plasma_pep)):
    df=pd.read_csv(plasma_pep[i]) 
    pep_number.append(len(remove_dul(tar_list=df['peptide'].tolist())))
pep_number   


# In[25]:


prolist=[254, 307, 393, 341, 348, 361]
peplist=[1897,2016, 2137, 2274, 2273, 2535]

name=['DD','NP1','NP2','NP3','NP4','NP5']


# In[26]:


+


# In[20]:


all_NPs=remove_dul(tar_list=np1_pro)+remove_dul(tar_list=np2_pro)+remove_dul(tar_list=np3_pro)+remove_dul(tar_list=np4_pro)+remove_dul(tar_list=np5_pro)


# In[21]:


len(remove_dul(tar_list=all_NPs))


# In[22]:


from matplotlib_venn import venn2
my_dpi=150
plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示

out=venn2(subsets =[set(remove_dul(tar_list=DD_pro)), set(remove_dul(tar_list=all_NPs))], 
      set_labels =('DD', 'all_NPs'),
      alpha=0.5,
      normalize_to=1)
for text in out.set_labels:
    text.set_fontsize(14)
for text in out.subset_labels:
    text.set_fontsize(12)
    text.set_color('black')
#plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\figures\\venn diagram\\45V-50V.svg', dpi=600) 
plt.show()


# In[27]:


K1=set(remove_dul(tar_list=np1_pro))& set(remove_dul(tar_list=np2_pro))
K2=K1&set(remove_dul(tar_list=np3_pro))
K3=K2&set(remove_dul(tar_list=np4_pro))
K4=K3&set(remove_dul(tar_list=np5_pro))

len(K4)



# In[28]:


all_NPs_shared=K4


# In[32]:


from matplotlib_venn import venn2
my_dpi=150
plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示

out=venn2(subsets =[set(remove_dul(tar_list=DD_pro)), K4], 
      set_labels =('DD', 'all_NPs_shared'),
      alpha=0.5,
      normalize_to=1)
for text in out.set_labels:
    text.set_fontsize(14)
for text in out.subset_labels:
    text.set_fontsize(12)
    text.set_color('black')
#plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\figures\\venn diagram\\45V-50V.svg', dpi=600) 
plt.show()


# In[33]:


labels = venn.get_labels([set(remove_dul(tar_list=np1_pro)), 
                          set(remove_dul(tar_list=np2_pro)), 
                          set(remove_dul(tar_list=np3_pro)), 
                          set(remove_dul(tar_list=np4_pro)), 
                          set(remove_dul(tar_list=np5_pro))], fill=['number'])
fig, ax = venn.venn5(labels, names=['NP1','NP2','NP3','NP4', 'NP5'])
fig.show()


# In[34]:


labels = venn.get_labels([set(remove_dul(tar_list=np1_pep)), 
                          set(remove_dul(tar_list=np2_pep)), 
                          set(remove_dul(tar_list=np3_pep)), 
                          set(remove_dul(tar_list=np4_pep)), 
                          set(remove_dul(tar_list=np5_pep))], fill=['number'])
fig, ax = venn.venn5(labels, names=['NP1','NP2','NP3','NP4', 'NP5'])
fig.show()


# In[35]:


labels = venn.get_labels([set(remove_dul(tar_list=DD_pro)),
                          set(remove_dul(tar_list=np1_pro)), 
                          set(remove_dul(tar_list=np2_pro)), 
                          set(remove_dul(tar_list=np3_pro)), 
                          set(remove_dul(tar_list=np4_pro)), 
                          set(remove_dul(tar_list=np5_pro))], fill=['number'])
fig, ax = venn.venn6(labels, names=['DD','NP1','NP2','NP3','NP4', 'NP5'])
fig.show()


# In[36]:


labels = venn.get_labels([set(remove_dul(tar_list=DD_pep)),
                          set(remove_dul(tar_list=np1_pep)), 
                          set(remove_dul(tar_list=np2_pep)), 
                          set(remove_dul(tar_list=np3_pep)), 
                          set(remove_dul(tar_list=np4_pep)), 
                          set(remove_dul(tar_list=np5_pep))], fill=['number'])
fig, ax = venn.venn6(labels, names=['DD','NP1','NP2','NP3','NP4', 'NP5'])
fig.show()


# In[76]:





# In[ ]:





# In[27]:


# protein identifications from different NPs 

sns.set_theme(style="ticks")
plt.figure(figsize=(9,5.5),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.43
plt.bar(x=name,height=prolist,color='firebrick',alpha=0.8,width=bar_width)
plt.xlabel('treatment',fontsize=24)
plt.ylabel('Number of proteins',fontsize=24)
plt.title('DISPA',fontsize=26)
plt.xticks(fontsize=20)
#plt.xticks(rotation=20)
plt.yticks(fontsize=22)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,400, 100))
#plt.xlim(7.8,23)
plt.ylim(0,400)
plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
#plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\robust and reproducibility\\CV analysis_protein.svg', dpi=600) 
plt.show()


# In[29]:


# peptide identifications from different NPs 

sns.set_theme(style="ticks")
plt.figure(figsize=(9,5.5),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.43
plt.bar(x=name,height=peplist,color='orange',alpha=0.8,width=bar_width)
plt.xlabel('treatment',fontsize=24)
plt.ylabel('Number of peptides',fontsize=24)
#plt.title('10 peptide quantities',fontsize=26)
plt.xticks(fontsize=20)
#plt.xticks(rotation=20)
plt.yticks(fontsize=22)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,2500, 500))
#plt.xlim(7.8,23)
plt.ylim(0,2600)
plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
#plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\robust and reproducibility\\CV analysis_protein.svg', dpi=600) 
plt.show()


# In[ ]:





# In[4]:


# protein IDS from DDA LC-MS/MS
proID_DDA=[269, 302, 421, 403, 418, 437]

name=['DD','NP1','NP2','NP3','NP4','NP5']

sns.set_theme(style="ticks")
plt.figure(figsize=(9,5.5),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.43
plt.bar(x=name,height=proID_DDA,color='firebrick',alpha=0.8,width=bar_width)
#plt.xlabel('treatment',fontsize=24)
plt.ylabel('Number of proteins',fontsize=24)
plt.title('DDA_LC_MS/MS',fontsize=26)
plt.xticks(fontsize=20)
#plt.xticks(rotation=20)
plt.yticks(fontsize=22)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,500, 100))
#plt.xlim(7.8,23)
plt.ylim(0,500)
plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\all_NPs_LC_DDA.svg', dpi=600) 
plt.show()



# In[3]:


df=pd.read_excel(r'F:\\Seer_project\\Supplementary Table 2.xlsx')


# In[4]:


np1=df[df['NP1'] !=1]["pro_names"].tolist()
np2=df[df['NP2'] !=1]["pro_names"].tolist()
np3=df[df['NP3'] !=1]["pro_names"].tolist()
np4=df[df['NP4'] !=1]["pro_names"].tolist()
np5=df[df['NP5'] !=1]["pro_names"].tolist()


# In[5]:


dd=df[df['DD'] !=1]["pro_names"].tolist()


# In[11]:


#plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示
from matplotlib_venn import venn2
my_dpi=150
plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示

out=venn2(subsets =[set(np1+np2+np3+np4+np5), set(dd)], 
      set_labels =('NPs', 'DD'),
      alpha=0.5,
      normalize_to=1)
for text in out.set_labels:
    text.set_fontsize(16)
for text in out.subset_labels:
    text.set_fontsize(16)
    text.set_color('black')
plt.savefig('F:\\Seer_project\\figures\\Venn diagram\\DISPA_TARGETED_VENN.svg', dpi=600) 
plt.show() 


# In[ ]:




