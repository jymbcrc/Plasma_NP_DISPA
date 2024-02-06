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
sample=[0,2,3,4]
def remove_dul(tar_list=sample):
    res=[]
    [res.append(x) for x in tar_list if x not in res]
    return(res)
def procount(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'):
    """calculate the protein IDs from all csodiaq files from a path"""
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
    """calculate the peptide IDs from all csodiaq files from a path"""
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


# In[3]:


plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]


# In[6]:


dfpro=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]


# In[7]:


def getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]):
    df1ug=dftar.iloc[:,:3]
    df05ug=dftar.iloc[:,3:6]
    df025ug=dftar.iloc[:,6:9]
    df0125ug=dftar.iloc[:,9:12]
    df11=pd.DataFrame(df1ug.mean(axis=1),index=df1ug.index,columns=['1ug'])
    df05=pd.DataFrame(df05ug.mean(axis=1),index=df1ug.index,columns=['05ug'])
    df025=pd.DataFrame(df025ug.mean(axis=1),index=df1ug.index,columns=['025ug'])
    df0125=pd.DataFrame(df0125ug.mean(axis=1),index=df1ug.index,columns=['0125ug'])
    result = pd.concat([df0125,df025,df05,df11], axis=1)
    maxdf=result.max(axis=1)
    norm_df = result.divide(maxdf, axis='index')
    return(norm_df)


# In[8]:


def plot_figure(filename=getmean(dftar=dfpro),title='protein quantities_NP4'):
    sns.set(font_scale=2)
    sns.set_style('white')
    plt.rcParams['figure.figsize'] = 6,6
    spep = filename
    for i in range(0, len(spep)):
        plt.plot([0.125,0.25,0.5,1], spep.iloc[i])   
    plt.plot([0.125,0.25,0.5,1], spep.mean(), color='black', linewidth=8,alpha=1)
    plt.xlim(0,1)
    plt.ylim(0,1)
    #plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=26)
    plt.ylabel('relative quantity, measured',fontsize=26)
    plt.title(title,fontsize=26)
    #plt.savefig('D:/Project 2 DIrect infusion shotgun proteome analysis/quantification analysis/significant peptides.svg', bbox_inches='tight')
    fig_path= 'E:\\project6_rapid profiling plasma peptides\\figures\\'
    plt.savefig(fig_path + "%s.svg" % title,dpi=600, bbox_inches='tight')
    


# In[9]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]),title='All protein quantities_NP4')
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[0]),title='All peptide quantities_NP4')


# In[8]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]),title='All protein quantities_NP5') 
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[0]),title='All peptide quantities_NP5')

#'''0 is pep,1 is pro'''


# In[10]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[1]),title='All protein quantities_NP1')
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[0]),title='All peptide quantities_NP1')


# In[11]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[1]),title='All protein quantities_NP3')
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[0]),title='All peptide quantities_NP3')


# In[12]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[1]),title='All protein quantities_NP2')
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[0]),title='All peptide quantities_NP2')


# In[13]:


plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]),title='All protein quantities_DD')
plot_figure(filename=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0]),title='All peptide quantities_DD')


# # compare the CV among replicates in different concentrations for DD and al NPs

# In[9]:


def calculate_pearson_r_squared(list1, list2):
    from scipy import stats
    """Compute Pearson R squared from two lists."""
    # Compute Pearson correlation coefficient
    correlation, _ = stats.pearsonr(list1, list2)
    
    # Compute R squared (the square of the correlation)
    r_squared = correlation**2

    return r_squared


# In[10]:


peptide_dd=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0]
protein_dd=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]

peptide_np1=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[0]
protein_np1=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[1]

peptide_np2=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[0]
protein_np2=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[1]

peptide_np3=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[0]
protein_np3=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[1]

peptide_np4=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[0]
protein_np4=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]

peptide_np5=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[0]
protein_np5=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]


# In[62]:



# procount(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")
# procount(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")
# procount(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")
# procount(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")
# procount(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")
# procount(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")


# In[ ]:





# In[11]:


def cal_corr_matrix(num=0,df=np.log2(protein_dd)):
    '''calculate pearson R square from replicates in each np and dd, output a matrix'''
    list1=df.iloc[:,num].tolist()
    list2=df.iloc[:,num+1].tolist()
    list3=df.iloc[:,num+2].tolist()
    lists=[list1,list2,list3]

    # 计算Pearson相关系数
    correlation_matrix = np.zeros((len(lists), len(lists)))

    for i in range(len(lists)):
        for j in range(i+1, len(lists)):
            corr = round(calculate_pearson_r_squared(lists[i], lists[j]),3)
            #corr=np.corrcoef(lists[i], lists[j])[0, 1]
            correlation_matrix[i][j] = corr
            correlation_matrix[j][i] = corr

    for i in range(len(lists)):
        corr = round(calculate_pearson_r_squared(lists[i], lists[i]),3)
        correlation_matrix[i][i] = corr

    return(correlation_matrix)


# In[57]:


# cal_corr_matrix(num=0,df=np.log2(protein_dd))


# In[58]:


# cal_corr_matrix(num=0,df=np.log2(protein_np1))


# In[59]:


# cal_corr_matrix(num=0,df=np.log2(protein_np2))


# In[60]:


# cal_corr_matrix(num=0,df=np.log2(protein_np3))


# In[61]:


# cal_corr_matrix(num=0,df=np.log2(protein_np4))


# In[62]:


# cal_corr_matrix(num=0,df=np.log2(protein_np5))


# In[34]:


# 创建热图
con_num=9
filename=protein_dd
plt.figure(figsize=(8, 6))
plt.rcParams['font.family'] = 'Arial'

plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
# Create a heatmap
ax = sns.heatmap(cal_corr_matrix(num=con_num,df=np.log2(filename)), annot=True, cmap='viridis', fmt=".3f", annot_kws={"size": 30},
            xticklabels=['Rep1','Rep2','Rep3'],yticklabels=['Rep1','Rep2','Rep3'])

# Access the colorbar and modify its font size
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=22)  # Change '14' to your desired font size

#cbar.set_ticks([0.94,0.96,0.98, 1])

fig_path= 'F:\\Seer_project\\figures\\Correlations of the protein\\'
figname= 'dd_9'
plt.savefig(fig_path + "%s.svg" % figname,dpi=800, bbox_inches='tight')  

plt.show()


# In[155]:


#fimat=np.zeros((18, 18))

fimat_pro = np.full((18, 18), np.nan)
fimat_pep = np.full((18, 18), np.nan)


# In[156]:


for i in [0,1,2]:
    for j in [0,1,2]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_dd))[i,j]
for i in [3,4,5]:
    for j in [3,4,5]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_np1))[i-3,j-3]
for i in [6,7,8]:
    for j in [6,7,8]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_np2))[i-6,j-6]
for i in [9,10,11]:
    for j in [9,10,11]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_np3))[i-9,j-9]
for i in [12,13,14]:
    for j in [12,13,14]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_np4))[i-12,j-12]
for i in [15,16,17]:
    for j in [15,16,17]:
        fimat_pep[i,j]=cal_corr_matrix(num=0,df=np.log2(peptide_np5))[i-15,j-15]
#fimat_pep


# In[157]:


for i in [0,1,2]:
    for j in [0,1,2]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_dd))[i,j]
for i in [3,4,5]:
    for j in [3,4,5]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_np1))[i-3,j-3]
for i in [6,7,8]:
    for j in [6,7,8]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_np2))[i-6,j-6]
for i in [9,10,11]:
    for j in [9,10,11]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_np3))[i-9,j-9]
for i in [12,13,14]:
    for j in [12,13,14]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_np4))[i-12,j-12]
for i in [15,16,17]:
    for j in [15,16,17]:
        fimat_pro[i,j]=cal_corr_matrix(num=0,df=np.log2(protein_np5))[i-15,j-15]
#fimat_pro


# In[159]:


# 创建热图
from matplotlib.patches import Rectangle

fimat=fimat_pep
plt.figure(figsize=(8, 6))
plt.rcParams['font.family'] = 'Arial'

ax =sns.heatmap(fimat, annot=False,fmt=".3f",cmap='viridis',linewidth=0.6, linecolor='grey',
               xticklabels=['R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3'],
               yticklabels=['R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3','R1','R2','R3'])
rect = Rectangle((0,0), fimat.shape[1], fimat.shape[0], linewidth=2, edgecolor='black', facecolor='none')
ax.add_patch(rect)
# Set x-tick and y-tick labels font properties
font_properties = {
    'size': 12,
    #'weight': 'bold',
    'style': 'normal',
    'color': 'black'
}
# Setting x-tick labels
for label in ax.get_xticklabels():
    label.set_fontproperties('Times New Roman')
    label.set(**font_properties)
# Setting y-tick labels
for label in ax.get_yticklabels():
    label.set_fontproperties('Times New Roman')
    label.set(**font_properties)
#plt.title("Pearson Correlation of Proteome")
#plt.xlabel("List Index")
#plt.ylabel("List Index")
#plt.savefig(r'E:\\project6_rapid profiling plasma peptides\\New folder\\figures\\protein_corre_heatmap.svg',dpi=800,bbox_inches='tight')
plt.savefig(r'E:\\project6_rapid profiling plasma peptides\\New folder\\figures\\peptide_corre_heatmap.svg',dpi=800,bbox_inches='tight')
plt.show()


# In[12]:


plt.rcParams['font.family'] = 'Arial'

def plot_pep_figure2(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0]), 
                     title='peptide quantities_NP4'):
    DD_pep_mean=pepfile.mean()
    DD_pep_std=pepfile.std()
    DD_pep_sem=pepfile.sem()
   
    plt.rcParams['figure.figsize'] = 6,6
    fig, ax = plt.subplots()
    #ax.plot([0.125,0.25,0.5,1], allpep.iloc[0:4,1],color='orange', linewidth=2,label='All peptides')
    ax.fill_between([0.125,0.25,0.5,1], DD_pep_mean+DD_pep_std, DD_pep_mean-DD_pep_std ,alpha=0.5, facecolor='dimgray')
    ax.plot([0.125,0.25,0.5,1], [0.125,0.25,0.5,1],color='blue', linewidth=3,label='X=Y')
    ax.errorbar([0.125,0.25,0.5,1], DD_pep_mean,yerr=np.array(DD_pep_sem), color='dimgray', linewidth=2,alpha=1,label='All peptides')
    #plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
    plt.ylabel('relative quantity, measured',fontsize=26)
    plt.title(title,fontsize=28)
    plt.legend(loc='lower right',fontsize=18)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    fig_path= 'E:\\project6_rapid profiling plasma peptides\\figures\\'
    plt.savefig(fig_path + "%s.svg" % title,dpi=600, bbox_inches='tight') 
    
def plot_pro_figure2(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]), 
                     title='protein quantities_NP4'):
 
    DD_pro_mean=profile.mean()
    DD_pro_std=profile.std()
    DD_pro_sem=profile.sem() 
    
    plt.rcParams['figure.figsize'] = 6,6
    fig, ax = plt.subplots()
    ax.fill_between([0.125,0.25,0.5,1], DD_pro_mean+DD_pro_std, DD_pro_mean-DD_pro_std ,alpha=0.5, facecolor='dimgray')
    ax.plot([0.125,0.25,0.5,1], [0.125,0.25,0.5,1],color='blue', linewidth=3,label='X=Y')
    ax.errorbar([0.125,0.25,0.5,1], DD_pro_mean,yerr=np.array(DD_pro_sem), color='dimgray', linewidth=2,alpha=1,label='All proteins')
    #plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
    plt.ylabel('relative quantity, measured',fontsize=28)
    plt.title(title,fontsize=26)
    plt.legend(loc='lower right',fontsize=18)
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    fig_path= 'E:\\project6_rapid profiling plasma peptides\\figures\\'
    plt.savefig(fig_path + "%s.svg" % title,dpi=600, bbox_inches='tight')    


# In[23]:


def plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0]), 
                     title='peptide quantities_NP4'):
    DD_pep_mean=pepfile.mean()
    DD_pep_std=pepfile.std()
    DD_pep_sem=pepfile.sem()

    plt.rcParams['figure.figsize'] = 6,6
    fig, ax = plt.subplots()
    #ax.plot([0.125,0.25,0.5,1], allpep.iloc[0:4,1],color='orange', linewidth=2,label='All peptides')
    ax.fill_between([0.05,0.1,0.2,0.4], DD_pep_mean+DD_pep_std, DD_pep_mean-DD_pep_std ,alpha=0.5, facecolor='#6387AB')
    ax.plot([0.05,0.1,0.2,0.4], [0.125,0.25,0.5,1],color='#E419A7', linewidth=3,label='Ideal')
    ax.errorbar([0.05,0.1,0.2,0.4], DD_pep_mean,yerr=np.array(DD_pep_sem), color='dimgray', linewidth=3,alpha=1,label='All peptides')
    #plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
    plt.ylabel('relative quantity',fontsize=26)
    plt.xlabel('µg/µl,injected',fontsize=26)
    plt.title(title,fontsize=28)
    plt.legend(loc='lower right',fontsize=18)
    plt.xlim(0,0.401)
    plt.xticks(np.arange(0, 0.4001, 0.1))
    plt.ylim(0,1)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    fig_path= 'F:\\Seer_project\\figures\\quantification\\'
    plt.savefig(fig_path + "%s.svg" % title,dpi=600, bbox_inches='tight')
    
def plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]), 
                     title='protein quantities_NP4'):
    DD_pro_mean=profile.mean()
    DD_pro_std=profile.std()
    DD_pro_sem=profile.sem() 
    
    plt.rcParams['figure.figsize'] = 6,6
    fig, ax = plt.subplots()
    ax.fill_between([0.05,0.1,0.2,0.4], DD_pro_mean+DD_pro_std, DD_pro_mean-DD_pro_std ,alpha=0.5, facecolor='#6387AB')
    ax.plot([0.05,0.1,0.2,0.4], [0.125,0.25,0.5,1],color='#E419A7', linewidth=3,label='Ideal')
    ax.errorbar([0.05,0.1,0.2,0.4], DD_pro_mean,yerr=np.array(DD_pro_sem), color='dimgray', linewidth=3,alpha=1,label='All proteins')
    #plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
    plt.ylabel('relative quantity',fontsize=28)
    plt.xlabel('µg/µl,injected',fontsize=26)
    plt.title(title,fontsize=26)
    plt.legend(loc='lower right',fontsize=18)
    plt.xlim(0,0.401)
    plt.xticks(np.arange(0, 0.4001, 0.1))
    plt.ylim(0,1)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    fig_path= 'F:\\Seer_project\\figures\\quantification\\'
    plt.savefig(fig_path + "%s.svg" % title,dpi=600, bbox_inches='tight')      


# In[15]:


df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0])


# In[16]:



def getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1])):
    from scipy.stats import pearsonr, spearmanr, linregress
    pearsonrs = []
    pearsonps= []
    spearmanrs = []
    spearmanps = []
    slopes = []
    r_square=[]
    for i in range(0, len(df_target)):
        x = np.array([0.05, 0.1, 0.2, 0.4])
        y = df_target.iloc[i]
        pr, pp = pearsonr(x, y)
        pearsonrs.append(pr)
        pearsonps.append(pp)
        sr, sp = spearmanr(x, y)  
        spearmanrs.append(sr)
        spearmanps.append(sp)
        result = linregress(x, y)
        slopes.append(result.slope)
        r_square.append(result.rvalue)
    regression_metrics = pd.DataFrame({'r':pearsonrs,'pvalp':pearsonps, 'rho':spearmanrs, 'pvals':spearmanps, 'slopes':slopes, "r_square":r_square}, index=df_target.index)
    return(regression_metrics)


# In[ ]:





# In[24]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0]), 
                     title='peptide quantities_DD')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]), 
                     title='protein quantities_DD')


# In[41]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1]))['r'].tolist()))


# In[40]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230118\\csodiaq\\DD")


# In[25]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[0]), 
                     title='peptide quantities_NP1')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[1]), 
                     title='protein quantities_NP1')


# In[42]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[1]))['r'].tolist()))


# In[39]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230129\\NP1\\csodiaq")


# In[26]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[0]), 
                     title='peptide quantities_NP2')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[1]), 
                     title='protein quantities_NP2')


# In[43]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[1]))['r'].tolist()))


# In[38]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230129\\NP2\\csodiaq")


# In[27]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[0]), 
                     title='peptide quantities_NP3')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[1]), 
                     title='protein quantities_NP3')


# In[44]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")[1]))['r'].tolist()))


# In[37]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq")


# In[28]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[0]), 
                     title='peptide quantities_NP4')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]), 
                     title='protein quantities_NP4')


# In[45]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir ="E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1]))['r'].tolist()))


# In[36]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230130\\NP4\\csodiaq")


# In[34]:



def cal_good_quan_ratio(path="E:\\yuming\\2023\\20230130\\NP4\\csodiaq"):
    pi=getpval(df_target=getmean(dftar=plasma_file(dir = path)[1]))
    return(len(pi[(pi['pvalp']<0.05) & (pi['slopes']>0)])/len(pi))


# In[29]:


plot_pep_figure(pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[0]), 
                     title='peptide quantities_NP5')
plot_pro_figure(profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]), 
                     title='protein quantities_NP5')


# In[46]:


print(np.mean(getpval(df_target=getmean(dftar=plasma_file(dir ="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[1]))['r'].tolist()))


# In[35]:


cal_good_quan_ratio(path="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")


# In[27]:


pepfile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[0])
profile=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230118\\csodiaq\\DD")[1])


# In[28]:


DD_pep_mean=pepfile.mean()
DD_pro_mean=profile.mean()
DD_pro_std=profile.std()
DD_pep_std=pepfile.std()
DD_pep_sem=pepfile.sem()
DD_pro_sem=profile.sem()


# In[29]:


plt.rcParams['figure.figsize'] = 6,6
fig, ax = plt.subplots()

#ax.plot([0.125,0.25,0.5,1], allpep.iloc[0:4,1],color='orange', linewidth=2,label='All peptides')
ax.fill_between([0.125,0.25,0.5,1], DD_pep_mean+DD_pep_std, DD_pep_mean-DD_pep_std ,alpha=0.5, facecolor='dimgray')

ax.plot([0.125,0.25,0.5,1], [0.125,0.25,0.5,1],color='blue', linewidth=3,label='X=Y')
ax.errorbar([0.125,0.25,0.5,1], DD_pep_mean,yerr=np.array(DD_pep_sem), color='dimgray', linewidth=2,alpha=1,label='All peptides')

#plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
plt.ylabel('relative quantity, measured',fontsize=28)
plt.title('Peptide quantities/DD',fontsize=28)
plt.legend(loc='lower right',fontsize=18)
plt.xlim(0,1)
plt.ylim(0,1)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
#plt.savefig('D:/Project 2 DIrect infusion shotgun proteome analysis/quantification analysis/peptides quantities.svg', bbox_inches='tight')


# In[ ]:





# In[ ]:





# In[8]:


dfpep=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP3\\csodiaq\\common_peptides.csv')
dfpro=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP3\\csodiaq\\common_proteins.csv')


# In[41]:


def getmean(dftar=dfpro):
    df1ug=dftar.iloc[:,:3]
    df05ug=dftar.iloc[:,3:6]
    df025ug=dftar.iloc[:,6:9]
    df0125ug=dftar.iloc[:,9:12]
    df11=pd.DataFrame(df1ug.mean(axis=1),index=df1ug.index,columns=['1ug'])
    df05=pd.DataFrame(df05ug.mean(axis=1),index=df1ug.index,columns=['05ug'])
    df025=pd.DataFrame(df025ug.mean(axis=1),index=df1ug.index,columns=['025ug'])
    df0125=pd.DataFrame(df0125ug.mean(axis=1),index=df1ug.index,columns=['0125ug'])
    result = pd.concat([df0125,df025,df05,df11], axis=1)
    maxdf=result.max(axis=1)
    norm_df = result.divide(maxdf, axis='index')
    return(norm_df)


# In[33]:


DD_pep_mean=getmean(dftar=dfpep).mean().tolist()


# In[34]:


DD_pro_mean=getmean(dftar=dfpro).mean()


# In[35]:


DD_pro_std=getmean(dftar=dfpro).std()


# In[36]:


DD_pep_std=getmean(dftar=dfpep).std()


# In[37]:


DD_pep_sem=getmean(dftar=dfpep).sem()


# In[38]:


DD_pro_sem=getmean(dftar=dfpro).sem()


# In[42]:


sns.set(font_scale=2)
sns.set_style('white')

plt.rcParams['figure.figsize'] = 6,6
spep = getmean(dftar=dfpep)
for i in range(0, len(spep )):
    plt.plot([0.125,0.25,0.5,1], spep.iloc[i])   
plt.plot([0.125,0.25,0.5,1], spep.mean(), color='black', linewidth=8,alpha=1)
plt.xlim(0,1)
plt.ylim(0,1)
#plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=26)
plt.ylabel('relative quantity, measured',fontsize=26)
plt.title('peptide quantities/DD',fontsize=26)
#plt.savefig('D:/Project 2 DIrect infusion shotgun proteome analysis/quantification analysis/significant peptides.svg', bbox_inches='tight')


# In[45]:


plt.rcParams['figure.figsize'] = 6,6
fig, ax = plt.subplots()

#ax.plot([0.125,0.25,0.5,1], allpep.iloc[0:4,1],color='orange', linewidth=2,label='All peptides')
ax.fill_between([0.125,0.25,0.5,1], DD_pep_mean+DD_pep_std, DD_pep_mean-DD_pep_std ,alpha=0.5, facecolor='dimgray')

ax.plot([0.125,0.25,0.5,1], [0.125,0.25,0.5,1],color='blue', linewidth=3,label='X=Y')
ax.errorbar([0.125,0.25,0.5,1], DD_pep_mean,yerr=np.array(DD_pep_sem), color='dimgray', linewidth=2,alpha=1,label='All peptides')

#plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
plt.ylabel('relative quantity, measured',fontsize=28)
plt.title('Peptide quantities/DD',fontsize=28)
plt.legend(loc='lower right',fontsize=18)
plt.xlim(0,1)
plt.ylim(0,1)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
#plt.savefig('D:/Project 2 DIrect infusion shotgun proteome analysis/quantification analysis/peptides quantities.svg', bbox_inches='tight')


# In[46]:


plt.rcParams['figure.figsize'] = 6,6
fig, ax = plt.subplots()

#ax.plot([0.125,0.25,0.5,1], allpep.iloc[0:4,1],color='orange', linewidth=2,label='All peptides')
ax.fill_between([0.125,0.25,0.5,1], DD_pro_mean+DD_pro_std, DD_pro_mean-DD_pro_std ,alpha=0.5, facecolor='dimgray')


ax.plot([0.125,0.25,0.5,1], [0.125,0.25,0.5,1],color='blue', linewidth=3,label='X=Y')
ax.errorbar([0.125,0.25,0.5,1], DD_pro_mean,yerr=np.array(DD_pro_sem), color='dimgray', linewidth=2,alpha=1,label='All peptides')

#plt.xlabel(u"\u03bcg/\u03bcL injected",fontsize=28)
plt.ylabel('relative quantity, measured',fontsize=28)
plt.title('Protein quantities/DD',fontsize=28)
plt.legend(loc='lower right',fontsize=18)
plt.xlim(0,1)
plt.ylim(0,1)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
#plt.savefig('D:/Project 2 DIrect infusion shotgun proteome analysis/quantification analysis/peptides quantities.svg', bbox_inches='tight')


# In[4]:


dfpep_np1=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_peptides.csv')
dfpro_np1=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP1\\csodiaq\\common_proteins.csv')






# In[6]:


dfpro_np1.set_index('Unnamed: 0',inplace=True)


# In[26]:


df=dfpro_np1.iloc[:,0:4]


# In[27]:


df


# In[59]:


dfpep_np2=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP2\\csodiaq\\common_peptides.csv')
dfpro_np2=pd.read_csv(r'E:\\yuming\\2023\\20230129\\NP2\\csodiaq\\common_proteins.csv')


# In[61]:


dfpep_np2.set_index('Unnamed: 0',inplace=True)
dfpro_np2.set_index('Unnamed: 0',inplace=True)


# In[62]:


dfpro_np2


# In[ ]:



    
    


# In[30]:


from scipy.stats import pearsonr


# In[34]:


oo=[1,2,3,5]
pp=[1,2,3,4.7]


# In[35]:


pearsonr(oo, pp)


# In[154]:


dfpep_np4=pd.read_csv(r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_peptides.csv')
dfpro_np4=pd.read_csv(r'E:\\yuming\\2023\\20230130\\NP5\\csodiaq\\common_proteins.csv')
dfpro_np4.set_index('Unnamed: 0',inplace=True)
dfpep_np4.set_index('Unnamed: 0',inplace=True)


# In[14]:


## now get the spearman and slope from each protein
def getpval(df_target=getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[1])):
    from scipy.stats import pearsonr, spearmanr, linregress
    pearsonrs = []
    pearsonps= []
    spearmanrs = []
    spearmanps = []
    slopes = []
    for i in range(0, len(df_target)):
        x = np.array([0.125, 0.25, 0.5, 1])
        y = df_target.iloc[i]
        pr, pp = pearsonr(x, y)
        pearsonrs.append(pr)
        pearsonps.append(pp)
        sr, sp = spearmanr(x, y)  
        spearmanrs.append(sr)
        spearmanps.append(sp)
        result = linregress(x, y)
        slopes.append(result.slope)
    regression_metrics = pd.DataFrame({'r':pearsonrs,'pvalp':pearsonps, 'rho':spearmanrs, 'pvals':spearmanps, 'slopes':slopes}, index=df_target.index)
    return(regression_metrics)

def rand_jitter(arr):
    stdev = .02 * (max(arr) - min(arr))
    return arr + np.random.randn(len(arr)) * stdev

def jitter(x, y, s=20, c='b', marker='o', cmap=None, norm=None, vmin=None, vmax=None, alpha=None, linewidths=None, verts=None, hold=None, **kwargs):
    return scatter(rand_jitter(x), rand_jitter(y), s=s, c=c, marker=marker, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, alpha=alpha, linewidths=linewidths, **kwargs)

from matplotlib.pyplot import scatter


# In[74]:


def typical_quan(path="E:\\yuming\\2023\\20230130\\NP5\\csodiaq",i=1):
    
    df=plasma_file(dir = path)[i]

    df_pval=getpval(df_target=getmean(dftar=df))
    wo=df_pval.sort_values(by='slopes',ascending=False)
    woai=wo.iloc[:20,:]
    namelist=woai.index.tolist()
    if i ==0: 
        for pepname in namelist:
            sap_pro=df.loc[pepname]
            k=sap_pro.to_frame()
            pp=k.iloc[:12,0].to_frame()
            groups=[0.4, 0.4, 0.4, 0.2, 0.2, 0.2,0.1,0.1,0.1,0.05,0.05,0.05]
            pp['groups'] = groups
            plt.rcParams['figure.figsize'] = 5,5
            categories = np.array([ 0, 0, 0, 1, 1, 1, 2, 2,2,3,3,3])
            colormap = np.array(['#5ec962', '#21918c', '#3b528b', '#440154'])
            jitter(pp.groups.values, pp.iloc[:12,0], s=100, c = colormap[categories])
            plt.xlabel(u"\u03bcg/\u03bcL injected")
            plt.ylabel('relative quantity, measured')
            plt.title(pepname)
            
            save_path= 'E:\\project6_rapid profiling plasma peptides\\figures\\quanti_typical\\'
            newname=pepname.replace('(UniMod:4)','')
            plt.savefig(save_path+"%s.svg"%str(newname),dpi=600, bbox_inches='tight') 
            plt.show() 
    if i==1:
        for proname in namelist:
            sap_pro=df.loc[proname]
            k=sap_pro.to_frame()
            pp=k.iloc[:12,0].to_frame()
            groups=[0.4, 0.4, 0.4, 0.2, 0.2, 0.2,0.1,0.1,0.1,0.05,0.05,0.05]
            pp['groups'] = groups
            plt.rcParams['figure.figsize'] = 5,5
            categories = np.array([ 0, 0, 0, 1, 1, 1, 2, 2,2,3,3,3])
            colormap = np.array(['#5ec962', '#21918c', '#3b528b', '#440154'])
            jitter(pp.groups.values, pp.iloc[:12,0], s=100, c = colormap[categories])
            plt.xlabel(u"\u03bcg/\u03bcL injected")
            plt.ylabel('relative quantity, measured')
            plt.title(proname)
            
            save_path= 'E:\\project6_rapid profiling plasma peptides\\figures\\quanti_typical\\'
            newname=proname.split("|")[-1]
            plt.savefig(save_path+"%s.svg"%str(newname),dpi=600, bbox_inches='tight') 
            plt.show() 




# In[94]:


namelist


# In[97]:


path="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2"

df=plasma_file(dir = path)[1]

df_pval=getpval(df_target=getmean(dftar=df))
wo=df_pval.sort_values(by='slopes',ascending=False)
woai=wo.iloc[:20,:]
namelist=woai.index.tolist()
namelist
sap_pro=df.loc['1/sp|P24593|IBP5_HUMAN']
k=sap_pro.to_frame()
pp=k.iloc[:12,0].to_frame()
groups=[0.4, 0.4, 0.4, 0.2, 0.2, 0.2,0.1,0.1,0.1,0.05,0.05,0.05]
pp['groups'] = groups
plt.rcParams['figure.figsize'] = 5,5
categories = np.array([ 0, 0, 0, 1, 1, 1, 2, 2,2,3,3,3])
colormap = np.array(['#5ec962', '#21918c', '#3b528b', '#440154'])
jitter(pp.groups.values, pp.iloc[:12,0], s=100, c = colormap[categories])
plt.xlabel(u"\u03bcg/\u03bcL injected")
plt.ylabel('relative quantity, measured')
plt.title('P24593|IBP5_HUMAN',fontsize = 26,loc='center')

plt.savefig('E:\\project6_rapid profiling plasma peptides\\figures\\quanti_typical\\sample.svg',dpi=600, bbox_inches='tight') 
plt.show() 



# In[29]:


dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist())
(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq").iloc[:,0].tolist())
dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2").iloc[:,0].tolist())
dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq").iloc[:,0].tolist())
dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq").iloc[:,0].tolist())
dir="E:\\yuming\\2023\\20230118\\csodiaq\\DD")


# In[84]:


typical_quan(path="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",i=0)


# In[194]:


sap_pro=df.loc['RPSEIVIGQC(UniMod:4)K']
k=sap_pro.to_frame()
pp=k.iloc[:12,0].to_frame()
groups=[0.4, 0.4, 0.4, 0.2, 0.2, 0.2,0.1,0.1,0.1,0.05,0.05,0.05]
pp['groups'] = groups
print(pp)
plt.rcParams['figure.figsize'] = 5,5
categories = np.array([ 0, 0, 0, 1, 1, 1, 2, 2,2,3,3,3])

colormap = np.array(['#5ec962', '#21918c', '#3b528b', '#440154'])

jitter(pp.groups.values, pp.iloc[:12,0], s=100, c = colormap[categories])

plt.xlabel(u"\u03bcg/\u03bcL injected")
plt.ylabel('relative quantity, measured')
#plt.title('EYWMDPEGEMKPGR')
#plt.savefig('../plots/240k_all_fragments_EYWMDPEGEMKPGR.svg', bbox_inches='tight')


# In[191]:




plt.rcParams['figure.figsize'] = 5,5
categories = np.array([ 0, 0, 0, 1, 1, 1, 2, 2,2,3,3,3])

colormap = np.array(['#5ec962', '#21918c', '#3b528b', '#440154'])

jitter(pp.groups.values, pp.iloc[:12,0], s=100, c = colormap[categories])

plt.xlabel(u"\u03bcg/\u03bcL injected")
plt.ylabel('relative quantity, measured')
#plt.title('EYWMDPEGEMKPGR')
#plt.savefig('../plots/240k_all_fragments_EYWMDPEGEMKPGR.svg', bbox_inches='tight')


# In[ ]:


# 'np5':  ['RPSEIVIGQC(UniMod:4)K',]


# In[65]:


dfpro_np4.loc['1/sp|P02745|C1QA_HUMAN']


# In[51]:


getmean(dftar=plasma_file(dir = "E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[0])


# In[ ]:




