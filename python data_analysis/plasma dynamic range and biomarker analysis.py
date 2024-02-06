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
def get_all_proteins(dir='E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2'):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    k=[]
    for i in range(len(plasma_pro)):  
        df_pro=pd.read_csv(plasma_pro[i])
        wolist=df_pro['protein'].tolist()
        k.extend(wolist)
        final=remove_dul(tar_list=k)
    return(final)


# In[3]:


#plot_sin_np(title='NP5',figname='pros_np5',listname=procount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist())
#plot_sin_np(title='NP4',figname='pros_np4',listname=procount(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq").iloc[:,0].tolist())
#plot_sin_np(title='NP3',figname='pros_np3',listname=procount(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2").iloc[:,0].tolist())
#plot_sin_np(title='NP2',figname='pros_np2',listname=procount(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq").iloc[:,0].tolist())
#plot_sin_np(title='NP1',figname='pros_np1',listname=procount(dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq").iloc[:,0].tolist())
#plot_sin_np(title='dd',figname='pros_dd',listname=procount(dir="E:\\yuming\\2023\\20230118\\csodiaq\\DD").iloc[:,0].tolist())


# In[4]:


procount(dir='E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2')


# In[5]:


len(get_all_proteins(dir='E:\\yuming\\2023\\20230129\\NP2\\csodiaq'))


# In[6]:


dic_all={'dd': get_all_proteins(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'),
        'np1': get_all_proteins(dir='E:\\yuming\\2023\\20230129\\NP1\\csodiaq'),
        'np2': get_all_proteins(dir='E:\\yuming\\2023\\20230129\\NP2\\csodiaq'),
        'np3': get_all_proteins(dir='E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2'),
        'np4': get_all_proteins(dir='E:\\yuming\\2023\\20230130\\NP4\\csodiaq'),
        'np5': get_all_proteins(dir='E:\\yuming\\2023\\20230130\\NP5\\csodiaq')
        }


# In[7]:


dfdf=pd.DataFrame.from_dict(dic_all,orient='index')


# In[8]:


dfdf


# In[9]:


right=dfdf.transpose()


# In[10]:


right


# In[ ]:





# In[ ]:





# In[9]:


#right.to_excel(r'C:\\Users\\jymbc\\Desktop\\right.xlsx')


# In[11]:


plasma_pro = []
dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq'
for file in os.listdir(dir):
    if file.endswith(".csv") and 'proteinFDR' in file:
        plasma_pro.append(os.path.join(dir, file))
plasma_pro


# In[12]:


k=[]
for i in range(len(plasma_pro)):  
    df_pro=pd.read_csv(plasma_pro[i])
    wolist=df_pro['protein'].tolist()
    k.extend(wolist)
print(len(k))


# In[13]:


len(remove_dul(tar_list=k))


# In[14]:


df_biomarker=pd.read_excel(r'F:\\Seer_project\\FDA_biomarkers.xlsx',sheet_name='Sheet3')

marklist=df_biomarker['biomarkers'].tolist()


# In[15]:


#procount(dir='E:\\yuming\\2022\\20221121\\csodiaq2')

#procount(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq')


# In[16]:


numlist=[0,1,2,3,4,5]


# In[17]:


df_dd=pd.read_csv(plasma_pro[0])
df_dd


# In[14]:


dicte={}

numlist=[0,1,2,3,4,5]
for i in numlist:
    df_dd=pd.read_csv(plasma_pro[i])
    
    alist=[]
    for item in df_dd['protein'].tolist():
        alist.append(item.split('|')[1])
    dicte[str(plasma_pro[i].split('\\')[-1])]=remove_dul(tar_list=alist) 
    dfdfnew=pd.DataFrame.from_dict(dicte, orient='index')             


# In[15]:


youy=dfdfnew.T


# In[16]:


#youy.to_excel(r'C:\\Users\\jymbc\\Desktop\\youy.xlsx')


# In[17]:


youy


# In[18]:


dd_only=remove_dul(youy.iloc[:,0].tolist())
allnps_toge=remove_dul(remove_dul(youy.iloc[:,1].tolist())+remove_dul(youy.iloc[:,2].tolist())+remove_dul(youy.iloc[:,3].tolist())+remove_dul(youy.iloc[:,4].tolist())+remove_dul(youy.iloc[:,5].tolist()))


# In[ ]:





# In[19]:


aa,dd=[],[]

for item in allnps_toge:
    if item in dd_only:
        aa.append(item)
    else:
        dd.append(item)


# In[20]:


common_val=list(set(allnps_toge) & set(dd_only))

newdic={}
newdic['common']=common_val
newdic['unique']=dd
dfdffinal=pd.DataFrame.from_dict(newdic, orient='index')


# In[21]:


yonique=dfdffinal.T
#yonique.to_excel(r'C:\\Users\\jymbc\\Desktop\\yonique.xlsx')


# In[22]:


df_pro=pd.read_csv(plasma_pro[i])
dfdfpro=df_pro.sort_values('ionCount', ascending=False).drop_duplicates(['protein'])
aa=[]


# In[22]:


def get_scatter(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",i=1):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    plasma_pro
    df_pro=pd.read_csv(plasma_pro[i])
    dfdfpro=df_pro.sort_values('ionCount', ascending=False).drop_duplicates(['protein'])
    aa=[]
    for item in dfdfpro['protein'].tolist():
        aa.append(item.split('|')[1])
    newlist=[]
    for value in dfdfpro['ionCount'].tolist():
        newlist.append(math.log10(value))
    dfdfpro['pro_id']=aa
    dfdfpro['inten']=newlist
    af=sorted(newlist)
    newdf=dfdfpro.sort_values(by='inten')
    newdf['rank']=list(range(len(newlist)))
    appended_data = []
    for name in marklist:
        data=newdf[newdf['pro_id'] == name]
        appended_data.append(data)
    appended_data = pd.concat(appended_data)
    return(appended_data,af)
def plot_scatter(figname='NP3',path="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",num=1):
    appended_data,af=get_scatter(dir=path,i=num)
    matplotlib.rcParams['font.family'] = "Arial"
    plt.rcParams['figure.figsize']=(3,9)
    sns.set_style("whitegrid")
    ax = plt.subplot()
    TK = plt.gca() #获取边框
    matplotlib.rcParams['font.family'] = "Arial"
    TK.spines['right'].set_visible(False)
    TK.spines['top'].set_visible(False)
    TK.spines['bottom'].set_visible(False)
    TK.spines['left'].set_visible(True)
    plt.yticks(np.arange(0,10, 1))
    plt.ylim(3,8)
    plt.xlabel(figname,fontsize=28)
    #plt.ylabel('Intensity [Log10]',fontsize=24)
    plt.scatter(range(len(af)),af,marker='o',s=35,linewidths=1,facecolors='none', edgecolors='grey')
    plt.scatter(appended_data['rank'].tolist(),appended_data['inten'].tolist(),marker='o',s=28,linewidths=0.5,c='r', edgecolors='red')
    plt.yticks(fontsize=28)
    ax.set_xticks([])
    #ax.set_yticks([])
    fig_path= 'F:\\Seer_project\\figures\\biomarker_scatter\\'
    plt.savefig(fig_path + "%s.svg" % figname,dpi=800, bbox_inches='tight')  
    plt.show()
def count_overlap(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",i=1):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    plasma_pro
    df_pro=pd.read_csv(plasma_pro[i])
    dfdfpro=df_pro.sort_values('ionCount', ascending=False).drop_duplicates(['protein'])
    aa=[]
    for item in dfdfpro['protein'].tolist():
        aa.append(item.split('|')[1])
    newlist=[]
    for value in dfdfpro['ionCount'].tolist():
        newlist.append(math.log10(value))
    dfdfpro['pro_id']=aa
    dfdfpro['inten']=newlist
    af=sorted(newlist)
    newdf=dfdfpro.sort_values(by='inten')
    newdf['rank']=list(range(len(newlist)))
    appended_data = []
    for name in marklist:
        data=newdf[newdf['pro_id'] == name]
        appended_data.append(data)
    appended_data = pd.concat(appended_data)
    #return(len(appended_data['inten'].tolist()))    
    return(appended_data)
    


# In[26]:


count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=0)


# In[48]:


ji=[]
for item in count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=2)['leadingProtein'].tolist():
    ji.append(item.split('|')[1])
len(remove_dul(tar_list=ji))


# In[43]:


alll=[]
for k in range(6):
    ji=[]
    for item in count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=k)['leadingProtein'].tolist():
        ji.append(item.split('|')[1])
    alll.append(ji)


# In[47]:


len(remove_dul(tar_list=[item for sublist in alll for item in sublist]))


# In[21]:


plot_scatter(figname='DD',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=0)

#count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=0)


# In[34]:


plot_scatter(figname='NP1',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=1)
count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=1)


# In[35]:


plot_scatter(figname='NP2',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=2)
count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=2)


# In[36]:


plot_scatter(figname='NP3',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=3)
count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=3)


# In[37]:


plot_scatter(figname='NP4',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=4)
count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=4)


# In[38]:


plot_scatter(figname='NP5',path='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',num=5)
count_overlap(dir='E:\\yuming\\2022\\20221121\\New folder\\csodiaq',i=5)


# In[27]:


plot_scatter(figname='NP3',path="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",num=1)

count_overlap(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2",i=1)


# In[ ]:





# In[28]:


plot_scatter(figname='DD',path='E:\\yuming\\2023\\20230118\\csodiaq\\DD',num=7)
count_overlap(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD',i=7)


# In[29]:


plot_scatter(figname='NP1',path='E:\\yuming\\2023\\20230129\\NP1\\csodiaq',num=1)
count_overlap(dir='E:\\yuming\\2023\\20230129\\NP1\\csodiaq',i=1)


# In[30]:


plot_scatter(figname='NP2',path='E:\\yuming\\2023\\20230129\\NP2\\csodiaq',num=1)

count_overlap(dir='E:\\yuming\\2023\\20230129\\NP2\\csodiaq',i=1)


# In[31]:


plot_scatter(figname='NP4',path='E:\\yuming\\2023\\20230130\\NP4\\csodiaq',num=2)

count_overlap(dir='E:\\yuming\\2023\\20230130\\NP4\\csodiaq',i=2)


# In[32]:


plot_scatter(figname='NP5',path='E:\\yuming\\2023\\20230130\\NP5\\csodiaq',num=0)

count_overlap(dir='E:\\yuming\\2023\\20230130\\NP5\\csodiaq',i=0)


# In[ ]:





# In[ ]:





# In[ ]:




