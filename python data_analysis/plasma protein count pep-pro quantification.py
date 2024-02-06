#!/usr/bin/env python
# coding: utf-8

# In[9]:


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


# In[10]:


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

def count(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'):
    ['''count the protein ID from a .csv file''']
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    pro_number=[]
    for i in range(len(plasma_pro)):
        df=pd.read_csv(plasma_pro[i]) 
        pro_number.append(len(remove_dul(tar_list=df['protein'].tolist())))
    pro_number 
    return(pro_number) 


def plot_sin_np(title='NP5',figname='pros_np5',listname=procount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist()):
    prolist_np=[average(listname[:3]),average(listname[3:6]),average(listname[6:9]),average(listname[9:12])]
    np_std=[np.std(listname[:3]),np.std(listname[3:6]),np.std(listname[6:9]),np.std(listname[9:12])]
    name2=['1','0.5','0.25','0.125']
    sns.set_theme(style="ticks")
    plt.figure(figsize=(6,4),edgecolor="#04253a")
    matplotlib.rcParams['font.family'] = "Arial"
    bar_width=0.5
    plt.bar(x=name2,height=prolist_np,yerr=np_std,color='firebrick',alpha=0.8,width=bar_width,capsize=5,label=title)
    #plt.errorbar(name2,prolist_np5,yerr=np5_std,color='firebrick',label='NP5', marker='o',ms = 5,linewidth=2,capsize=6,linestyle='dotted')
    plt.xlabel('Concentration',fontsize=24)
    plt.ylabel('Number of proteins',fontsize=24)
    plt.title(title,fontsize=24)
    plt.xticks(fontsize=20)
    #plt.xticks(rotation=20)
    plt.yticks(fontsize=22)
    #plt.xticks(np.arange(0, 24, 4))
    plt.yticks(np.arange(0,400, 100))
    #plt.xlim(7.8,23)
    plt.ylim(0,400)
    plt.legend(fontsize=20)
    plt.grid()
    plt.subplots_adjust(left=0.25,bottom=0.2)
    #plt.tight_layout()
    #plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
    #plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\robust and reproducibility\\CV analysis_protein.svg', dpi=600) 
    fig_path= 'E:\\project6 rapid profiling plasma peptides\\figures\\'
    plt.savefig(fig_path + "%s.svg" % figname,dpi=600, bbox_inches='tight')  
    plt.show()

def plot_pep_np(title='NP5',figname='pros_np5',listname=pepcount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist()):
    prolist_np=[average(listname[:3]),average(listname[3:6]),average(listname[6:9]),average(listname[9:12])]
    np_std=[np.std(listname[:3]),np.std(listname[3:6]),np.std(listname[6:9]),np.std(listname[9:12])]
    name2=['1','0.5','0.25','0.125']
    sns.set_theme(style="ticks")
    plt.figure(figsize=(6,4),edgecolor="#04253a")
    matplotlib.rcParams['font.family'] = "Arial"
    bar_width=0.5
    plt.bar(x=name2,height=prolist_np,yerr=np_std,color='firebrick',alpha=0.8,width=bar_width,capsize=5,label=title)
    #plt.errorbar(name2,prolist_np5,yerr=np5_std,color='firebrick',label='NP5', marker='o',ms = 5,linewidth=2,capsize=6,linestyle='dotted')
    plt.xlabel('Concentration',fontsize=24)
    plt.ylabel('Number of proteins',fontsize=24)
    plt.title(title,fontsize=24)
    plt.xticks(fontsize=20)
    #plt.xticks(rotation=20)
    plt.yticks(fontsize=22)
    #plt.xticks(np.arange(0, 24, 4))
    plt.yticks(np.arange(0,1500, 500))
    #plt.xlim(7.8,23)
    plt.ylim(0,1500)
    plt.legend(fontsize=20)
    plt.grid()
    plt.subplots_adjust(left=0.25,bottom=0.2)
    #plt.tight_layout()
    #plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
    #plt.savefig('D:\\Project 2 DIrect infusion shotgun proteome analysis\\robust and reproducibility\\CV analysis_protein.svg', dpi=600) 
    fig_path= 'E:\\project6 rapid profiling plasma peptides\\figures\\'
    plt.savefig(fig_path + "%s.svg" % figname,dpi=600, bbox_inches='tight')  
    plt.show()


# In[21]:


#procount(dir="E:\\yuming\\2023\\20230226\\csodiaq")


# In[20]:


#procount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")


# In[11]:


# targeted results among different NPs

listnp5=procount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist()
listnp4=procount(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq").iloc[:,0].tolist()
listnp3=procount(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2").iloc[:,0].tolist()
listnp2=procount(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq").iloc[:,0].tolist()
listnp1=procount(dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq").iloc[:,0].tolist()
listdd=procount(dir="E:\\yuming\\2023\\20230118\\csodiaq\\DD").iloc[:,0].tolist()


# In[14]:


listnp4


# In[19]:



def prolist(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD'):
    plasma_pro = []
    for file in os.listdir(dir):
        if file.endswith(".csv") and 'proteinFDR' in file:
            plasma_pro.append(os.path.join(dir, file))
    dic={}
    for i in range(len(plasma_pro)):
        df=pd.read_csv(plasma_pro[i]) 
        yummy=remove_dul(tar_list=df['protein'].tolist())
        dic[str(i)]=yummy
    return(dic)


# In[26]:


prolist(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")['0']


# In[32]:


namelist=['0','1','2','3','4','5','6']

haha=[]

for nameq in namelist:

    pron_np5=prolist(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq")[nameq]
    pron_np4=prolist(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq")[nameq]
    pron_np3=prolist(dir="E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2")[nameq]
    pron_np2=prolist(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq")[nameq]
    pron_np1=prolist(dir="E:\\yuming\\2023\\20230129\\NP1\\csodiaq")[nameq]
    pron_dd=prolist(dir="E:\\yuming\\2023\\20230118\\csodiaq\\DD")[nameq]


    final_count= len(remove_dul(pron_np5+pron_np4+pron_np3+pron_np2+pron_np1))
    haha.append(final_count)
haha


# In[ ]:





# In[ ]:





# In[ ]:





# In[4]:


name=['DD','NP1','NP2','NP3','NP4','NP5']


quan_error2=[np.std(listdd[:3]),np.std(listnp1[:3]),np.std(listnp2[:3]),np.std(listnp3[:3]),np.std(listnp4[:3]),np.std(listnp5[:3])]

pro_quan2=[average(listdd[:3]),average(listnp1[:3]),average(listnp2[:3]),average(listnp3[:3]),average(listnp4[:3]),average(listnp5[:3])]


# In[5]:


pro_quan2


# In[5]:


sns.set_theme(style="ticks")
plt.figure(figsize=(9,6),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.55
plt.bar(x=name,height=pro_quan2,yerr=quan_error2,color='firebrick',alpha=0.8,width=bar_width,capsize=5)
#plt.xlabel('Treatment',fontsize=24)
plt.ylabel('Number of proteins',fontsize=24)
plt.title('DISPA',fontsize=26)
plt.xticks(fontsize=22)
#plt.xticks(rotation=20)
plt.yticks(fontsize=22)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,400, 100))
#plt.xlim(7.8,23)
plt.ylim(0,400)
#plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
#plt.savefig('D:\\project6_rapid profiling plasma peptides\\figures\\Pro_quan_targeted.svg', dpi=800) 
plt.show()


# In[6]:


procount(dir='E:\\yuming\\2023\\20230214\\csodiaq')


# In[7]:


procount(dir='E:\\yuming\\2023\\20230116\\csodiaq')


# In[8]:


procount(dir='E:\\yuming\\2023\\20230117\\csodiaq')


# In[9]:


procount(dir='E:\\yuming\\2022\\20221121\\csodiaq2')


# In[53]:


NP1=[273,307,266]
NP2=[393,353,320]
NP3=[341,346,319]
NP4=[348,322,313]
NP5=[352,361,298]
DD=[254,262,245]
all_NPs=[459, 466, 456]
         
name=['NP1','NP2','nNP3','NP4','NP5','all_NPs','DD']

quan_error=[np.std(NP1),np.std(NP2),np.std(NP3),np.std(NP4),np.std(NP5),np.std(all_NPs),np.std(DD)]

pro_quan=[average(NP1),average(NP2),average(NP3),average(NP4),average(NP5),average(all_NPs),average(DD)]


# In[54]:


# protein identifications from different NPs  ---scouting experiments

sns.set_theme(style="ticks")
plt.figure(figsize=(9.5,6),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.55
plt.bar(x=name,height=pro_quan,yerr=quan_error,color='grey',alpha=0.7,width=bar_width,capsize=5)

#plt.xlabel('Treatment',fontsize=24)
plt.ylabel('Number of proteins',fontsize=26)
#plt.title('DISPA',fontsize=26)
#plt.xticks(fontsize=22)
#plt.xticks(rotation=20)
plt.xticks(name,fontsize=22,rotation=20)
plt.yticks(fontsize=24)
#plt.xticks(np.arange(0, 24, 4))

plt.yticks(np.arange(0,500, 100))
#plt.xlim(7.8,23)
plt.ylim(0,500)
#plt.legend(fontsize=20,frameon=False)
plt.grid()
#plt.gca().invert_yaxis()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('F:\\Seer_project\\figures\\Pro_quan_scout2.svg', dpi=800) 
plt.show()


# In[36]:


# protein IDS from DDA LC-MS/MS
proID_DDA=[269, 302, 421, 403, 418, 437]

name=['DD','NP1','NP2','NP3','NP4','NP5']

sns.set_theme(style="ticks")
plt.figure(figsize=(9,6),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.55
plt.bar(x=name,height=proID_DDA,color='firebrick',alpha=0.8,width=bar_width)
#plt.xlabel('treatment',fontsize=24)
plt.ylabel('Number of proteins',fontsize=26)
#plt.title('DDA_LC_MS/MS',fontsize=26)
plt.xticks(fontsize=22)
#plt.xticks(rotation=20)
plt.yticks(fontsize=24)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,500, 100))
#plt.xlim(7.8,23)
plt.ylim(0,500)
plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
plt.gca().invert_yaxis()
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('F:\\Seer_project\\figures\\Pro_quan_LC.svg', dpi=800) 
plt.show()


# In[6]:


pep_NP1=[1788,1917,1929]
pep_NP2=[2042,2074,2133]
pep_NP3=[1976,2315,1965]
pep_NP4=[2029,2202,1863]
pep_NP5=[2183,2426,2430]
pep_DD=[1760,1846,1821]


pep_quan_error=[np.std(pep_DD),np.std(pep_NP1),np.std(pep_NP2),np.std(pep_NP3),np.std(pep_NP4),np.std(pep_NP5)]


pep_quan=[average(pep_DD),average(pep_NP1),average(pep_NP2),average(pep_NP3),average(pep_NP4),average(pep_NP5)]


# In[7]:


pep_quan


# In[14]:


# peptide identifications from different NPs  ---scouting experiments

sns.set_theme(style="ticks")
plt.figure(figsize=(9,5.5),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
bar_width=0.55
plt.bar(x=name,height=pep_quan,yerr=pep_quan_error,color='purple',alpha=0.7,width=bar_width,capsize=5)
#plt.xlabel('Treatment',fontsize=24)
plt.ylabel('Number of peptides',fontsize=24)
#plt.title('DISPA',fontsize=26)
plt.xticks(fontsize=22)
#plt.xticks(rotation=20)
plt.yticks(fontsize=22)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,2800, 500))
#plt.xlim(7.8,23)
plt.ylim(0,2800)
#plt.legend(fontsize=20,frameon=False)
plt.grid()
plt.subplots_adjust(left=0.25,bottom=0.2)
#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('F:\\Seer_project\\figures\\pep_quan_scout.svg', dpi=800) 
plt.show()


# In[15]:


procount(dir='E:\\yuming\\2023\\20230214\\csodiaq')


# In[16]:


procount(dir='E:\\yuming\\2023\\20230116\\csodiaq')


# In[27]:


procount(dir='E:\\yuming\\2022\\20221121\\new folder\\csodiaq')


# In[8]:



prolist_np1=procount(dir='E:\\yuming\\2023\\20230129\\NP1\\csodiaq').iloc[:,0].tolist()  #NP1    

prolist_np2=procount(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq").iloc[:,0].tolist()    #NP2

prolist_np3=procount(dir='E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2').iloc[:,0].tolist()  #NP3  

prolist_np4=procount(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq").iloc[:,0].tolist()   # NP4

prolist_np5=procount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist()   # NP5

prolist_dd=procount(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD').iloc[:,0].tolist()    # DD

peplist_np1=pepcount(dir='E:\\yuming\\2023\\20230129\\NP1\\csodiaq').iloc[:,0].tolist()  #NP1 peptides    

peplist_np2=pepcount(dir="E:\\yuming\\2023\\20230129\\NP2\\csodiaq").iloc[:,0].tolist()    #NP2 peptides

peplist_np3=pepcount(dir='E:\\yuming\\2023\\20230129\\NP3\\New folder\\csodiaq2').iloc[:,0].tolist()  #NP3 peptides 

peplist_np4=pepcount(dir="E:\\yuming\\2023\\20230130\\NP4\\csodiaq").iloc[:,0].tolist()   # NP4 peptides

peplist_np5=pepcount(dir="E:\\yuming\\2023\\20230130\\NP5\\csodiaq").iloc[:,0].tolist()   # NP5 peptides

peplist_dd=pepcount(dir='E:\\yuming\\2023\\20230118\\csodiaq\\DD').iloc[:,0].tolist()    # DD peptides


# In[9]:


prolist=[average(prolist_dd[3:6]),average(prolist_np1[:3]),average(prolist_np2[:3]),average(prolist_np3[:3]),average(prolist_np4[:3]),average(prolist_np5[:3])]

pro_error=[np.std(prolist_dd[3:6]),np.std(prolist_np1[:3]),np.std(prolist_np2[:3]),np.std(prolist_np3[:3]),np.std(prolist_np4[:3]),np.std(prolist_np5[:3])]

peplist=[average(peplist_dd[3:6]),average(peplist_np1[:3]),average(peplist_np2[:3]),average(peplist_np3[:3]),average(peplist_np4[:3]),average(peplist_np5[:3])]

pep_error=[np.std(peplist_dd[3:6]),np.std(peplist_np1[:3]),np.std(peplist_np2[:3]),np.std(peplist_np3[:3]),np.std(peplist_np4[:3]),np.std(peplist_np5[:3])]



# In[12]:


peplist


# In[209]:


ug1=[average(prolist_dd[0:3]),average(prolist_np1[:3]),average(prolist_np2[:3]),average(prolist_np3[:3]),average(prolist_np4[:3]),average(prolist_np5[:3])]
ug2=[average(prolist_dd[3:6]),average(prolist_np1[3:6]),average(prolist_np2[3:6]),average(prolist_np3[3:6]),average(prolist_np4[3:6]),average(prolist_np5[3:6])]
ug3=[average(prolist_dd[6:9]),average(prolist_np1[6:9]),average(prolist_np2[6:9]),average(prolist_np3[6:9]),average(prolist_np4[6:9]),average(prolist_np5[6:9])]
ug4=[average(prolist_dd[9:12]),average(prolist_np1[9:12]),average(prolist_np2[9:12]),average(prolist_np3[9:12]),average(prolist_np4[9:12]),average(prolist_np5[9:12])]



ug1_std=[np.std(prolist_dd[0:3]),np.std(prolist_np1[:3]),np.std(prolist_np2[:3]),np.std(prolist_np3[:3]),np.std(prolist_np4[:3]),np.std(prolist_np5[:3])]
ug2_std=[np.std(prolist_dd[3:6]),np.std(prolist_np1[3:6]),np.std(prolist_np2[3:6]),np.std(prolist_np3[3:6]),np.std(prolist_np4[3:6]),np.std(prolist_np5[3:6])]
ug3_std=[np.std(prolist_dd[6:9]),np.std(prolist_np1[6:9]),np.std(prolist_np2[6:9]),np.std(prolist_np3[6:9]),np.std(prolist_np4[6:9]),np.std(prolist_np5[6:9])]
ug4_std=[np.std(prolist_dd[9:12]),np.std(prolist_np1[9:12]),np.std(prolist_np2[9:12]),np.std(prolist_np3[9:12]),np.std(prolist_np4[9:12]),np.std(prolist_np5[9:12])]







# In[265]:


pepug1=[average(peplist_dd[0:3]),average(peplist_np1[:3]),average(peplist_np2[:3]),average(peplist_np3[:3]),average(peplist_np4[:3]),average(peplist_np5[:3])]
pepug2=[average(peplist_dd[3:6]),average(peplist_np1[3:6]),average(peplist_np2[3:6]),average(peplist_np3[3:6]),average(peplist_np4[3:6]),average(peplist_np5[3:6])]
pepug3=[average(peplist_dd[6:9]),average(peplist_np1[6:9]),average(peplist_np2[6:9]),average(peplist_np3[6:9]),average(peplist_np4[6:9]),average(peplist_np5[6:9])]
pepug4=[average(peplist_dd[9:12]),average(peplist_np1[9:12]),average(peplist_np2[9:12]),average(peplist_np3[9:12]),average(peplist_np4[9:12]),average(peplist_np5[9:12])]



pepug1_std=[np.std(peplist_dd[0:3]),np.std(peplist_np1[:3]),np.std(peplist_np2[:3]),np.std(peplist_np3[:3]),np.std(peplist_np4[:3]),np.std(peplist_np5[:3])]
pepug2_std=[np.std(peplist_dd[3:6]),np.std(peplist_np1[3:6]),np.std(peplist_np2[3:6]),np.std(peplist_np3[3:6]),np.std(peplist_np4[3:6]),np.std(peplist_np5[3:6])]
pepug3_std=[np.std(peplist_dd[6:9]),np.std(peplist_np1[6:9]),np.std(peplist_np2[6:9]),np.std(peplist_np3[6:9]),np.std(peplist_np4[6:9]),np.std(peplist_np5[6:9])]
pepug4_std=[np.std(peplist_dd[9:12]),np.std(peplist_np1[9:12]),np.std(peplist_np2[9:12]),np.std(peplist_np3[9:12]),np.std(peplist_np4[9:12]),np.std(peplist_np5[9:12])]


# In[ ]:





# In[ ]:





# In[259]:


import matplotlib
sns.set_theme(style="ticks")
plt.figure(figsize=(12,6),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
listname1=['DD','NP1','NP2','NP3','NP4','NP5']
bar_width=0.2
plt.bar(x=range(len(listname1)),height=ug1,yerr=ug1_std, color='tab:cyan',label='0.4μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+bar_width,height=ug2,yerr=ug2_std,color='blue',label='0.2μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+2*bar_width,height=ug3,yerr=ug3_std, color='purple',label='0.1μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+3*bar_width,height=ug4,yerr=ug4_std,  color='darkgray',label='0.05μg/μl',alpha=0.8,width=bar_width,capsize=5)
#plt.bar(x=np.arange(len(listname1))+4*bar_width,height=np5,yerr=np5_std,  color='darkgray',label='240k',alpha=0.8,width=bar_width)

plt.xticks(np.arange(len(listname1))+0.31, listname1)
#plt.xlabel('Concentration, μg/µl',fontsize=24)
plt.ylabel('Number of proteins',fontsize=26) 
#plt.title('10 peptide quantities',fontsize=26)
plt.xticks(fontsize=26)
#plt.xticks(rotation=60)
plt.yticks(fontsize=26)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,300, 100))
#plt.xlim(7.8,23)
plt.ylim(0,350)
plt.legend(fontsize=20,loc='upper center', bbox_to_anchor=(0.5, 1.16),
          ncol=4, fancybox=True, shadow=True)
plt.grid()

#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('E:\\project6 rapid profiling plasma peptides\\figures\\Pro_quan_all_target.svg', dpi=800) 
plt.show()


# In[269]:


import matplotlib
sns.set_theme(style="ticks")
plt.figure(figsize=(12,6),edgecolor="#04253a")
matplotlib.rcParams['font.family'] = "Arial"
listname1=['DD','NP1','NP2','NP3','NP4','NP5']
bar_width=0.2
plt.bar(x=range(len(listname1)),height=pepug1,yerr=pepug1_std, color='tab:cyan',label='0.4μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+bar_width,height=pepug2,yerr=pepug2_std,color='blue',label='0.2μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+2*bar_width,height=pepug3,yerr=pepug3_std, color='purple',label='0.1μg/μl',alpha=0.8,width=bar_width,capsize=5)
plt.bar(x=np.arange(len(listname1))+3*bar_width,height=pepug4,yerr=pepug4_std,  color='darkgray',label='0.05μg/μl',alpha=0.8,width=bar_width,capsize=5)
#plt.bar(x=np.arange(len(listname1))+4*bar_width,height=np5,yerr=np5_std,  color='darkgray',label='240k',alpha=0.8,width=bar_width)

plt.xticks(np.arange(len(listname1))+0.31, listname1)
#plt.xlabel('Concentration, μg/µl',fontsize=24)
plt.ylabel('Number of peptides',fontsize=26)
#plt.title('10 peptide quantities',fontsize=26)
plt.xticks(fontsize=26)
#plt.xticks(rotation=60)
plt.yticks(fontsize=26)
#plt.xticks(np.arange(0, 24, 4))
plt.yticks(np.arange(0,1300, 500))
#plt.xlim(7.8,23)
plt.ylim(0,1300)
plt.legend(fontsize=20,loc='upper center', bbox_to_anchor=(0.5, 1.16),
          ncol=4, fancybox=True, shadow=True)
plt.grid()

#plt.tight_layout()
#plt.legend(bbox_to_anchor=(1.02,1), borderaxespad=0,frameon=False)
plt.savefig('E:\\project6 rapid profiling plasma peptides\\figures\\Pep_quan_all_target.svg', dpi=800) 
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




