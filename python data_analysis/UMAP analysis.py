#!/usr/bin/env python
# coding: utf-8

# In[50]:


import umap
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from sklearn.datasets import load_iris
import pandas as pd
import numpy as np


# # Doing UMAP analysis for the scouting results from all NPs and neat plasma  

# In[3]:


dfdf=pd.read_excel(r'F:\\Seer_project\\NPs_scout.xlsx',sheet_name='Sheet1')
dfdf.set_index('Unnamed: 0',inplace=True)


# In[4]:


dfdf=dfdf.fillna(1)
dfdf.index=range(1,622)


# In[5]:


dfdf


# In[6]:


i=1

np_only=dfdf[(dfdf['DD']==i)&((dfdf['NP1']!=i)|(dfdf['NP2']!=i)|(dfdf['NP3']!=i)|(dfdf['NP4']!=i)|(dfdf['NP5']!=i))]

dd_only=dfdf[(dfdf['DD']!=i)&((dfdf['NP1']==i)&(dfdf['NP2']==i)&(dfdf['NP3']==i)&(dfdf['NP4']==i)&(dfdf['NP5']==i))]

share=dfdf[(dfdf['DD']!=i)&((dfdf['NP1']!=i)|(dfdf['NP2']!=i)|(dfdf['NP3']!=i)|(dfdf['NP4']!=i)|(dfdf['NP5']!=i))]


# In[7]:


np=dfdf[((dfdf['NP1']!=i)|(dfdf['NP2']!=i)|(dfdf['NP3']!=i)|(dfdf['NP4']!=i)|(dfdf['NP5']!=i))]

dd=dfdf[(dfdf['DD']!=i)]


# In[8]:


# Mark rows based on their index names being in the three lists
list1,list2=np.index,dd.index


def mark_rows(idx):
    if idx in list1:
        return "yes"
    else:
        return "no"
def mark_rows2(idx):
    if idx in list2:
        return "yes"
    else:
        return "no"

dfdf['np'] = dfdf.index.map(mark_rows)
dfdf['dd'] = dfdf.index.map(mark_rows2)   


# In[9]:


dfdf


# In[10]:


final=dfdf.iloc[:,:6]


# In[11]:



reducer = umap.UMAP(random_state=42,n_neighbors=10, min_dist=1)

embedding = reducer.fit_transform(final)


# In[12]:


df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'])
df['np'] = dfdf['np'].tolist()
df['dd'] = dfdf['dd'].tolist()

df['biomarkers'] = dfdf['Biomarkers'].tolist()


# In[13]:


df


# In[14]:


knp=df[df['np']=='yes']
kdd=df[df['dd']=='yes']
kbio=df[df['biomarkers']=='yes']


# In[ ]:





# In[152]:


#DC4A51
matplotlib.rcParams['font.family'] = "Arial"
plt.figure(figsize=(6, 5.5))
sns.set_style("whitegrid")
plt.scatter(knp['UMAP1'], knp['UMAP2'], c='#428DDD', s=90, marker=".", alpha=0.7,label='NPs')
plt.scatter(kdd['UMAP1'], kdd['UMAP2'],  c='#F3DB3D', s=90,marker='.', alpha=0.6,label='DD')
plt.scatter(kbio['UMAP1'], kbio['UMAP2'],  c='#DC4A51', s=70,marker='^', alpha=0.7,label='Biomarker')
plt.legend(fontsize=16)
plt.xlabel('UMAP 2',fontsize=20)
plt.ylabel('UMAP 1',fontsize=20)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
#plt.legend(loc='upper left', bbox_to_anchor=(1, 1),fontsize=16)
plt.savefig(r'F:\\Seer_project\\figures\\scatter_biomarker_NP_scout.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[58]:


dfdf.iloc[:,:6]


# In[57]:


# do reverse log10 to the dataframe


df_ratio=np.power(10,dfdf.iloc[:,:6])
df_ratio


# In[64]:



ratio=[]
for i in range(len(df_ratio)):
    bilv=np.log10(max(df_ratio.iloc[i,1:6])/df_ratio.iloc[i,0])
    ratio.append(bilv)
len(ratio)


# In[61]:


len([i for i in ratio if i>0])


# In[65]:


df_ratio['ratio']=ratio


# In[36]:


#dfdf.to_excel(r'F:\\Seer_project\\scout_all_log.xlsx')


# In[78]:


#ratio


# In[153]:


from matplotlib.colors import TwoSlopeNorm

color_gradient = ratio
matplotlib.rcParams['font.family'] = "Arial"
plt.figure(figsize=(6, 5.5))
sns.set_style("whitegrid")

# Set 0 as the center for the colormap
norm = TwoSlopeNorm(vcenter=0, vmin=np.array(color_gradient).min(), vmax=np.array(color_gradient).max())

scatter=plt.scatter(df['UMAP1'], df['UMAP2'], c=color_gradient, s=70, marker=".",  cmap="viridis_r", norm=norm)
#twilight_shifted
#
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel('UMAP 2',fontsize=20)
plt.ylabel('UMAP 1',fontsize=20)
# Define position and size for colorbar axes
cax_position = [0.25, 0.31, 0.04, 0.2]  # Modify these values as needed
cax = plt.gcf().add_axes(cax_position)

# Draw the colorbar in the specified axes
cbar = plt.colorbar(scatter, cax=cax)
cbar.set_ticks([-3, 0, 4])
#cbar.set_label('Ratio', rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=16)

plt.savefig(r'F:\\Seer_project\\figures\\scatter_biomarker_NP_ratio_scout.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[19]:


# doing umap analysis for five different NPs
or_df=dfdf.iloc[:,1:6]
or_df


# In[ ]:





# # Doing UMAP analysis for the targeted results(1.4min) from all NPs and neat plasma

# In[ ]:





# In[121]:


da=pd.read_excel(r'F:\\Seer_project\\NPs_DD_Highcon_only.xlsx')


# In[122]:


da


# In[123]:


da.set_index('Unnamed: 0',inplace=True)


# In[124]:


da=da.fillna(1)
da.index=range(1,444)


# In[125]:


da


# In[126]:


i=1

np_only_tar=da[(da['DD']==i)&((da['NP1']!=i)|(da['NP2']!=i)|(da['NP3']!=i)|(da['NP4']!=i)|(da['NP5']!=i))]

dd_only_tar=da[(da['DD']!=i)&((da['NP1']==i)&(da['NP2']==i)&(da['NP3']==i)&(da['NP4']==i)&(da['NP5']==i))]

share_tar=da[(da['DD']!=i)&((da['NP1']!=i)|(da['NP2']!=i)|(da['NP3']!=i)|(da['NP4']!=i)|(da['NP5']!=i))]


# In[127]:


np_tar=da[((da['NP1']!=i)|(da['NP2']!=i)|(da['NP3']!=i)|(da['NP4']!=i)|(da['NP5']!=i))]

dd_tar=da[(da['DD']!=i)]


# In[128]:


# Mark rows based on their index names being in the three lists
list1,list2=np_tar.index,dd_tar.index

def mark_rows(idx):
    if idx in list1:
        return "yes"
    else:
        return "no"
def mark_rows2(idx):
    if idx in list2:
        return "yes"
    else:
        return "no"

da['np'] =da.index.map(mark_rows)
da['dd'] = da.index.map(mark_rows2)   


# In[43]:


#da.to_csv(r'F:\\Seer_project\\NPs_DD_Highcon_only.csv')


# In[129]:


final_tar=da.iloc[:,:6]    # make a final_Tar for umap fit 


# In[130]:


final_tar


# In[131]:


reducer_tar = umap.UMAP(random_state=42,n_neighbors=10, min_dist=1)

embedding_tar = reducer_tar.fit_transform(final_tar)


# In[132]:


dada = pd.DataFrame(embedding_tar, columns=['UMAP1', 'UMAP2'])
dada['np'] = da['np'].tolist()
dada['dd'] = da['dd'].tolist()

dada['biomarkers'] = da['Biomarkers'].tolist()


# In[133]:


dada


# In[134]:


knp_tar=dada[dada['np']=='yes']
kdd_tar=dada[dada['dd']=='yes']
kbio_tar=dada[dada['biomarkers']=='yes']


# In[154]:


matplotlib.rcParams['font.family'] = "Arial"
plt.figure(figsize=(6, 5.5))
sns.set_style("whitegrid")
plt.scatter(knp_tar['UMAP1'], knp_tar['UMAP2'], c='#428DDD', s=90, marker=".", alpha=0.75,label='NPs')
plt.scatter(kdd_tar['UMAP1'], kdd_tar['UMAP2'],  c='#F3DB3D', s=90,marker='.', alpha=0.6,label='DD')
plt.scatter(kbio_tar['UMAP1'], kbio_tar['UMAP2'],  c='#DC4A51', s=70,marker='^', alpha=0.7,label='Biomarker')
#plt.legend(fontsize=16)
plt.xlabel('UMAP 2',fontsize=20)
plt.ylabel('UMAP 1',fontsize=20)
plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
#plt.legend(loc='upper left', bbox_to_anchor=(1, 1),fontsize=16)
plt.legend(fontsize=17)
plt.savefig(r'F:\\Seer_project\\figures\\scatter_biomarker_NP_tar.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[147]:



df_ratio2=np.power(10,da.iloc[:,:6])
ratio_tar=[]
for i in range(len(df_ratio2)):
    bilv2=np.log10(max(df_ratio2.iloc[i,1:6])/df_ratio2.iloc[i,0])
    ratio_tar.append(bilv2)
len(ratio_tar)


# In[148]:


da['ratio']=ratio_tar


# In[149]:


da


# In[155]:


from matplotlib.colors import TwoSlopeNorm

color_gradient_tar = ratio_tar
matplotlib.rcParams['font.family'] = "Arial"
plt.figure(figsize=(6, 5.5))
sns.set_style("whitegrid")

# Set 0 as the center for the colormap
norm = TwoSlopeNorm(vcenter=0, vmin=np.array(color_gradient_tar).min(), vmax=np.array(color_gradient_tar).max())
scatter=plt.scatter(dada['UMAP1'], dada['UMAP2'], c=color_gradient_tar, s=70, marker=".",  cmap="viridis", norm=norm)

plt.yticks(fontsize=18)
plt.xticks(fontsize=18)
plt.xlabel('UMAP 2',fontsize=20)
plt.ylabel('UMAP 1',fontsize=20)
# Define position and size for colorbar axes
cax_position = [0.25, 0.29, 0.04, 0.2]  # Modify these values as needed
cax = plt.gcf().add_axes(cax_position)

# Draw the colorbar in the specified axes
cbar = plt.colorbar(scatter, cax=cax)
cbar.set_ticks([-4,0,5])
#cbar.set_label('Ratio', rotation=270, labelpad=15)
cbar.ax.tick_params(labelsize=16)

plt.savefig(r'F:\\Seer_project\\figures\\scatter_biomarker_NP_ratio_tar.svg', dpi=800,bbox_inches='tight')
plt.show()


# In[ ]:




