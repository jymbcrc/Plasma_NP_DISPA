#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio.SeqUtils.ProtParam import ProteinAnalysis


# In[2]:


#print("%0.2f" % X.instability_index())


# In[3]:


#print("%0.2f" % X.gravy())


# In[4]:


from Bio import SeqIO
record_dict = SeqIO.index("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\DD\\protein.fas", "fasta")

record_dd = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\DD\\protein.fas", "fasta"))
record_np1 = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\NP1\\protein.fas", "fasta"))
record_np2 = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\NP2\\protein.fas", "fasta"))
record_np3 = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\NP3\\protein.fas", "fasta"))
record_np4 = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\NP4\\protein.fas", "fasta"))
record_np5 = list(SeqIO.parse("E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\NP5\\protein.fas", "fasta"))





len(record_np1)

#print(record_dict.get_raw("gi|1348917|gb|G26685|G26685").decode())


# In[ ]:





# In[5]:


record_dd[1].name.split('|')[1]


# In[ ]:





# In[31]:


record=record_dd
record[122].seq
 
#X = ProteinAnalysis(str(record[122].seq))
#X.gravy()


# In[34]:


len(record[121].seq)


# In[35]:


X = ProteinAnalysis(str(record[121].seq))
X.isoelectric_point()


# In[7]:


#X.gravy()


# In[36]:


def get_gravy_values(record = record_dd):
    '''get Gravy values of a given protein sequence'''
    value = []
    for i in range(len(record)):
        try:
            X = ProteinAnalysis(str(record[i].seq))
            value.append(X.gravy())
        except Exception as e:  # Catch all exceptions (better to specify the exact one if known)
                #print(f"An error occurred with item {i}: {e}. But we're continuing.")
                continue
    return(value)
def get_IP_values(record = record_dd):
    '''get isoelectric_point of a given protein sequence'''
    value = []

    for i in range(len(record)):
        try:
            X = ProteinAnalysis(str(record[i].seq))
            value.append(X.isoelectric_point())
        except Exception as e:  # Catch all exceptions (better to specify the exact one if known)
                #print(f"An error occurred with item {i}: {e}. But we're continuing.")
                continue
    return(value)


# In[37]:


get_IP_values(record = record_np2)


# In[44]:


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(9, 6), edgecolor='#04253a')
# Convert the lists to a single DataFrame
data = {
    'DD': get_gravy_values(record = record_dd),
    'NP1': get_gravy_values(record = record_np1),
    'NP2': get_gravy_values(record = record_np2),
    'NP3': get_gravy_values(record = record_np3),
    'NP4': get_gravy_values(record = record_np4),
    'NP5': get_gravy_values(record = record_np5),
}

# Convert the dict to a list of tuples, then to a DataFrame
long_data = [(key, val) for key, values in data.items() for val in values]
df_gravy = pd.DataFrame(long_data, columns=["list", "value"])

# Plotting the strip plot
sns.stripplot(x="list", y="value", data=df_gravy,palette='viridis',jitter=True)
plt.axhline(y=0, color='red', linestyle='--')
plt.xlabel('')
plt.ylabel("Gravy Value",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.show()


# In[51]:


# Plotting the strip plot

plt.figure(figsize=(8, 6), edgecolor='#04253a')
sns.violinplot(x="list", y="value", data=df_gravy,palette='viridis',jitter=True)
plt.axhline(y=0, color='red', linestyle='--')
plt.xlabel('')
plt.ylabel("Gravy Value",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.show()


# In[48]:


plt.rcParams['font.family'] = 'Arial'
plt.figure(figsize=(9, 6), edgecolor='#04253a')
# Convert the lists to a single DataFrame
data = {
    'DD': get_IP_values(record = record_dd),
    'NP1': get_IP_values(record = record_np1),
    'NP2': get_IP_values(record = record_np2),
    'NP3': get_IP_values(record = record_np3),
    'NP4': get_IP_values(record = record_np4),
    'NP5': get_IP_values(record = record_np5),
}

# Convert the dict to a list of tuples, then to a DataFrame
long_data = [(key, val) for key, values in data.items() for val in values]
df_ip = pd.DataFrame(long_data, columns=["list", "value"])

# Plotting the strip plot
sns.stripplot(x="list", y="value", data=df_ip,palette='viridis',jitter=True)
plt.axhline(y=7.4, color='red', linestyle='--')
plt.xlabel('')
plt.ylabel("Isoelectric_point",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.show()


# In[52]:


plt.figure(figsize=(8, 6), edgecolor='#04253a')

# Plotting the strip plot
sns.violinplot(x="list", y="value", data=df_ip,palette='viridis',jitter=True)
plt.axhline(y=7.4, color='red', linestyle='--')
plt.xlabel('')
plt.ylabel("Isoelectric_point",fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)
plt.show()


# In[16]:


[len(get_gravy_values(record = record_dd)),
len([item for item in get_gravy_values(record = record_dd) if item>0])]


# In[17]:


[len(get_gravy_values(record = record_np1)),
len([item for item in get_gravy_values(record = record_np1) if item>0])]


# In[18]:


[len(get_gravy_values(record = record_np2)),
len([item for item in get_gravy_values(record = record_np2) if item>0])]


# In[19]:


[len(get_gravy_values(record = record_np3)),
len([item for item in get_gravy_values(record = record_np3) if item>0])]


# In[20]:


[len(get_gravy_values(record = record_np4)),
len([item for item in get_gravy_values(record = record_np4) if item>0])]


# In[21]:


[len(get_gravy_values(record = record_np5)),
len([item for item in get_gravy_values(record = record_np5) if item>0])]


# In[ ]:





# In[14]:


#print(record_dict["sp|A0A075B6H9|LV469_HUMAN"].format("fasta"))


# In[13]:


from Bio import SeqIO
filename = "E:\\yuming\\2022\\20221110\\New_folder\\fragpipe2\\DD\\protein.fas"
for record in SeqIO.parse(filename, "fasta"):
    print("ID %s" % record.id)
    print("Sequence length %i" % len(record))
    print("Sequence alphabet %s" % record.seq)


# In[ ]:


#[record_dd[i].name.split('|')[1] for i in range(len(record_dd))]


# # plot Venn diagram between scouting DISPA methods and LC results 

# In[32]:


record_dd[4].name


# In[6]:


list_dd = [record_dd[i].name.split('|')[1] for i in range(len(record_dd))]
list1 = [record_np1[i].name.split('|')[1] for i in range(len(record_np1))]
list2 = [record_np2[i].name.split('|')[1] for i in range(len(record_np2))]
list3 = [record_np3[i].name.split('|')[1] for i in range(len(record_np3))]
list4 = [record_np4[i].name.split('|')[1] for i in range(len(record_np4))]
list5 = [record_np5[i].name.split('|')[1] for i in range(len(record_np5))]


# In[7]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import pandas as pd
import venn


# In[8]:


lc_all=list_dd+list1+list2+list3+list4+list5
len(lc_all)


# In[9]:


sample=[0,2,3,4]
def remove_dul(tar_list=sample):
    res=[]
    [res.append(x) for x in tar_list if x not in res]
    return(res)


# In[10]:


len(remove_dul(tar_list=lc_all))


# In[33]:


# import the scouting DISPA results here 

dfdf=pd.read_excel(r'F:\\Seer_project\\NPs_scout.xlsx',sheet_name='Sheet1')
dfdf.set_index('Unnamed: 0',inplace=True)
#dfdf


# In[34]:


dfdf


# In[45]:


# try to find out all uniprot IDs for all potential protein IDs

import re
nillist=[]
pattern = r'sp\|(\w{6})'
for item in dfdf['old_index'].tolist():
    matches = re.findall(pattern, item)
    nillist.append(matches)
flatted=[item for sublist in nillist for item in sublist]
len(flatted)


# In[26]:


#dfdf['pro_names'].tolist()


# In[53]:


# plot venn diagram

from matplotlib_venn import venn2
matplotlib.rcParams['font.family'] = "Arial"
my_dpi=150
plt.figure(figsize=(600/my_dpi, 600/my_dpi), dpi=my_dpi)#控制图尺寸的同时，使图高分辨率（高清）显示

out=venn2(subsets =[set(dfdf['pro_names'].tolist()), set(remove_dul(tar_list=lc_all))], 
      set_labels =('DISPA', 'LC-MS'),set_colors=('#8098C1', '#E49994'),
      alpha=0.9,
      normalize_to=1)
for text in out.set_labels:
    text.set_fontsize(16)
for text in out.subset_labels:
    text.set_fontsize(16)
    text.set_color('black')
fig_path= 'F:\\Seer_project\\figures\\' 
figname='venn_scout_lc'
plt.savefig(fig_path + "%s.svg" % figname,dpi=800, bbox_inches='tight')  
plt.show()


# In[54]:


list33=[item for item in dfdf['pro_names'].tolist() if item not in remove_dul(tar_list=lc_all)]
list44=[item for item in remove_dul(tar_list=lc_all) if item not in dfdf['pro_names'].tolist()]
list55=[item for item in dfdf['pro_names'].tolist() if item in remove_dul(tar_list=lc_all)]


# In[58]:


len(remove_dul(list33))


# In[60]:


# Calculate max length and fill with NaN if necessary
list11=remove_dul(tar_list=lc_all)
list22=dfdf['pro_names'].tolist()
list333=remove_dul(list33)
list444=remove_dul(list44)
list555=remove_dul(list55)


max_len = max(len(list11), len(list22), len(list333), len(list444), len(list555))
list11 = list11 + [None]*(max_len - len(list11))
list22 = list22 + [None]*(max_len - len(list22))
list333 = list333 + [None]*(max_len - len(list333))
list444 = list444 + [None]*(max_len - len(list444))
list555 = list555 + [None]*(max_len - len(list555))
df = pd.DataFrame({
    'LC': list11,
    'DISPA': list22,
    'shared':list555,
    'dispa_only':list333,
    'LC_only':list444
})

df


# In[62]:


df.to_excel(r'F:\\Seer_project\\scout_LC.xlsx')


# In[46]:


#df.to_csv(r'F:\\Seer_project\\scout_LC.txt')

