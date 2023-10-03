#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


data = pd.read_excel('us_treasury_data.xlsx')
data


# In[3]:


data['Maturity'] = pd.to_datetime(data['Maturity'])
data['Issue Date'] = pd.to_datetime(data['Issue Date'])


# In[4]:


data = data[data['Maturity']>pd.to_datetime('2021-11-02')]
data.to_excel('HW1.xlsx')


# In[ ]:





# In[ ]:





# In[ ]:




