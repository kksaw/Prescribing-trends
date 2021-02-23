# -*- coding: utf-8 -*-
"""
Created on Sun Jan 3 23:57:05 2021

@author: khaik
"""

import sys
import os
import json
import requests

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import geopandas as gpd
import pandas as pd
import fingertips_py as ftp

import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std

import functions as fn
import pickle


with open('druglist.json', 'r') as f:
    druglist = json.load(f)

#estimate population -- need to make this more flexible so that it works when OP adds new data
#fpopn needs to be array
def getPopn():
    popn = [54786300, 55268100, 55619400, 55977200, 56287000]
    months = np.arange(0,12*len(popn),12)
    f = interp1d(months, popn, fill_value='extrapolate')
    estpopn = f(np.arange((len(popn)+1)*12))
    return estpopn[5:65]/1000
fpopn = getPopn()


def getAPI(url):
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception('oops')
    return response.json()

def jprint(data):
    print(json.dumps(data, sort_keys=True, indent=4))
    
test = getAPI("https://openprescribing.net/api/1.0/spending_by_ccg/?code=0403010B0&format=json")


'''
TIME SERIES by drug name -- ok
1. Number of prescription/1000 popn
2. Cost/1000 popn normalised with inflation
3. (Maybe) quantity...? whatever this means..
'''
by_drugs = []
for BNFid in druglist:    
    url = "https://openprescribing.net/api/1.0/spending/?code=" + BNFid + "&format=json"
    data = getAPI(url)
    
    by_drugs.append({'id':BNFid, 'spending':[sub['actual_cost'] for sub in data],
                     'items':[sub['items'] for sub in data], 'quantity':[sub['quantity'] for sub in data]})

nmonths = len(data)

#check strange results eg. old drug
baddrug = [a['id'] for a in by_drugs if len(a['items'])!=nmonths]
print('strange results from these drugs!'+ str(baddrug))
for b in baddrug:
    by_drugs = [item for item in by_drugs if item['id']!=b]

#ITEM normalise by 1000 popn   and COST normalise by 1000 popn + inflation
yItems, yCost, yQuantity = [], [], []
for a in by_drugs:
    yItems.append(np.array(a['items']))
    yCost.append(np.array(a['spending']))
    yQuantity.append(np.array(a['quantity']))
yItems = np.vstack(yItems)
yCost = np.vstack(yCost)
yQuantity = np.vstack(yQuantity)

np.savetxt("by_year_items.csv", yItems, delimiter=",")
np.savetxt("by_year_cost.csv", yCost, delimiter=",")
np.savetxt("by_year_quantity.csv", yQuantity, delimiter=",")


def plotTimeSeries(y):
    x = np.arange(nmonths)

    col = np.concatenate([fn.linear_gradient('#393b79', '#9c9ede', n=12)['hex'], 
                        fn.linear_gradient('#637939', '#cedb9c', n=4)['hex'],
                        fn.linear_gradient('#8c6d31', '#e7cd94', n=7)['hex'], 
                        fn.linear_gradient('#843c39', '#e7969c', n=10)['hex']])
    labels = ['Jan 16', 'Jul 16', 'Jan 17', 'Jul 17', 'Jan 18', 'Jul 18', 'Jan 19', 'Jul 19', 'Jan 20', 'Jul 20']
    xticks = np.arange(2,nmonths,6)
    legends = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'TCAs', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 
               '_nolegend_', '_nolegend_', 'MAOIs', '_nolegend_', 
               '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'SSRIs', '_nolegend_', 
               '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'Others', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', ]
    
    fig, ax = plt.subplots()
    ax.stackplot(x,y, colors=col, labels=legends)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel('Year')
    ax.set_ylabel('Number of prescriptions')
    ax.set_title('Number of antidepressant prescriptions/1000 popn')
    ax.legend(loc='upper left')
    plt.show()

for y in [yItems, yCost, yQuantity]:
    plotTimeSeries(y)
  

'''
CCG PROFILE
1. Quantity/1000 pts
2. Cost/1000 pts

how to deal with mergers before 2020, ignore and use current ones?
'''
drugcode = "0403"

with open('ccgdict.json', 'r') as f:
    ccgdict = json.load(f)

fp = "CCG_boundaries\Clinical_Commissioning_Groups__April_2020__EN_BUC.shp"
map_df = gpd.read_file(fp)
map_df['ccg20nm'] = map_df['ccg20nm'].str.upper()
map_df['ccgid'] = map_df['ccg20nm'].map(ccgdict)

url = "https://openprescribing.net/api/1.0/org_details/?org_type=ccg&keys=total_list_size&format=json"
data = getAPI(url)
listsize_ccg = {}
for sub in data:
    if sub['row_id'] in listsize_ccg:
        listsize_ccg[sub['row_id']].append(sub['total_list_size'])
    else:
        listsize_ccg[sub['row_id']] = [sub['total_list_size']]

by_ccgs_general = []
for CCGid in list(listsize_ccg.keys()):
    url = "https://openprescribing.net/api/1.0/spending_by_ccg/?code=" + drugcode + "&org=" + CCGid + "&format=json"
    data = getAPI(url)
    ls = np.array(listsize_ccg[CCGid])
    a, b, c = [sub['items'] for sub in data], [sub['actual_cost'] for sub in data], [sub['quantity'] for sub in data]
    by_ccgs_general.append({'id':CCGid, 'name':data[0]['row_name'], 
                            'items':a, 'items_norm':a/ls,
                            'spending':b, 'spending_norm':b/ls,
                            'quantity':c, 'quantity_norm':c/ls})
    
#n months average [items, spending, quantity, list size]
n = 6
by_ccgs_dict = {}
for sub in by_ccgs_general:
    by_ccgs_dict[sub['id']] = [ np.mean(sub['items_norm'].reshape(-1,n),axis=1),
                                np.mean(sub['spending_norm'].reshape(-1,n),axis=1),
                                np.mean(sub['quantity_norm'].reshape(-1,n),axis=1),
                                np.mean(np.array(listsize_ccg[sub['id']]).reshape(-1,n),axis=1)]    

index = -1      #which months/years?
map_df['items'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][0][index])
map_df['spending'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][1][index])
map_df['quantity'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][2][index])
map_df['list_size'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][3][index])

column = 'items'
def plotCcgPerf(column):
    vmin, vmax = min(map_df[column]), max(map_df[column])
    fig, ax = plt.subplots(1, figsize=(10,6))
    map_df.plot(column=column, cmap='Blues', linewidth=0.5, ax=ax, edgecolor='0.1')
    ax.axis('off')
    ax.set_title(str("CCG "+column))
    sm = plt.cm.ScalarMappable(cmap='Blues', norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm._A = []
    cbar = fig.colorbar(sm)
    fig.savefig(str(column+"_map_export.png"),dpi=200)

for column in ['items', 'spending', 'quantity', 'list_size']:
    plotCcgPerf(column)
    
    
'''
detailed by specific drug -- not done
by_ccgs = [
{'name': CCGname, 
'list_size':[array_by_year]
'0403010B0':[[items_by_year], [spending], [quantity]],
'0403010D0':[[items], [spending], [quantity]],
...]
#filtered_list = [elem for elem in dicts if elem.get('z') == 'apple' and elem.get('aa') == 'banana']    
#ccglist = list(set(row['row_id'] for row in data))
'''
by_ccgs = [{'id':ccglist[i]} for i in range(135)]
for BNFid in druglist:
    url = "https://openprescribing.net/api/1.0/spending_by_ccg/?code=" + BNFid + "&format=json"
    data = getAPI(url)
    for row in by_ccgs:
        row[BNFid] = [sub['actual_cost'] for sub in data]
    
    by_ccgs.append({})


'''
CORRELATION ANALYSIS data (by practice) --  ok

FINGERTIPS
ftp.get_all_profiles() -- match ID to profiles (PHOF (19), GP (20), GP_supporting (21))
ftp.get_all_ages -- match ID to ages (65+ = 27)
ftp.get_all_areas -- match ID to areas (GP = 7)
ftp.get_areas_for_area_type(7) -- get all GPid (# = 6709)
ftp.get_metadata_for_profile_as_dataframe(20) -- returns indicator IDs
ftp.retrieve_data.get_data_by_indicator_ids(ind[i],area_type_id=7,profile_id=21) -- returns indicator value
#(indicator_ids, area_type_id=7, parent_area_type_id=15, profile_id=None, include_sortable_time_periods=None, is_test=False)

20 - [295,93553,336,351,848,90646]
21 - [114]

Variables (profile, indicator) -- fingertips
1. QOF (Apr 19 - Mar 20) -- fingertips (20,295) 6698
2. IMD -- fingertips (20, 93553) 6705
3. Practice list size  -- fingertips (21, 114) 6698
4. % patients over 65 yo -- fingertips (20,336) 6705. 2020 data is available
5. % patients with long-term health disease -- fingertips (20, 351) 6701. 2020 data available
6. % with diagnosis of depression (?) -- fingertips (20, 848) 6698
7. % with new diagnosis of depression -- fingertips (20, 90646) 6698

prescription/1000 pts -- OpenPrescribing (2019 #6775) !! Quantity or Items or Cost???

Quintile range, median prescription per 1000 pts, univariate IRR, multivariate IRR, p value

Let's try 1 year first and exclude GPs that joined/closed during the period
#'2015-12-01'

average of 2019 data
openprescribing: 7017 GPs -> list size 6825 excluded 192 (less than 12 values) -> quantity 6775 excluded 50 (less than 12 values)
fingertips: 6698-6705 GPs

use latest data from fingertips - not sure which date range this is from
'''
year = 2019
def getData(year, depvar='quantity', mode=0):
    ystr=str(year)
    datelist = ['01','02','03','04','05','06','07','08','09','10','11','12']
    drugcode = '0403'
    
    if 2015<year<=2019:
        imd_year=2019
    elif 2010<year<=2015:
        imd_year=2015
    elif year<=2010:
        imd_year=2010
    else:
        sys.exit('Invalid year!')
    
    url = "https://openprescribing.net/api/1.0/org_details/?org_type=practice&keys=total_list_size&format=json"
    data = getAPI(url)
    listsize_gp = {}
    for sub in data:
        if sub['date'][:4] == ystr:
            if sub['row_id'] in listsize_gp:
                listsize_gp[sub['row_id']].append(sub['total_list_size'])
            else:
                listsize_gp[sub['row_id']] = [sub['total_list_size']]
    
    exclude = []            
    for key, value in listsize_gp.items():
        if len(value)<12:
            exclude.append(key)
    [listsize_gp.pop(key,None) for key in exclude];
    
    Qdict = dict.fromkeys(listsize_gp,0)
    Cdict = dict.fromkeys(listsize_gp,0)        
    for date in datelist:
        url = "https://openprescribing.net/api/1.0/spending_by_practice/?code=" + drugcode + "&date=2019-" + date + "-01&format=json"
        data = getAPI(url)
        for sub in data:
            if sub['row_id'] in list(Qdict.keys()):
                if Qdict[sub['row_id']] == 0:
                    Qdict[sub['row_id']] = [sub[depvar]]
                    Cdict[sub['row_id']] = sub['ccg']
                else:
                    Qdict[sub['row_id']].append(sub[depvar])
    
    exclude1 = []
    for key, value in Qdict.items():
        if value == 0:
            exclude1.append(key)
            continue
        if len(value)==12:
            Qdict[key] = round(np.mean(value/np.array(listsize_gp[key])),3)
        else:
            exclude1.append(key)
    
    [Qdict.pop(key,None) for key in exclude1];
    
    exc = [exclude, exclude1]
    
    ind_df = pd.DataFrame(list(Qdict.items()), columns=['ID', 'Quantity'])
    ind_df['CCG'] = ind_df['ID'].map(Cdict)
    ind_df = ind_df[['ID', 'CCG', 'Quantity']]
    
    ind = [295,93553,114,336,351,848,90646]
    keys = ['QOF', 'IMD', 'List Size', '% 65+', '% Long-term Health Conditions', '% Depression', '% New Depression']
    tp = ystr+'/'+str(year+1)[2:]
    
    for i in range(len(ind)):
        if ind[i] == 114:
            df = ftp.retrieve_data.get_data_by_indicator_ids(ind[i],area_type_id=7,profile_id=21)
            df2 = df.loc[(df['Area Type'] == 'GP') & (df['Time period'] == tp)][['Area Code','Value']]
    
        else:
            df = ftp.retrieve_data.get_data_by_indicator_ids(ind[i],area_type_id=7,profile_id=20)
            if ind[i] in [336,351]:
                df2 = df.loc[(df['Area Type'] == 'GP') & (df['Time period'] == year)][['Area Code','Value']]
            elif ind[i] == 93553:
                df2 = df.loc[(df['Area Type'] == 'GP') & (df['Time period'] == imd_year)][['Area Code','Value']]
            else:
                df2 = df.loc[(df['Area Type'] == 'GP') & (df['Time period'] == tp)][['Area Code','Value']]
    
        ind_df[keys[i]] = ind_df['ID'].map(df2.set_index('Area Code')['Value'])
    
    ind_df.columns = ['id', 'ccg', 'quant', 'qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
    
    desc_stat = summary(ind_df)
    
    if mode==0:
        save_path = 'data_'+ystr+'/'
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        filename=os.path.join(save_path,'data.pickle')
        with open(filename, 'wb') as f:
            pickle.dump(ind_df, f)
        filename=os.path.join(save_path,'exclude.pickle')
        with open(filename, 'wb') as f:
            pickle.dump(exc, f)
        filename=os.path.join(save_path,'Qdict.pickle')
        with open(filename, 'wb') as f:
            pickle.dump(Qdict, f)
        filename=os.path.join(save_path,'summary.csv')
        desc_stat.to_csv(filename, index=False)    
    
    elif mode==1:
        return ind_df, exclude, Qdict, desc_stat

ind_df = getData(2019, depvar='quantity', mode=1)[0]
    
def summary(data):
    #75-86 have missing values from fingertips for 2019
    new_df = pd.DataFrame(columns = ['Variable', 'Median', 'P025', 'P075', 'Missing'])
    new_df['Variable'] = data.columns[2:]
    
    n=0
    for col in data.columns[2:]:
        x = data[col]
        if x.isnull().values.any():
            print("%s: %d" %(col, x.isnull().values.sum()), end=' ')
            new_df.iloc[n,4] = x.isnull().values.sum()
        new_df.iloc[n,1] = np.nanmedian(x)
        new_df.iloc[n,2:4] = np.nanpercentile(x, [25 ,75])
        n+=1
    return new_df

def misc(data, var):
    qr = np.zeros((5,2))
    data['q']=pd.qcut(data[var], 5, labels=False)
    for i in range(5):
        df = data[var].loc[data['q']==i]
        qr[i] = min(df), max(df)
    return qr

def get_total_quantity(year, drugcode, depvar='quantity'):
    ystr = str(year)
    url = "https://openprescribing.net/api/1.0/spending/?code="+drugcode+"&format=json"
    data = getAPI(url)
    total = []
    for sub in data:
        if sub['date'][:4] == ystr:
            total.append(sub[depvar])
    total = sum(total)
    return total

total = []
for i in range(2016,2021):
    total.append(get_total_quantity(i, '0403'))
    
#might not be necessary since nbrm ignores rows with nan
def check_inddf(ind_df, keys):
    key1 = list(ind_df['ID'])
    key2 = list(df2['Area Code'])
    keyd = list(set(key1).symmetric_difference(set(key2)))
    dfprob = ind_df.loc[ind_df['ID'].isin(keyd)]
    
    for key in keys:
        dfprob = dfprob.append(ind_df.loc[ind_df[key].isna()])    
    key3 = list(set(dfprob['ID']))
    
    ind_df = ind_df.loc[-ind_df['ID'].isin(key3)]
    
    return dfprob, ind_df


'''LINEAR REGRESSION'''
with open('ind_df_ccg.pickle', 'rb') as f:
    data = pickle.load(f)
keys = ['id', 'ccg', 'quant', 'qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
data.columns = keys

def linreg(data, mode=0):
    data = data.dropna()
    X = np.zeros((len(data),6))
    
    if mode==0:
        for i in range(6):
            X[:,i] = pd.qcut(data[keys[3+i]], 5, labels=False)
    elif mode==1:
        for i in range(3):
            X[:,i] = pd.qcut(data[keys[3+i]], 5, labels=False)
        for i in range(3):
            X[:,3+i] = data[keys[6+i]]
    elif mode==2:
        for i in range(5):
            X[:,i] = pd.qcut(data[keys[3+i]], 5, labels=False)
        X[:,-1] = data[keys[-1]]
    elif mode==3:
        for i in range(3):
            X[:,i] = pd.qcut(data[keys[3+i]], 5, labels=False)
        for i in range(2):
            X[:,3+i] = data[keys[6+i]]
        X[:,-1] = data[keys[-1]]
    
    y = np.array(data['quant'])    
    model = sm.OLS(y,X)
    results = model.fit()
    a = results.params
    b = results.summary()
    return a, b

coeff = pd.DataFrame(index=keys[3:9], columns=['all_quint_dep','quint_cont_dep', 'all_quint_newdep', 'quint_cont_newdep'])
summary = []
for i in range(4):
    a,b = linreg(data,mode=i)
    coeff.iloc[:,i]=a
    summary.append(b)
    
coeff.to_csv('linear regression.csv', index=False)    


