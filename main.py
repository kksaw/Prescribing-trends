# -*- coding: utf-8 -*-
"""
Created on Sun Jan 3 23:57:05 2021

@author: khaik
"""

import sys, os
import json, requests, pickle
from datetime import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from scipy.interpolate import interp1d
import geopandas as gpd
import pandas as pd
import fingertips_py as ftp

from patsy import dmatrices
import statsmodels.api as sm

os.chdir("C:\\Users\\khaik\\Documents\\camb\\4th year\\SSC\\prescribing trends\\data")
import functions as fn

#estimate population -- need to make this more flexible so that it works when OP adds new data
#fpopn needs to be array
def getPopn(n, save=False):
    popn = [54786300, 55268100, 55619400, 55977200, 56287000]
    months = np.arange(0,12*len(popn),12)
    f = interp1d(months, popn, fill_value='extrapolate', kind='quadratic')
    estpopn = f(np.arange((len(popn)+1)*12))
    fpopn = estpopn[n:n+60]

    labels = ['Jul 15', 'Jan 16', 'Jul 16', 'Jan 17', 'Jul 17', 'Jan 18', 'Jul 18', 'Jan 19', 'Jul 19', 'Jan 20', 'Jul 20', 'Jan 21']

    fig, ax = plt.subplots()
    ax.plot(estpopn, label='Extrapolated population')
    ax.plot(np.arange(n,n+60),fpopn, label='Time period used in study')
    ax.plot(months, popn,'o', label='Actual population estimates')
    ax.legend(loc='lower right')
    ax.set_xlabel('Time Period')
    ax.set_xticks(np.arange(0,len(estpopn),6))
    ax.set_xticklabels(labels, rotation=45)
    ax.set_ylabel('Population')
    ax.set_title('Extrapolated population from England Mid-Year Population Estimates')
    plt.show()
    
    if save:
        fn.savecsv('1 time series/', ['popn.csv', 'popn_full.csv'], [fpopn, estpopn])
            
    return fpopn/1000


'''
TIME SERIES by drug name -- ok
1. Number of prescription/1000 popn
2. Cost/1000 popn normalised with inflation
3. (Maybe) quantity...? whatever this means..
'''
def plotTimeSeries(y, nmonths, xstart=0, ylabel='', title='', exclude_amtt=False):
    x = np.arange(nmonths)
    ntca = 11 if exclude_amtt else 12        
    col = np.concatenate([fn.linear_gradient('#393b79', '#9c9ede', n=ntca)['hex'], 
                        fn.linear_gradient('#637939', '#cedb9c', n=4)['hex'],
                        fn.linear_gradient('#8c6d31', '#e7cd94', n=7)['hex'], 
                        fn.linear_gradient('#843c39', '#e7969c', n=10)['hex']])
    labels = ['Jan 16', 'Jul 16', 'Jan 17', 'Jul 17', 'Jan 18', 'Jul 18', 'Jan 19', 'Jul 19', 'Jan 20', 'Jul 20']
    xticks = np.arange(xstart,nmonths,6)
    if exclude_amtt:
        legends = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'TCAs', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 
                   '_nolegend_', '_nolegend_', 'MAOIs', '_nolegend_', 
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'SSRIs', '_nolegend_', '_nolegend_',
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'Others', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', ]    
    else:
            legends = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'TCAs', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 
                   '_nolegend_', '_nolegend_', 'MAOIs', '_nolegend_', 
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'SSRIs', '_nolegend_', '_nolegend_',
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'Others', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', ]

    fig, ax = plt.subplots()
    ax.stackplot(x,y, colors=col, labels=legends)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel('Date')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='upper left')
    fig.savefig(str("1 time series/"+ylabel+".png"))

def plotTimeSeries_line(y, nmonths, xstart=0, ylabel='', title='', exclude_amtt=False):
    x = np.arange(nmonths)
    ntca = 11 if exclude_amtt else 12        
    col = np.concatenate([fn.linear_gradient('#393b79', '#9c9ede', n=ntca)['hex'], 
                        fn.linear_gradient('#637939', '#cedb9c', n=4)['hex'],
                        fn.linear_gradient('#8c6d31', '#e7cd94', n=7)['hex'], 
                        fn.linear_gradient('#843c39', '#e7969c', n=10)['hex']])
    labels = ['Jan 16', 'Jul 16', 'Jan 17', 'Jul 17', 'Jan 18', 'Jul 18', 'Jan 19', 'Jul 19', 'Jan 20', 'Jul 20']
    xticks = np.arange(xstart,nmonths,6)
    if exclude_amtt:
        legends = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'TCAs', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 
                   '_nolegend_', '_nolegend_', 'MAOIs', '_nolegend_', 
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'SSRIs', '_nolegend_', '_nolegend_',
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'Others', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', ]    
    else:
            legends = ['_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'TCAs', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 
                   '_nolegend_', '_nolegend_', 'MAOIs', '_nolegend_', 
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'SSRIs', '_nolegend_', '_nolegend_',
                   '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', 'Others', '_nolegend_', '_nolegend_', '_nolegend_', '_nolegend_', ]

    fig, ax = plt.subplots()
    for i in range(len(y)):
        ax.plot(x,y[i], color=col[i], label=legends[i])
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_xlabel('Date')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='upper left')
    fig.savefig(str("1 time series/"+ylabel+"_line.png"))

    
def time_series_by_drug(n=6, save=False, exclude_amtt=False):
    fpopn = getPopn(n, save=save)
    with open('druglist.json', 'r') as f:
        druglist = json.load(f)
    
    by_drugs = []
    for BNFid in druglist:    
        url = "https://openprescribing.net/api/1.0/spending/?code=" + BNFid + "&format=json"
        data = fn.getAPI(url)
        
        by_drugs.append({'id':BNFid, 'spending':[sub['actual_cost'] for sub in data],
                         'items':[sub['items'] for sub in data], 'quantity':[sub['quantity'] for sub in data]})
    
    #check strange results eg. old drug
    nmonths = len(data)
    baddrug = [a['id'] for a in by_drugs if len(a['items'])!=nmonths]
    print('strange results from these drugs!'+ str(baddrug))
    
    if exclude_amtt:
        baddrug.append('0403010B0')
        
    for b in baddrug:
        by_drugs = [item for item in by_drugs if item['id']!=b]
    
    #ITEM normalise by 1000 popn   and COST normalise by 1000 popn + inflation
    with open('CPI_index.pickle', 'rb') as f:
        cpi_index = pickle.load(f)
    cpi_index = np.array([x/cpi_index[0] for x in cpi_index])
            
    yItems, yCost, yQuantity = [], [], []
    for a in by_drugs:
        yItems.append(np.array(a['items']))
        yCost.append(np.array(a['spending']))
        yQuantity.append(np.array(a['quantity']))
        
    yItems = np.vstack(yItems)
    yCost = np.vstack(yCost)
    yQuantity = np.vstack(yQuantity)
    nItems = yItems/fpopn
    nCost = (yCost/fpopn)/cpi_index
    nQuantity = yQuantity/fpopn
    
    var = [nItems, nCost, nQuantity]
    ylabels = ['Number of prescriptions issued','Cost','Amount of drug prescribed']
    titles = [str(y+'/1000 popn') for y in ylabels]
    
    [plotTimeSeries(var[i], nmonths, ylabel=ylabels[i], title=titles[i], exclude_amtt=exclude_amtt) for i in range(3)]
    
    if save:
        filenames = ['cpi_index.csv', 'items_original.csv', 'cost_original.csv', 'quantity_original.csv',
                     'items_norm.csv', 'cost_norm.csv', 'quantity_norm.csv']
        files = [cpi_index, yItems, yCost, yQuantity,
                 nItems, nCost, nQuantity]
        fn.savecsv('1 time series/', filenames, files, ftype='array')

time_series_by_drug(n=6, save=True, exclude_amtt=True)  


'''
CCG PROFILE
1. Quantity/1000 pts
2. Cost/1000 pts

how to deal with mergers before 2020, ignore and use current ones?
'''

def plotCcgPerf(column, map_df, timeperiod):
    fig, ax = plt.subplots(1, figsize=(10,10))
    map_df.plot(column=column, cmap='Blues', scheme='QUANTILES', k=5, linewidth=0.5, ax=ax, edgecolor='0.1', legend=True)
    ax.set_title(str("CCG " + column + ' ' + timeperiod[-1]))
    ax.axis('off')

    save_path = '2 geographical variation/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    filename=os.path.join(save_path,str(column+'_map_export.png'))
    fig.savefig(filename,dpi=200)

def ccg_profile():
    drugcode = "0403"
    
    with open('ccgdict.json', 'r') as f:
        ccgdict = json.load(f)
    
    fp = "CCG_boundaries\Clinical_Commissioning_Groups__April_2020__EN_BUC.shp"
    map_df = gpd.read_file(fp)
    map_df['ccg20nm'] = map_df['ccg20nm'].str.upper()
    map_df['ccgid'] = map_df['ccg20nm'].map(ccgdict)
    
    url = "https://openprescribing.net/api/1.0/org_details/?org_type=ccg&keys=total_list_size&format=json"
    data = fn.getAPI(url)
    listsize_ccg = {}
    for sub in data:
        if sub['row_id'] in listsize_ccg:
            listsize_ccg[sub['row_id']].append(sub['total_list_size'])
        else:
            listsize_ccg[sub['row_id']] = [sub['total_list_size']]
    
    by_ccgs_general = []
    for CCGid in list(listsize_ccg.keys()):
        url = "https://openprescribing.net/api/1.0/spending_by_ccg/?code=" + drugcode + "&org=" + CCGid + "&format=json"
        data = fn.getAPI(url)
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
    
    timeperiod = ['2016 Q1', '2016 Q2', '2016 Q3', '2016 Q4',
                  '2017 Q1', '2017 Q2', '2017 Q3', '2017 Q4',
                  '2018 Q1', '2018 Q2', '2018 Q3', '2018 Q4',
                  '2019 Q1', '2019 Q2', '2019 Q3', '2019 Q4',
                  '2020 Q1', '2020 Q2', '2020 Q3', '2020 Q4']

    for column in ['items', 'spending', 'quantity', 'list_size']:
        plotCcgPerf(column, map_df, timeperiod)
        
ccg_profile()

'''by drugtype'''
def map_date(datestr, firstdate=datetime(2016,1,1)):
    d1 = firstdate
    d2 = datetime.strptime(datestr,'%Y-%m-%d')
    return (d2.year - d1.year) * 12 + d2.month - d1.month

def plotCcgPerf_anim(by_ccgs_dict, map_df, column, ind, timeperiod=timeperiod):
    fig, ax = plt.subplots(1, figsize=(10,10))
    def plotMap(i):
        ax.clear()
        ax.axis('off')
        map_df[column] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][ind][i])
        ax.set_title(str("CCG " + column + ' ' + timeperiod[i]))
        map_df.plot(column=column, cmap='Blues', scheme='QUANTILES', k=5, linewidth=0.5, ax=ax, edgecolor='0.1', legend=True)
        
    animator = ani.FuncAnimation(fig, plotMap, interval=20)
    writergif = ani.PillowWriter(fps=5) 
    save_path = '2 geographical variation/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    filename=os.path.join(save_path,str(column+'.gif'))
    animator.save(filename, writer=writergif)
    
def plotCcgHist(by_ccgs_dict, cols, ind, timeperiod=timeperiod):
    a = []
    fig, axs = plt.subplots(2,2, figsize=(12,10))
    axs = axs.ravel()
    for i in range(4):
        a.append(axs[i].hist([row[i][-1] for row in list(by_ccgs_dict.values())], bins=15))
        axs[i].set_xlabel(cols[i])
        axs[i].set_ylabel("Number of CCGs")
    fig.suptitle(str("CCG histogram " + timeperiod[ind]))
    fig.subplots_adjust(top=0.93, bottom=0.085, hspace=0.3, wspace=0.2)
        
    save_path = '2 geographical variation/'
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    filename=os.path.join(save_path,'CCG histogram.png')
    fig.savefig(filename)
    
    return a


def ccg_profile_1(n=3, animation=False, save=False, excl_amtt=True, hist_return=False):
    with open('druglist.json', 'r') as f:
        druglist = json.load(f)
    if excl_amtt:
        [druglist.remove(x) for x in ['0403030Y0', '0403010B0']];
    else:
        [druglist.remove(x) for x in ['0403030Y0']];
    
    with open('ccgdict.json', 'r') as f:
        ccgdict = json.load(f)
        
    with open('CPI_index.pickle', 'rb') as f:
        cpi_index = pickle.load(f)
    cpi_index = np.array([x/cpi_index[0] for x in cpi_index])

    
    fp = "CCG_boundaries\Clinical_Commissioning_Groups__April_2020__EN_BUC.shp"
    map_df = gpd.read_file(fp)
    map_df['ccg20nm'] = map_df['ccg20nm'].str.upper()
    map_df['ccgid'] = map_df['ccg20nm'].map(ccgdict)
    
    url = "https://openprescribing.net/api/1.0/org_details/?org_type=ccg&keys=total_list_size&format=json"
    data = fn.getAPI(url)
    listsize_ccg = {}
    for sub in data:
        if sub['row_id'] in listsize_ccg:
            listsize_ccg[sub['row_id']].append(sub['total_list_size'])
        else:
            listsize_ccg[sub['row_id']] = [sub['total_list_size']]
    
    #4320 elements (135ccgs * 32 drugs)
    by_ccgs1 = []
    
    by_ccgs_general = []
    total = np.zeros((3,60))
    for CCGid in list(listsize_ccg.keys()):
        ls = np.array(listsize_ccg[CCGid])
        for drugcode in druglist:
            url = "https://openprescribing.net/api/1.0/spending_by_ccg/?code=" + drugcode + "&org=" + CCGid + "&format=json"
            data = fn.getAPI(url)
            
            temp = np.zeros((3,60))
            ind = list(map(map_date, np.array([sub['date'] for sub in data])))
            temp[:,ind]=np.array([sub['items'] for sub in data]), np.array([sub['actual_cost'] for sub in data]), np.array([sub['quantity'] for sub in data])
            by_ccgs1.append(temp)
            total+=temp
            
        by_ccgs_general.append({'id':CCGid, 'name':data[0]['row_name'], 
                         'items':total[0], 'items_norm':total[0]/ls,
                         'spending':total[1], 'spending_norm':(total[1]/ls)/cpi_index,
                         'quantity':total[2], 'quantity_norm':total[2]/ls})
    
    if save:
        save_path = '2 geographical variation/'
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        filename=os.path.join(save_path,'by_ccgs_by_drugs.pickle')
        with open(filename, 'wb') as f:
            pickle.dump(by_ccgs1, f)
        filename=os.path.join(save_path,'by_ccgs_general.pickle')
        with open(filename, 'wb') as f:
            pickle.dump(by_ccgs_general, f)

    #n months average [items, spending, quantity, list size]
    by_ccgs_dict = {}
    for sub in by_ccgs_general:
        by_ccgs_dict[sub['id']] = [ np.mean(sub['items_norm'].reshape(-1,n),axis=1),
                                    np.mean(sub['spending_norm'].reshape(-1,n),axis=1),
                                    np.mean(sub['quantity_norm'].reshape(-1,n),axis=1),
                                    np.mean(np.array(listsize_ccg[sub['id']]).reshape(-1,n),axis=1)]        
    
    timeperiod = ['2016 Q1', '2016 Q2', '2016 Q3', '2016 Q4',
                  '2017 Q1', '2017 Q2', '2017 Q3', '2017 Q4',
                  '2018 Q1', '2018 Q2', '2018 Q3', '2018 Q4',
                  '2019 Q1', '2019 Q2', '2019 Q3', '2019 Q4',
                  '2020 Q1', '2020 Q2', '2020 Q3', '2020 Q4']

    cols = ['items', 'spending', 'quantity', 'list_size']
    
    #histogram - what to do with outliers?
    hist_summary = plotCcgHist(by_ccgs_dict, cols, -1)
    
    #plot map
    if animation:    
        for i in range(4):
            plotCcgPerf_anim(by_ccgs_dict, map_df, cols[i], i, timeperiod)
        
    else:
        index = -1      #which months/years?
        map_df['items'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][0][index])
        map_df['spending'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][1][index])
        map_df['quantity'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][2][index])
        map_df['list_size'] = map_df['ccgid'].map(lambda x: by_ccgs_dict[x][3][index])
        
        for column in cols:
            plotCcgPerf(column, map_df, timeperiod)
    
    if hist_return:
        return hist_summary


def check_SSRI(by_ccgs1, listsize_ccg, Cdict, Qdict, drugcode=21, ndrugs=32):        
    indices = np.array(np.arange(drugcode,4320,ndrugs))
    sel_drug = np.asarray([by_ccgs1[i][1] for i in indices])
    
    nmonths = np.shape(sel_drug)[1]
    p = np.argmax(sel_drug)
    r = int(p/nmonths)
    c = int(p - (r*nmonths))
    max_ccg = list(listsize_ccg.keys())[r]
    maxv = sel_drug[(r,c)]
    print("%s: %d, month: %d" %(max_ccg, maxv, c), end=' ')
    
    gps = [k for k,v in Cdict.items() if v == max_ccg]
    gpv = [Qdict[a][1][drugcode] for a in gps]
    Qdict[gps]

    return 
        
    
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
def summary(data, depvar, dropna=False):
    new_df = pd.DataFrame(columns = ['Variable', 'Median', 'P025', 'P075', 'Missing'])
    new_df['Variable'] = data.columns[1:]
    
    data = data[data.iloc[:,4].notna()] if dropna else data
    
    n=0
    for col in data.columns[1:]:
        x = data[col]
        if x.isnull().values.any():
            print("%s: %d" %(col, x.isnull().values.sum()), end=' ')
            new_df.iloc[n,4] = x.isnull().values.sum()
        new_df.iloc[n,1] = np.nanmedian(x)
        new_df.iloc[n,2:4] = np.nanpercentile(x, [25 ,75])
        n+=1
    return new_df

def quintile_summary(data, depvar):
    varcol = ['Variable', 'Qmin', 'Qmax']+[str('Med_'+var) for var in depvar]
    qs_df = pd.DataFrame(columns = varcol)
    qs_df['Variable'] = [col+str(i+1) for col in data.columns[4:] for i in range(5)]

    n=0
    for col in data.columns[4:11]:
        data['q']=pd.qcut(data[col], 5, labels=False)
        for i in range(5):
            df = data[depvar+[col]].loc[data['q']==i]
            qs_df.iloc[n+i,1:3] = min(df[col]), max(df[col])
            qs_df.iloc[n+i,3:] = np.nanmedian(df[depvar], axis=0)
        n+=5
        
    return qs_df

col='list_size'
data['q']=pd.qcut(data[col], 5, labels=False)
df = data[depvar+[col]].loc[data['q']==0]

def check_Qdict(Qdict):
    a = []
    for key, value in Qdict.items():
        #mean of all drugs in each month (3x12 values), flag if long period of 0? not necessary since already checked list size
        check1 = np.nanmean(value[0], axis=0)
        a.append(np.count_nonzero(check1==0))
        
    a = np.array(a)
    ind = np.where(a!=0)[0]
    nzero = a[ind]
    print("Potential GP closures: %d" %(len(ind)), end='\n')
    return np.array((ind, nzero))

# need to adjust for inflation?? since it's only within 1 year. if adjust_inflation==True, adjust with ref to 2016
year=2019
def getData(year, depvar=['items', 'actual_cost', 'quantity'], 
            save=True, excl_amtt=True, adjust_inflation=False,
            dropna=False):
    
    ystr=str(year)
    datelist = ['01','02','03','04','05','06','07','08','09','10','11','12']
    
    with open('druglist.json', 'r') as f:
        druglist = json.load(f)
    if excl_amtt:
        [druglist.remove(x) for x in ['0403030Y0', '0403010B0']];
    else:
        [druglist.remove(x) for x in ['0403030Y0']];
        
    with open('CPI_index.pickle', 'rb') as f:
        cpi_index = pickle.load(f)
    cpi_index = np.array([x/cpi_index[0] for x in cpi_index])
    ind = np.arange(0,12)+((year-2016)*12)
    cpi_index = cpi_index[ind]
    
    if 2015<year<=2020:
        imd_year=2019
    elif 2010<year<=2015:
        imd_year=2015
    elif year<=2010:
        imd_year=2010
    else:
        sys.exit('Invalid year!')
    
    url = "https://openprescribing.net/api/1.0/org_details/?org_type=practice&keys=total_list_size&format=json"
    data = fn.getAPI(url)
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
    
    nvar, ndate, ndrug, kgp = len(depvar), len(datelist), len(druglist), list(listsize_gp.keys())
    Qdict = {k : np.zeros((nvar, ndrug, ndate)) for k in kgp}
    Cdict = {k : 0 for k in kgp}
    drugcount=1
    for i in range(ndrug):
        drugcode = druglist[i]
        for j in range(ndate):
            date = datelist[j]
            url = "https://openprescribing.net/api/1.0/spending_by_practice/?code=" + drugcode + "&date=" + ystr + "-" + date + "-01&format=json"
            data = fn.getAPI(url)
            for sub in data:
                key = sub['row_id']
                if key in kgp:
                    Qdict[key][:,i,j] = np.array([sub.get(k) for k in depvar])
                    
                    if Cdict[key]==0:
                        Cdict[key] = sub['ccg']
        drugcount+=1
        print(drugcount, end=' ')
        
    print('Qdict done', end=' ')
    
    Mdict = {k : np.zeros(3) for k in kgp}
    for key, value in Qdict.items():
        norm_value = (value/listsize_gp[key])*1000
        norm_value[1] = norm_value[1]/cpi_index if adjust_inflation else norm_value[1]            
        Mdict[key] = np.nanmean(norm_value, axis=(1,2))
        
    z = check_Qdict(Qdict)      #remove??
    
    ind_df = pd.DataFrame.from_dict(Mdict, orient='index', columns=depvar)
    ind_df['ccg'] = ind_df.index.to_series().map(Cdict)
    ind_df = ind_df[['ccg']+depvar]
    
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
    
        ind_df[keys[i]] = ind_df.index.to_series().map(df2.set_index('Area Code')['Value'])
    
    ind_df.columns = ['ccg']+depvar+['qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
    
    print('ind_df done', end='\n')
    
    desc_stat = summary(ind_df, depvar, dropna=dropna)
    desc_quint = quintile_summary(ind_df, depvar)   #doesnt assign nan to quintiles so don't need to dropna
    
    if save:
        save_path = 'data_'+ystr+'/'
        
        filenames = ['data.pickle', 'exclude.pickle', 'Qdict.pickle', 'Cdict.pickle', 'listsize.pickle']
        files = [ind_df, exclude, Qdict, Cdict, listsize_gp]    
        fn.savepickle(save_path, filenames, files)

        filenames = ['summary.csv', 'quintile_summary.csv']
        files = [desc_stat, desc_quint]
        fn.savecsv(save_path, filenames, files, ftype='dataframe')
   
    else:
        return {'data': ind_df, 'exclude': exclude, 'Qdict': Qdict, 'Cdict': Cdict, 
                'listsize': listsize_gp, 'summary': desc_stat, 'quintile_summary': desc_quint}

getData(2019, depvar=['items', 'actual_cost', 'quantity'], save=True, excl_amtt=True)
results_dict = getData(2019, depvar=['items', 'actual_cost', 'quantity'], save=False, excl_amtt=True)
    
[getData(i) for i in range(2017,2020)]

def get_total_quantity(year, depvar=['item', 'actual_cost', 'quantity'], excl_amtt=True):
    with open('druglist.json', 'r') as f:
        druglist = json.load(f)
    if excl_amtt:
        [druglist.remove(x) for x in ['0403030Y0', '0403010B0']];
    else:
        [druglist.remove(x) for x in ['0403030Y0']];

    ystr = str(year)

    total = np.zeros(len(depvar))
    for drugcode in druglist:
        url = "https://openprescribing.net/api/1.0/spending/?code="+drugcode+"&format=json"
        data = fn.getAPI(url)
        
        for sub in data:
            if sub['date'][:4] == ystr:
                total += np.array([sub.get(k) for k in depvar])

    return total

total = []
for i in range(2016,2021):
    total.append(get_total_quantity(i))

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


'''REGRESSION MODEL
Modes
0: all quintile, % depression
1: qof, imd, list_size quintile; the rest continuous, % depression
2: all quintile, % new depression
3: qof, imd, list_size quintile; the rest continuous, % new depression
'''
with open('ind_df_ccg.pickle', 'rb') as f:
    data = pickle.load(f)
keys = ['id', 'ccg', 'quant', 'qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
data.columns = keys

def regmodel(data, mode=0):
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
    return X, y

X, y = regmodel(data)

model = sm.OLS(y,X)
results = model.fit()
a = results.params
b = results.summary()

coeff = pd.DataFrame(index=keys[3:9], columns=['all_quint_dep','quint_cont_dep', 'all_quint_newdep', 'quint_cont_newdep'])
summary = []
for i in range(4):
    a,b = linreg(data,mode=i)
    coeff.iloc[:,i]=a
    summary.append(b)
    
coeff.to_csv('linear regression.csv', index=False) 


def preprocess_data_for_stata(year):
    ystr = str(year)
    keys = ['qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
    
    input_file = 'data_' + ystr + '/data.pickle'
    with open(input_file, 'rb') as f:
        data = pickle.load(f)
    for i in range(len(keys)):
        data[str(keys[i]+'_5')] = pd.qcut(data[keys[i]], 5, labels=False)
    
    output_file = 'data_' + ystr + '.csv'
    data.to_csv(os.path.join('stata data', output_file), index=True, index_label='gp')


