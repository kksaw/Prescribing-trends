# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 19:00:54 2021

@author: khaik
"""

'''
NEGATIVE BINOMIAL REGRESSION -- ok

formula = 'quant ~ qof + imd + list_size + elderly_popn + health_cond + dep + new_dep'
response, predictors = dmatrices(formula, data, return_type='dataframe')
po_results = sm.GLM(response, predictors, family=sm.families.Poisson()).fit()
print(po_results.summary())

def ct_response(row):
    y = row['quant']
    m = row['mu']
    return ((y-m)**2 - y)/m

#ct_data = data.copy()
ct_data = predictors.copy()
ct_data = ct_data.drop('Intercept', 1)
ct_data.insert(loc=0, column='quant', value=response)

ct_data['mu'] = po_results.mu
ct_data['ct_resp'] = ct_data.apply(ct_response, axis=1)
ct_results = smf.ols('ct_resp ~ mu-1', ct_data).fit()
alpha_ci95 = ct_results.conf_int(0.05).loc['mu']
print('\nC-T dispersion test: alpha = {:5.3f}, 95% CI = ({:5.3f}, {:5.3f})'.format(ct_results.params[0], alpha_ci95.loc[0], alpha_ci95.loc[1]))
#why is alpha negative...?

#alpha default 1
nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=0.15)).fit()
print(nb_results.summary())
nb_irr = np.exp(nb_results.params.values[1:])
nb_confinf = np.exp(nb_results.conf_int())
nb_rsquared =  1 - nb_results.llf/nb_results.llnull
'''
import pandas as pd
import numpy as np
from patsy import dmatrices
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats.distributions import chi2
import pickle

with open('ind_df_ccg.pickle', 'rb') as f:
    data = pickle.load(f)

cols = ['ccg', 'quant', 'qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
keys = ['qof', 'imd', 'list_size', 'elderly_popn', 'health_cond', 'dep', 'new_dep']
data = ind_df[cols]
dof = len(cols)-2

def nbr(data, n, a=1, dof=dof, keys=keys):
    var = keys[n]
    new_df = pd.DataFrame(columns=['Variable', 'Vmin', 'Vmax', 'Qmed', 'uIRR', 'uCI025', 'uCI975', 'mIRR', 'mCI025', 'mCI975'])
    data[str(var+'5')] = pd.qcut(data[var], 5, labels=False)
    f1 = 'quant ~ ' + str(var)
    f2 = 'quant ~ ccg + qof + imd + list_size + elderly_popn + health_cond + dep + new_dep'
    
    for i in range(5):
        df = data.loc[data[str(var+'5')]==i]
        
        response, predictors = dmatrices(f1, df, return_type='dataframe')
        nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=a)).fit()
        u_irr = np.exp(nb_results.params.values[1])      #save
        u_confint = np.exp(nb_results.conf_int().values[1])    #save
        #u_rsquared =  1 - nb_results.llf/nb_results.llnull
                
        response, predictors = dmatrices(f2, df, return_type='dataframe')
        nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=a)).fit()
        m_irr = np.exp(nb_results.params.values[n-len(keys)])      #save
        m_confint = np.exp(nb_results.conf_int().values[n-len(keys)])    #save
        #m_rsquared =  1 - nb_results.llf/nb_results.llnull
        
        if i == 0:
            ref = {'uIRR':u_irr, 
                   'uCI025':u_confint[0], 
                   'uCI975':u_confint[1], 
                   'mIRR':m_irr, 
                   'mCI025':m_confint[0], 
                   'mCI975':m_confint[1]}
            
        results = pd.DataFrame({'Variable':var+str(i), 
                 'Vmin':np.min(df[var]), 
                 'Vmax':np.max(df[var]), 
                 'Qmed':np.median(df['quant']), 
                 'uIRR':u_irr/ref['uIRR'], 
                 'uCI025':u_confint[0]/ref['uCI025'], 
                 'uCI975':u_confint[1]/ref['uCI975'], 
                 'mIRR':m_irr/ref['mIRR'], 
                 'mCI025':m_confint[0]/ref['mCI025'], 
                 'mCI975':m_confint[1]/ref['mCI975']}, index=[0])
        
        new_df = pd.concat([new_df, results])
            
    response, predictors = dmatrices(f1, data, return_type='dataframe')
    nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=a)).fit()
    ull = nb_results.llf
    response, predictors = dmatrices(f2, data, return_type='dataframe')
    nb_results = sm.GLM(response, predictors, family=sm.families.NegativeBinomial(alpha=a)).fit()
    mll = nb_results.llf
    LR = 2*(mll-ull)
    p = chi2.sf(LR, dof)            #save; pvalue is without quintile??

    return new_df, p

results_df = pd.DataFrame(columns=['Variable', 'Vmin', 'Vmax', 'Qmed', 'uIRR', 'uCI025', 'uCI975', 'mIRR', 'mCI025', 'mCI975'])
results_p = []
for n in range(len(keys)):
    new_df, pvalue = nbr(data, n)
    results_df = pd.concat([results_df, new_df])
    results_p.append(pvalue)

results_df.to_csv('nbr_results.csv', index=False)    

'''  NHS BSA
import urllib.request
a = 'https://opendata.nhsbsa.net/api/3/action/datastore_search?resource_id=EPD_202010&limit=5&q=title:jones'  
with urllib.request.urlopen(a) as url:
    fileobj = url.read()
    print(fileobj)

fileobj = urllib.urlopen(url)
print(fileobj.read())
'''




'''inflation'''
index = [101.9, 101.4, 102.1, 102.5, 102.9, 103.3, 103.3, 103.2, 103.8, 104.1, 104.2, 104.6,
104.9, 104.4, 104.9, 105, 105.4, 105.8, 105.8, 105.8, 106.5, 106.6, 106.7, 107, 107.1,
106.3, 106.8, 107, 107.6, 107.9, 107.9, 107.9, 108.4, 108.5, 108.3, 108.5, 108.5, 108.2,
108.6, 108.6, 108.5, 108.5, 108.6, 109.1, 108.6, 109.1, 109.1, 108.9, 109.2, 109, 109.1, 
109.4, 110.1, 110.8, 111.3, 111.3, 112.1, 112.4, 113.6, 114.5]
with open('CPI_index.pickle', 'wb') as f:
    pickle.dump(index, f)

index_old = [99.9, 100.1, 100.4, 100.6, 100.8, 101, 100.9, 101.2, 101.5, 101.6, 
101.8, 102.2, 101.8, 102.4, 102.7, 103.2, 103.5, 103.5, 103.5, 104, 104.3, 104.4,
104.7, 105, 104.5, 104.9, 105.1, 105.5, 105.9, 105.9, 105.9, 106.5, 106.6, 106.7, 
106.9, 107.1, 106.4, 106.8, 107, 107.6, 107.9, 107.9, 108, 108.3, 108.4, 108.3, 
108.5, 108.5, 108.3, 108.6, 108.6, 108.6, 108.6, 108.8, 109.2, 108.8, 109.2, 109.2,
109.1, 109.4]
