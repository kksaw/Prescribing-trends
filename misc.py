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