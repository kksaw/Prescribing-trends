# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 22:00:28 2021

@author: khaik
"""
'''
DRUG DICTIONARY - 34

TCA: '0403010x0' where x = B, C, F, J, L, N, R, S, T, V, X, Y
MAOIS: '0403020x0' where x = H, K, M, Q
SSRIs: '0403030x0' where x = D, Z, Y, X, E, L, P, Q
Others: '0403040x0' where x = Z, Y, F, X, T, R, U, S, W, AB

'''

import json
from main import getAPI

druglist = []
a = ['B', 'C', 'F', 'J', 'L', 'N', 'R', 'S', 'T', 'V', 'X', 'Y']
b = ['H', 'K', 'M', 'Q']
c = ['D', 'Z', 'Y', 'X', 'E', 'L', 'P', 'Q']
d = ['Z', 'Y', 'F', 'X', 'T', 'R', 'U', 'S', 'W', 'AB']
for letters in a:
    druglist.append('0403010'+letters+'0')
for letters in b:
    druglist.append('0403020'+letters+'0')
for letters in c:
    druglist.append('0403030'+letters+'0')
for letters in d[:-1]:
    druglist.append('0403040'+letters+'0')
druglist.append('0403040'+d[-1])
    
with open('druglist.json', 'w') as f:
    json.dump(druglist,f)


url = "https://openprescribing.net/api/1.0/org_details/?org_type=ccg&keys=total_list_size&format=json"
data = getAPI(url)
ccgdict = {}
for sub in data:
    if sub['row_name'] in ccgdict:
        continue
    else:
        ccgdict[sub['row_name']] = sub['row_id']

with open('ccgdict.json', 'w') as f:
    json.dump(ccgdict,f)

