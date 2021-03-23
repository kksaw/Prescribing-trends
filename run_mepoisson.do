cd "C:\Users\khaik\Documents\camb\4th year\SSC\prescribing trends\data\stata data"

import delimited using data_2019

mepoisson items qof_5 imd_5 list_size_5 elderly_popn_5 health_cond_5 dep_5 || ccg:

mepoisson items qof_5 imd_5 list_size_5 elderly_popn health_cond dep || ccg:

mepoisson items qof_5 imd_5 list_size_5 elderly_popn_5 health_cond_5 new_dep_5 || ccg:

mepoisson items qof_5 imd_5 list_size_5 elderly_popn health_cond new_dep || ccg:



