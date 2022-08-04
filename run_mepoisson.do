cd "C:\Users\khaik\Documents\camb\4th year\SSC\prescribing trends\data\stata data"

import delimited using data_2019

mepoisson quantity i.rul_5 i.qof_5 i.imd_5 i.list_size_5 i.elderly_popn_5 i.health_cond_5 i.dep_5 || ccg:

mepoisson items rul qof_5 imd_5 list_size_5 elderly_popn_5 health_cond_5 dep_5 || ccg:

mepoisson items rul qof_5 imd_5 list_size_5 elderly_popn health_cond dep || ccg:

mepoisson items rul qof_5 imd_5 list_size_5 elderly_popn_5 health_cond_5 new_dep_5 || ccg:

mepoisson items rul qof_5 imd_5 list_size_5 elderly_popn health_cond new_dep || ccg:
