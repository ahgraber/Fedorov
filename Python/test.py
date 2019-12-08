#%% import
# import numpy as np
# import random
from Python.FedorovDesignClass import DesignMatrix

#%% create test design
dm = DesignMatrix(n=8)

#%% test add_attribute
dm.add_attribute(name="age", levels=3, dist=[25,50,25])
dm.add_attribute(name="sex", levels=2, dist=[50,50])
dm.add_attribute(name="bmi", levels=3, dist=[33,33,34])

#%%
dm.generate_design()


#%%
dm.dslacks
#%%
dm.islacks
#%%
dm.X

#%% test add_row
dm.add_row([3,3,3]) # note: throws an error because 3 is not in any of the levels


#%%
dm.dslacks
#%%
dm.X

#%% try again with valid vals
dm.add_row([1,1,1]) 


#%%
dm.dslacks
#%%
dm.X

#%% test del_row
dm.del_row(1)
dm.update_dslacks()


#%%
dm.dslacks
#%%
dm.X
#%%
doptimality(dm)




# %%
