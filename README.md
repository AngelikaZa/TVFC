# TVFC

### Included files:

##### NetworkControl.py: 
Includes functions used to calculate controllability Gramian, and extract dwell time and transition times from TVFC data. 
Please note that the function to calculate minimal control energy was not used in this paper (this was extracted using the matlab functions derived from https://github.com/jastiso/NetworkControl).

##### connectomeFunctions.py:
Includes several functions used to derive connectome-level metrics.

##### demographics.py:
Includes functions to automatically derive and compare group-level demographics and clinical measures

##### optim_fun.m
This function estimates optimal energy required to drive the system from a given initial to a given final state. 
Adjusted from from https://github.com/jastiso/NetworkControl.

##### EnergyCal_Function.m
Calculates the optimal energy of a specific transition using optim_fun and exports the relevant results.

##### CalculateControlEnergy.m
Example script to calculate control energy of a specific transition over the whole cohort.
