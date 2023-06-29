# + Load simulation data
# + Load measurements
# + Compute scaling factor and chi-squared for each parameter set
# + Save/record the scaling factor and chi-squared to Pandas DataFrame (pickle)

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import offset_simulation as off

## Read in the template data
pklName = 'templates/resets_BIG_templates.pkl'
sim = pd.read_pickle(pklName)
print(f"mux:   {sim['mux'].unique()}")
print(f"sigx:  {sim['sigx'].unique()}")
print(f"theta: {sim['theta'].unique()}")
print(f"phi:   {sim['phi'].unique()}")
print(f"diff:  {sim['diff'].unique()}")
print()

# read in the pressure scan data
pp, rst, rstRaw = off.read_pressure_data("data/pressure_scan/20230519/pressure_scan_data.root", fmt='numpy')
pStrs = [f'{p:.0f}' for p in pp]
print(pStrs)

fout = pklName.replace('.pkl', '_chisq.pkl')
print(f'will save chisq values in: {fout}')


if __name__ == '__main__':

    chans = np.arange(16) # zero indexed...?

    # skip channel 1
    ids = np.arange(1, 15, dtype=int)  
    
    nd = len(pp)  # number of datasets (pressures)
    ns = len(sim)  # number of simulations

    print(f'nd, ns = {nd}, {ns}')

    # FIXME: is loop faster if we swap (outer loop is over sim parameters?)
    #        or if you avoid looping altogether and use array-based calculation?
    AA = np.zeros((ns, nd), dtype=float)
    chisq = np.zeros((ns, nd), dtype=float)
    dd = {} # will save output calculations
    for ip in range(nd):  # loop over pressures (data)
        pres = pp[ip]
        print('pres = '+pStrs[ip])
        for ii in range(ns): # loop over simulation parameters
            AA[ii,ip] = off.optimal_scaling_factor(sim['rst'][ii], rst[ip], errors=None, ids=ids)
            chisq[ii,ip] = np.sum( (rst[ip] - sim['rst'][ii]*AA[ii,ip])**2 )

        dd[pStrs[ip]+"_scaleFactor"] = [x for x in AA[:,ip]]
        dd[pStrs[ip]+"_chisq"] = [x for x in chisq[:,ip]]

    df = pd.DataFrame(dd)

    df['mux'] = sim['mux']
    df['sigx'] = sim['sigx']
    df['theta'] = sim['theta']
    df['phi'] = sim['phi']
    df['diff'] = sim['diff']
    df['rst'] = sim['rst']
    
    df.to_pickle(fout)

    print(df)
