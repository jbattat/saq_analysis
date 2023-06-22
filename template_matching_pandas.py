# load in the simulation data
# load in the measurements
# "fit" by scaling the measurements to best match the data
# save/record the scaled data and the chi-squared and the abs-chi in a "bdat" tree

import sys
import numpy as np
import pandas as pd
import ROOT
import uproot
import matplotlib.pyplot as plt
#import awkward as ak

data_fname = "data/pressure_scan/20230519/pressure_scan_data.root"
    
def optimal_scaling_factor(model, data, errors=None, ids=None):
    if ids is None:
        ids = np.arange(len(data), dtype=int)
        
    if errors is None:
        return np.sum(data*model) / np.sum(model*model)
    else:
        return np.sum(data*model/errors) / np.sum(model*model/errors)

def read_pressure_data(fname):
    tree = uproot.open(fname)["tree"]    
    pressures = tree['pressure'].array()
    rst = tree['rst'].array()
    rstRaw = tree['rstRaw'].array()

    return pressures, rst, rstRaw


def read_template_data(fnames):
    # fnames can be a string (single file) or a list of strings (multiples files)
    rdf = ROOT.RDataFrame("tree", fnames)

    # entries that are numpy arrays need special treatment
    # so first, get all the float/double entries
    cols = ['mux', 'sigx', 'theta', 'phi', 'diff']
    dd  = rdf.AsNumpy(cols)
    
    # Then handle the RVec entries separately
    dd2 = rdf.AsNumpy(['rst'])
    dd2 = [list(x) for x in dd2['rst']]
    #print(dd2)
    dd['rst'] = dd2
    
    return pd.DataFrame(dd)
    
if __name__ == '__main__':

    sim_fname = sys.argv[1:] # FIXME: add check for valid name
    print(f'sim_fname = {sim_fname}')
    foutChisq = 'resets_chisq.pkl'
    foutTemplates = 'resets_templates.pkl'
    print(f"chisq to be saved as pandas DataFrame in: {foutChisq}")
    print(f"merged templates to be saved as pandas DataFrame in: {foutTemplates}")
    
    # read in the pressure scan data
    pp, rst, rstRaw = read_pressure_data(data_fname)
    presStrings = [f'{xx:.1f}' for xx in pp]  # used for Pandas column names

    # read in the template data
    sim = read_template_data(sim_fname)
    #print(sim)
    #print(sim['rst'][0])
    #sys.exit()
    
    chans = np.arange(16)+1

    # skip channel 1
    ids = np.arange(1, 15, dtype=int)  
    
    nd = len(pp)         # number of datasets (pressures)
    ns = len(sim.index)  # number of simulations
    AA = np.zeros((ns, nd), dtype=float)
    chisq = np.zeros((ns, nd), dtype=float)
    for ip in range(nd):  # loop over pressures (data)
        pres = pp[ip]
        print('pres = '+presStrings[ip])
        for ii in range(ns): # loop over simulation parameters
            templ = np.array(sim['rst'][ii])
            AA[ii,ip] = optimal_scaling_factor(templ, rst[ip], errors=None, ids=ids)
            chisq[ii,ip] = np.sum( (rst[ip] - templ*AA[ii,ip])**2 )
            #AA[ii,ip] = optimal_scaling_factor(sim['rst'][ii], rst[ip], errors=None, ids=ids)
            #chisq[ii,ip] = np.sum( (rst[ip] - sim['rst'][ii]*AA[ii,ip])**2 )

    ## compile the data
    dd = {}

    # each pressure gets two columns (scale factor and chisquared value)
    for ip, ps in enumerate(presStrings):
        dd[ps+"_scaleFactor"] = AA[:,ip]
        dd[ps+"_chisq"] = chisq[:,ip]
    
    df2 = pd.DataFrame.from_dict(dd)

    #print(sim)
    #print(df2)
    df2.to_pickle(foutChisq)
    sim.to_pickle(foutTemplates)
