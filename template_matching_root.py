# load in the simulation data
# load in the measurements
# "fit" by scaling the measurements to best match the data
# save/record the scaled data and the chi-squared and the abs-chi in a "bdat" tree

import numpy as np
import uproot
import matplotlib.pyplot as plt
import awkward as ak

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


def read_template_data(fname):
    tree = uproot.open(fname)["tree"]
    dd = tree.arrays(["mux", "sigx", "theta", "phi", "diff", "rst"])
    return dd

    #"{pp['mux'][ii]}, {pp['sigx'][ii]}, {pp['theta'][ii]}, {pp['phi'][ii]}, {pp['diff'][ii]}"
    
if __name__ == '__main__':

    data_fname = "data/pressure_scan/20230519/pressure_scan_data.root"
    
    # read in the pressure scan data
    pp, rst, rstRaw = read_pressure_data(data_fname)
    presStrings = [f'{xx:.1f}' for xx in pp]  # used for TTree branch names?

    # template data
    #fname = 'resets_9.root'
    sim_fname = 'resets_small.root'
    bdat_fname = sim_fname.replace(".root", "_bdat.root")
    sim = read_template_data(sim_fname)

    chans = np.arange(16)+1

    ids = np.arange(15, dtype=int)+1 # skip channel 1
    
    nd = len(pp)         # number of datasets (pressures)
    ns = len(sim['mux'])  # number of simulations
    AA = np.zeros((ns, nd), dtype=float)
    chisq = np.zeros((ns, nd), dtype=float)
    for ip in range(nd):  # loop over pressures (data)
        pres = pp[ip]
        print('pres = '+presStrings[ip])
        for ii in range(ns): # loop over simulation parameters
            AA[ii,ip] = optimal_scaling_factor(sim['rst'][ii], rst[ip], errors=None, ids=ids)
            chisq[ii,ip] = np.sum( (rst[ip] - sim['rst'][ii]*AA[ii,ip])**2 )

    # compile the TTree data
    dd = {'mux':ak.Array(sim['mux']),
          'sigx':ak.Array(sim['sigx']),
          'theta':ak.Array(sim['theta']),
          'phi':ak.Array(sim['phi']),
          'diff':ak.Array(sim['diff'])
          }
    # each pressure gets two columns (scale factor and chisquared value)
    for ip, ps in enumerate(presStrings):
        dd[ps+"_scaleFactor"] = ak.Array(AA[:,ip])
        dd[ps+"_chisq"] = ak.Array(chisq[:,ip])
    
    ff = uproot.recreate(bdat_fname)
    ff['bdat'] = dd
    
    # Plot a single pressure data set with a bunch of templates overlaid
    #ip = 7
    #plt.plot(chans[ids], rst[ip][ids], 'ko')
    #colors = np.linspace(0.2,1, num=ns)
    #for ii in range(ns):
    #    plt.plot(chans[ids], sim['rst'][ii][ids]*AA[ii,ip], color=str(colors[ii]), linestyle=":")
    #plt.title(f"Pressure = {pp[ip]} Torr")    
    #plt.show()
