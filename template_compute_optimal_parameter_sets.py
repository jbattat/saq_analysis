# + Load simulation data and chisquared data
# + loop over each set of nuissance parameters and find the best combination of diffusions for that parameter set
# + Save/record the optimal sets

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import offset_simulation as off

# read in the chisq data (also contains simulation parameters and simulated rst values)
pklName = 'templates/resets_BIG_templates_chisq.pkl'
chi = pd.read_pickle(pklName)
print(chi)

outfile = pklName.replace(".pkl", "_best.pkl")

cols = [x for x in list(chi.columns) if "_chisq" in x and "NoErr" not in x]
cols.sort()  # ['200_chisq', '250_chisq', ..., '750_chisq']

colsNoErr = [x for x in list(chi.columns) if "_chisqNoErr" in x]
colsNoErr.sort()  # ['200_chisqNoErr', ..., '750_chisqNoErr']

out = {'mux':[],
       'sigx':[],
       'theta':[],
       'phi':[],
       'chisq':[],
       'chisqNoErr':[]
       }
diffs = []
diffsNoErr = []

if __name__ == '__main__':

    print(f"mux:   {chi['mux'].unique()}")
    print(f"sigx:  {chi['sigx'].unique()}")
    print(f"theta: {chi['theta'].unique()}")
    print(f"phi:   {chi['phi'].unique()}")
    print(f"diff:  {chi['diff'].unique()}")

    nsim = len(chi) # number of parameter sets

    print(f"mux =", end='', flush=True)
    for ii, mux in enumerate(chi['mux'].unique()):
        print(f" {mux}", end='', flush=True)
        for jj, sigx in enumerate(chi['sigx'].unique()):
            for kk, theta in enumerate(chi['theta'].unique()):
                for ww, phi in enumerate(chi['phi'].unique()):
                    df = chi.loc[ (chi['mux']==mux) & (chi['sigx']==sigx) &
                                  (chi['theta']==theta) & (chi['phi']==phi) ]
                    #print(df[cols])
                    ids = off.ids_of_min_chisq(df, cols)
                    bestChi = off.get_min_chisq(df, cols, ids)
                    diffs.append(off.get_best_diffs(df, cols, ids))

                    idsNoErr = off.ids_of_min_chisq(df, colsNoErr)
                    bestChiNoErr = off.get_min_chisq(df, colsNoErr, idsNoErr)
                    diffsNoErr.append(off.get_best_diffs(df, colsNoErr, idsNoErr))
                    
                    #print(ids)
                    #print(f"mux, sigx, theta, phi = {mux}, {sigx}, {theta}, {phi}")
                    #diffStr = ''
                    #for idx in ids:
                    #    diffStr += f" {df['diff'][idx]:.2f}"
                    #print(f"                 diff = "+diffStr)
                    #print(f"              bestChi = {bestChi}")
                    out['mux'].append(mux)
                    out['sigx'].append(sigx)
                    out['theta'].append(theta)
                    out['phi'].append(phi)
                    out['chisq'].append(bestChi)
                    out['chisqNoErr'].append(bestChiNoErr)
                    #sys.exit()
                    
    out['diffs'] = [x for x in diffs]
    out['diffsNoErr'] = [x for x in diffsNoErr]
    dfOut = pd.DataFrame(out)
    dfOut.to_pickle(outfile)

    # convert the rst data to numpy arrays in the dataframe
    #tmp = [np.array(x) for x in sim['rst']]
    #sim['rst'] = [x for x in tmp]
    
