import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import offset_simulation as off

# read in the pressure scan data
pres, rstData, rstDataRaw = off.read_pressure_data("data/pressure_scan/20230519/pressure_scan_data.root", fmt='numpy')
pStrs = [f'{p:.0f}' for p in pres]

# read in the chisq data (also contains simulation parameters and simulated rst values)
pklName = 'templates/resets_BIG_templates_chisq.pkl'
sim = pd.read_pickle(pklName)
print(sim.columns)

# read in the chisq data (also contains simulation parameters and simulated rst values)
pklName = 'templates/resets_BIG_templates_chisq_best.pkl'
chi = pd.read_pickle(pklName)
print(chi)

print(f"min chisq = {chi['chisq'].min()}")
idx = chi['chisq'].idxmin()
print(chi.loc[idx])
#print(chi['diffs'][idx])
diffs = chi['diffs'][idx]
print(diffs)

# get the values of the nuissance parameters
mux   = chi['mux'][idx]
sigx  = chi['sigx'][idx]
theta = chi['theta'][idx]
phi   = chi['phi'][idx]
# extract the rows of the full simulation set that have those nuissance parameters
df = sim.loc[ (sim['mux']==mux) & (sim['sigx']==sigx) &
              (sim['theta']==theta) & (sim['phi']==phi) ]
print(df)
rstSim = []
scale = []
for ii, diff in enumerate(chi['diffs'][idx]):
    rstSim.append(df.loc[df['diff'] == diff, 'rst'])  # FIXME: extract values,  not series...?
    colName = f"{pres[ii]:.0f}_scaleFactor"
    print(colName)
    scale.append(df.loc[df['diff'] == diff, colName])

# convert to arrays (FIXME could do this directly in the loop above...)
rstSim = [x.values[0] for x in rstSim]
scale = [x.values[0] for x in scale]

chan = np.arange(16)
ids = np.arange(1, 16) # ignore first channel
fig, axs = plt.subplots(3,4, figsize=(12,8))
for ii, pp in enumerate(pres):
    row = int(ii/4)
    col = ii % 4
    #print(f"ii, row, col = {ii}, {row}, {col}")
    axs[row, col].plot(chan[ids], rstData[ii][ids], 'ko')
    axs[row, col].plot(chan[ids], rstSim[ii][ids]*scale[ii], 'r:')
    titleStr = f"{pp} Torr; {diffs[ii]:.2f}"
    axs[row, col].set_title(titleStr)
plt.savefig("junk.pdf")
