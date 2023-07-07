import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import offset_simulation as off

# read in the pressure scan data
data = off.load_pressure_data()
pres = data['pres'].to_numpy()
rstData = data['rst'].to_numpy()
err = data['rstErr'].to_numpy()
pStrs = [f'{p:.0f}' for p in pres]
#pres, rstData, rstDataRaw = off.read_pressure_data("data/pressure_scan/20230519/pressure_scan_data.root", fmt='numpy')

# read in the chisq data (also contains simulation parameters and simulated rst values)
pklName = 'templates/resets_BIG_templates_chisq.pkl'
sim = pd.read_pickle(pklName)
print(sim.columns)

# read in the chisq data (also contains simulation parameters and simulated rst values)
pklName = 'templates/resets_BIG_templates_chisq_best.pkl'
chi = pd.read_pickle(pklName)
#print(chi)
print(chi.columns)

print(f"min chisq = {chi['chisq'].min()}")
idx = chi['chisq'].idxmin()
idxNoErr = chi['chisqNoErr'].idxmin()
print(chi.loc[idx])
diffs = chi['diffs'][idx]
print(diffs)
#print(chi['diffs'][idx])
print()
print(f"min chisqNoErr = {chi['chisqNoErr'].min()}")
print(chi.loc[idxNoErr])
diffsNoErr = chi['diffsNoErr'][idxNoErr]
print(diffsNoErr)

# get the values of the nuissance parameters
mux   = chi['mux'][idx]
sigx  = chi['sigx'][idx]
theta = chi['theta'][idx]
phi   = chi['phi'][idx]
# extract the rows of the full simulation set that have those nuissance parameters
df = sim.loc[ (sim['mux']==mux) & (sim['sigx']==sigx) &
              (sim['theta']==theta) & (sim['phi']==phi) ]
print(df)
print()

# get the values of the nuissance parameters
muxNoErr   = chi['mux'][idxNoErr]
sigxNoErr  = chi['sigx'][idxNoErr]
thetaNoErr = chi['theta'][idxNoErr]
phiNoErr   = chi['phi'][idxNoErr]
# extract the rows of the full simulation set that have those nuissance parameters
dfNoErr = sim.loc[ (sim['mux']==muxNoErr) & (sim['sigx']==sigxNoErr) &
                   (sim['theta']==thetaNoErr) & (sim['phi']==phiNoErr) ]
print(dfNoErr)

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
    #axs[row, col].plot(chan[ids], rstData[ii][ids], 'ko', label=f"{pp:.0f} Torr")
    axs[row, col].errorbar(chan[ids], rstData[ii][ids], yerr=err[ii][ids], fmt='ko', label=f"{pp:.0f} Torr")
    axs[row, col].plot(chan[ids], rstSim[ii][ids]*scale[ii], 'r:', label=f"{diffs[ii]:.2f}mm")
    #titleStr = f"{pp} Torr; {diffs[ii]:.2f}"
    #axs[row, col].set_title(titleStr)
    if ii == 0:
        axs[row, col].set_title(f"$\mu_x$, $\sigma_x$, $\\theta$, $\phi$ = {mux}mm, {sigx}mm, {theta}$^o$, {phi}$^o$")
        #axs[row, col].set_xlabel(f"$\mu_x$, $\sigma_x$, $\\theta$, $\phi$ = {mux}, {sigx}, {theta}, {phi}")
    axs[row, col].legend()

outfile = 'junk.pdf'
plt.savefig(outfile)
plt.savefig(outfile.replace(".pdf", ".png"))

plt.clf()
#chi.hist(column='chisq', bins=100)
#bins = 200
bins = np.linspace(chi['chisq'].min(), 200, num=100)
chi['chisq'].plot(kind='hist', bins=bins, logy=True, xlabel='Chisq', ylabel='N')
plt.savefig("junk2.pdf")


#def plot_pressure_scan(df, chi=None):
#    # df as returned by load_pressure_data() (data)
#    # chi a dataframe of chisq and template parameters
#    npres = len(df)
#    nrows = 3
#    ncols = 4
#    fig, axs = plt.subplots(nrows, ncols, figsize=(12,8))
#    chans = np.arange(16)
#    idx = np.arange(15)+1 # skip channel 1
#    for ii in range(npres):
#        irow = int(np.floor(ii/ncols))
#        icol = ii % ncols
#        pres = f"{df['pres'][ii]:.0f} Torr"
#        axs[irow, icol].errorbar(chans[idx], df['rst'][ii][idx], yerr=df['rstErr'][ii][idx], fmt='ko', label=pres)
#        axs[irow, icol].legend()
#    plt.show()
#
