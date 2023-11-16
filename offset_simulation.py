import sys
import itertools
import os.path
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
from scipy.optimize import curve_fit
import awkward as ak
import pandas as pd
import uproot

# gaussian parameters
#sigx = [1, 1.5, 2, 2.5, 3] # mm (diffusion width)
#sigx = [0.25, 0.5, 0.75, 1]
#sigx = [1.5, 1.6, 1.7]
sigx = [4.6, 4.6]
sigy = sigx
#mux = 0
mux = [3.3, 3.3]
#mux = np.array([0, 1, 2, 3, 4, 5])  # mm
#mux = [6, 6.8 , 7]
muy = 0  # mm

xmin = -50  # mm
xmax = 50
ymin = xmin
ymax = xmax

nx = 2000
ny = 2000

# Geometry of anode rings
chans = np.arange(16)

# Measured using GerbView by James Battat on 2023 June 8
# inner and outer radii (to midpoints of gaps between channels) in millimeters
annuli = [ [ 0.000,  0.670], [ 0.670,  1.386], [ 1.386,  2.106], [ 2.106,  2.981],
           [ 2.981,  4.005], [ 4.005,  4.994], [ 4.994,  6.015], [ 6.015,  7.481],
           [ 7.481,  9.994], [ 9.994, 12.497], [12.497, 15.023], [15.023, 19.996],
           [19.996, 24.962], [24.962, 30.026], [30.026, 39.977], [39.977, 50.065] ]

# min and max radii of each annulus
# OLD VALUES -- USE AT YOUR OWN RISK
#annuli = [ [ 0, 0.64], [ 0.64, 1.4775], [ 1.4775, 2.205], [ 2.205, 3.095],
#           [ 3.095, 4.105], [ 4.105, 5.1], [ 5.1, 6.115], [ 6.115, 7.595],
#           [ 7.595, 9.085], [ 9.0855, 11.590], [11.590, 16.505], [16.505, 21.460],
#           [21.460, 26.460], [26.460, 31.495], [31.495, 41.440], [41.440, 51.275] ]

def compute_areas():
    areas = np.zeros(16)
    for ii, chan in enumerate(chans):
        areas[chan]   = (annuli[ii][1]**2 - annuli[ii][0]**2)*np.pi
    return areas

def compute_midpoints():
    midpoints = np.zeros(16)
    for ii in range(16):
        #midpoints[ii] = (annuli[ii][1] - annuli[ii][0])/2 + annuli[ii][0]
        midpoints[ii] = 0.5*(annuli[ii][1] + annuli[ii][0])
    #print(f'midpoints = {midpoints}')
    return midpoints

def compute_masks2(xx, yy):
    nx = len(xx)
    ny = len(yy)
    masks = {}
    distSq = xx**2 + yy**2

    for ii, chan in enumerate(chans):
        #print(f"ii, chan = {ii}, {chan}")
        masks[chan] = np.zeros((ny, nx), dtype=int)
        rminSq = annuli[ii][0]**2
        rmaxSq = annuli[ii][1]**2
        #print(f"rmin, rmax = {np.sqrt(rminSq)}, {np.sqrt(rmaxSq)}")
        ids = np.where( ((distSq >= rminSq) & (distSq < rmaxSq)) )
        masks[chan][ids] = 1
        
    return masks

def make_diffusion_kernel(xx, yy, sigma, norm=1.0):
    # xx: meshgrid x values
    # yy: meshgrid y values
    # sigma: mm (diffusion amount from drift)
    return norm*np.exp( - 0.5*(xx**2 + yy**2)/sigma**2 )

def apply_drift_diffusion(gg, kernel):
    # standard convolution is ***way*** too slow (scipy.signal.convolve2d)
    return scipy.signal.fftconvolve(gg, kernel, mode='same')

def make_gaussian(xx, yy, sigx, mux, theta, phi):
    # theta [degrees]  is angle from the normal to the cathode
    # theta=0 gives symmetric gaussian
    # phi [deg.] angle of fiber in the cathode plane
    # phi=0 points to ... =pi/2 points to...
    
    sigy = sigx/np.cos(theta*np.pi/180.)
    muy = 0  # mm

    # see e.g. astropy documentation
    # https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Gaussian2D.html
    cpsq = np.cos(phi*np.pi/180.)**2 # cos squared
    spsq = np.sin(phi*np.pi/180.)**2 # sin squared
    s2p  = np.sin(2*phi*np.pi/180.)  # sin 2phi
    a = 0.5*cpsq/sigx**2 + 0.5*spsq/sigy**2
    b = 0.5*s2p/sigx**2 - 0.5*s2p/sigy**2
    c = 0.5*spsq/sigx**2 + 0.5*cpsq/sigy**2
    
    norm = 1
    #norm = 1/(2.0*np.pi *sigma**2)
    #gg = norm*np.exp( -(xx-mux)**2/(2*sigx**2) - (yy-muy)**2/(2*sigy**2) )
    gg = norm*np.exp( -a*(xx-mux)**2 - b*(xx-mux)*(yy-muy) -c*(yy-muy)**2 )
    return gg


def compute_masks(xx, yy):
    ##### DEFUNCT! see compute_masks2()
    nx = len(xx)
    ny = len(yy)
    # (inner,outer) in mm
#    annuli = [[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10],[10,11],[11,12],[12, 13], [13, 14], [14,15],[15,16]]
#    annuli = [ [ 0, 0.51], [ 0.77, 1.41], [ 1.535, 2.175], [ 2.235, 3.055],
#               [ 3.135, 4.075], [ 4.135, 5.075], [ 5.125, 6.085], [ 6.145, 7.565],
#               [ 7.625, 9.045], [ 9.125, 11.545], [11.635, 16.46], [16.55, 21.375],
#               [21.545, 26.375], [26.545, 31.375], [31.615, 41.265], [41.625, 51.275]
#               ]
    annuli = [ [ 0, 0.64], [ 0.64, 1.4775], [ 1.4775, 2.205], [ 2.205, 3.095],
               [ 3.095, 4.105], [ 4.105, 5.1], [ 5.1, 6.115], [ 6.115, 7.595],
               [ 7.595, 9.085], [ 9.0855, 11.590], [11.590, 16.505], [16.505, 21.460],
               [21.460, 26.460], [26.460, 31.495], [31.495, 41.440], [41.440, 51.275]
               ]
    nchans = 16
    chans = np.arange(nchans)
    masks = {}
    area = np.zeros(16)
    distSq = xx**2 + yy**2
    midpoint = np.zeros(16)
    print(f"nx = {nx}")
    print(f"ny = {ny}")
    for ii, chan in enumerate(chans):
        masks[chan] = np.zeros((ny, nx), dtype=int)
        rminSq = annuli[ii][0]**2
        rmaxSq = annuli[ii][1]**2
        area[chan]   = (annuli[ii][1]**2 - annuli[ii][0]**2)*np.pi
        midpoint[chan] = (annuli[ii][1] - annuli[ii][0])/2 + annuli[ii][0]
        ids = np.where( ((distSq >= rminSq) & (distSq < rmaxSq)) )
        masks[chan][ids] = 1
        
    return masks, area, midpoint

def gauss(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def chi_squared(data, model):

    chi_squared = np.sum((data - model)**2/model)
    
    return chi_squared

    
if __name__ == '__main__':
    
    xx, yy = np.meshgrid( np.linspace(xmin, xmax, nx),
                          np.linspace(ymin, ymax, ny))

    norm = 1
    #norm = 1/(2.0*np.pi *sigma**2)
    gg = np.exp( -(xx-mux[0])**2/(2*sigx[0]**2) - (yy-muy)**2/(2*sigy[0]**2) )


    masks, area, midpoint  = compute_masks(xx, yy)
    nn = len(masks)
   
    # loop over all masks and integrate
#    nn = len(masks)
#    yy = []
#    for ii in range(nn):
#        yy.append(np.sum(gg*masks[ii]))
#    plt.plot(yy, 'ko')
#    plt.show()
#    plt.savefig('fluxVsRadius.pdf')

#    plt.clf()

#    fig, axs = plt.subplots(4,4)
#    for ii in range(nn):
#        irow = ii//4
#        icol = ii % 4
#        axs[irow, icol].imshow(gg*masks[ii])
#    plt.show()
    

#    plt.imshow(gg)
    plt.show()
    data = [0.          ,0.         ,0.         ,0.         ,0.18001703 ,0.38114421
            ,0.86250856 ,1.         ,0.73741708 ,0.02527634 ,0.         ,0.00579908
            ,0.         ,0.         ,0.         ,0.        ]
    
    #make a graph with multiple subplots
    fig, axs = plt.subplots(len(mux), len(sigx), layout = 'constrained')

    #go through for each combination of mu and sigma to make the gaussian and then simulate the resets
    for ii, diff in enumerate(sigx):
        for jj, offset in enumerate(mux):
            ll = []
            gg = np.exp( -(xx-offset)**2/(2*diff**2) - (yy-offset)**2/(2*diff**2) )
            #determine which subplot to graph on 
            irow = jj % len(mux)
            icol = ii % len(sigx)

            #count the number of resets on each mask
            for kk in range(nn):
                ll.append(np.sum(gg*masks[kk]))

            #normalize based on area of ring
            ll = ll/area
            ll = ll/max(ll)
            print(ll)
            #chi_2 = chi_squared(data, ll)
            #print(chi_2)
            
            #fit each diffusion graph to a gaussian 
            mean = sum(midpoint*ll)/sum(ll)                  
            sigma = sum(ll*(midpoint-mean)**2)/sum(ll)

            popt, pcov = curve_fit(gauss, midpoint, ll, p0=[max(ll),mean,sigma])
            #print("diffusion = ", diff, " offset = ", offset, "center = ", popt[1], "sigma = ",  popt[2])

            #graph the model and the simulation
            xm = np.linspace(0, 20, 100)
            ym = gauss(xm, popt[0], popt[1], popt[2])
            
            axs[irow, icol].plot(midpoint, ll, 'o')
            axs[irow, icol].plot(xm, ym)
            axs[irow, icol].set_title(f'mean = {round(popt[1], 2)}, sigma = {round(popt[2], 2)}')
            axs[irow, icol].set_xlim(0, 20)
    plt.show()

def plot_pressure_scan(df, chi=None, radii=None, plot_uncal=True, plot_cal=False, plot_bestFit=False):
    # df as returned by load_pressure_data() (data)
    # chi a dataframe of chisq and template parameters
    npres = len(df)
    nrows = 3
    ncols = 4
    fig, axs = plt.subplots(nrows, ncols, figsize=(12,8))
    chans = np.arange(16)
    idx = np.arange(15)+1 # skip channel 1
    for ii in range(npres):
        irow = int(np.floor(ii/ncols))
        icol = ii % ncols
        pres = f"{df['pres'][ii]:.0f} Torr"
        if radii is not None:
            xvals = radii[:]
        else:
            xvals = chans[idx][:]
        if plot_uncal:
            axs[irow, icol].errorbar(xvals[idx], df['rst'][ii][idx], yerr=df['rstErr'][ii][idx], fmt='ko', label=pres)
            #axs[irow, icol].errorbar(chans[idx], df['rst'][ii][idx], yerr=df['rstErr'][ii][idx], fmt='ko', label=pres)
        if plot_cal:  #'rstCal' in df.columns:
            axs[irow, icol].plot(xvals[idx], df['rstCal'][ii][idx], 'bo', label=pres)
            
        #if 'bestFit' in df.columns:
        if plot_bestFit:
            axs[irow, icol].plot(xvals[idx], df['bestFit'][ii][idx], 'r:')
            
        if chi is None:
            pass
        else: 
            nchi  = chi.index
            for jj in nchi:
                axs[irow, icol].plot(chans, chi['rst'][jj], label = jj)
        axs[irow, icol].set_xlim([0,20])
        axs[irow, icol].legend()
    plt.show()

def load_pressure_data():
    # read all pressure data
    # return a single Pandas DataFrame
    # along with errors
    froot = 'data/pressure_scan/20230519/'
    fnames = ["05_19_2023_08_41_26_data.csv", "05_19_2023_09_31_24_data.csv",
              "05_19_2023_10_17_29_data.csv", "05_19_2023_11_03_17_data.csv",
              "05_19_2023_11_50_02_data.csv", "05_19_2023_12_36_24_data.csv",
              "05_19_2023_13_23_33_data.csv", "05_19_2023_14_10_43_data.csv",
              "05_19_2023_14_57_48_data.csv", "05_19_2023_15_45_23_data.csv",
              "05_19_2023_16_33_02_data.csv", "05_19_2023_17_21_28_data.csv"
              ]
    fnames = [os.path.join(froot, x) for x in fnames]
    pressures = np.array([200., 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750])

    # Read in the CSV files
    rsts = []
    rstsRaw = []
    for ff in fnames:
        chans, rst, rstRaw = np.loadtxt(ff, unpack=True, delimiter=',', skiprows=1)
        rsts.append(rst)
        rstsRaw.append(rstRaw)
    
    dd = {}
    dd['pres'] = pressures
    df = pd.DataFrame(dd)
    # enter the reset data into the dataframe as np arrays
    df['rst'] = [x for x in rsts]
    df['rstRaw'] = [x for x in rstsRaw]

    # invent errors
    #err = np.array([5, 0.25, 0.25, 0.25,
    #                2, 0.25, 1, 1,
    #                1, 1, 1, 1,
    #                1, 1, 1, 1])
    err = np.array([0.1]*16) # tiny errors...
    for ii in range(4,6):
        err[ii] = err[ii]/10
    errors = [err]*len(pressures)
    df['rstErr'] = [x for x in errors]
    return df
    
def read_pressure_data(fname, fmt='awkward'):
    if fmt.upper()=='CSV':
        df = load_pressure_data()
        return np.array(df['pres']), np.array(df['rst']), np.array(df['rstRaw'])
        
    tree = uproot.open(fname)["tree"]    
    pressures = tree['pressure'].array()
    rst = tree['rst'].array()
    rstRaw = tree['rstRaw'].array()

    if fmt == 'numpy':
        pressures = ak.to_numpy(pressures)
        rst = ak.to_numpy(rst)
        rstRaw = ak.to_numpy(rstRaw)
        
    return pressures, rst, rstRaw

def optimal_scaling_factor(model, data, errors=None, ids=None):
    if ids is None:
        ids = np.arange(len(data), dtype=int)
        
    if errors is None:
        return np.sum(data*model) / np.sum(model*model)
    else:
        return np.sum(data*model/errors**2) / np.sum(model*model/errors**2)

def ids_of_min_chisq(df, cols):
    # find the diffusion parameter that gives the minimum chi squared for each dataset
    # df is a subset of the whole simulation set for which the nuissance parameters are
    #    all the same (mux, sigx, theta, phi). But the diffusion parameter can vary
    #    across datasets
    ids = []
    for col in cols:
        ids.append(df[col].idxmin())
    return ids

def get_min_chisq(df, cols, ids):
    # get the chisq value of the best fit
    chisq = 0
    for ii, col in enumerate(cols):
        chisq += df[col][ids[ii]]
    return chisq

def get_best_diffs(df, cols, ids):
    # Return an array of diffusion values that give the best fit
    diffs = []
    for ii, idx in enumerate(ids):
        diffs.append(df['diff'][idx])
    return np.array(diffs)


##### DEFUNCT (not needed)
#def find_smallest_chisq_sum(df):
#    sums = []
#    bestCombo = None
#    minChi = 1e99
#    cols = [x for x in list(df.columns) if "_chisq" in x]
#    cols.sort()
#    for combination in itertools.product(range(len(df)), repeat=len(cols)):
#        sums.append(0)
#        print(combination)
#        for ii, colname in enumerate(cols):
#            sums[-1] += df[colname][df.index[combination[ii]]]
#        if sums[-1] < minChi:
#            minChi = sums[-1]
#            bestCombo = combination
#            #print(f"Found new winner: {bestCombo}, {minChi}")
#    return bestCombo, minChi
#
#    
#def compute_all_possible_sums(df):
#    sums = []
#    for combination in itertools.product(range(len(df)), repeat=len(df.columns)):
#        #print(f"combination = {combination}")
#        sums.append(0)
#        for ii, colname in enumerate(df.columns):
#            sums[-1] += df[colname][combination[ii]]
#    return sums
    
