import os
import sys
import numpy as np
import awkward as ak
import uproot

froot = 'data/pressure_scan/20230519/'
fnames = ["05_19_2023_08_41_26_data.csv", "05_19_2023_09_31_24_data.csv",
          "05_19_2023_10_17_29_data.csv", "05_19_2023_11_03_17_data.csv",
          "05_19_2023_11_50_02_data.csv", "05_19_2023_12_36_24_data.csv",
          "05_19_2023_13_23_33_data.csv", "05_19_2023_14_10_43_data.csv",
          "05_19_2023_14_57_48_data.csv", "05_19_2023_15_45_23_data.csv",
          "05_19_2023_16_33_02_data.csv", "05_19_2023_17_21_28_data.csv"
          ]
fnames = [os.path.join(froot,x) for x in fnames]
pressures = np.array([200., 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750])

if __name__ == '__main__':

    nn = 16  # Number of SAQ readout channels
    nfiles = len(fnames)
    data = np.zeros( (nfiles, 16), dtype=float)
    dataRaw = np.zeros( (nfiles, 16), dtype=float)
    for ii, fname in enumerate(fnames):
        _, data[ii], dataRaw[ii] = np.loadtxt(fname, delimiter=",", unpack=True, skiprows=1)
        
    dd = {'pressure':pressures,
          'rst':data,
          'rstRaw':dataRaw}
    ff = uproot.recreate(os.path.join(froot,"pressure_scan_data.root"))
    ff['tree'] = dd
