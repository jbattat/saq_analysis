import sys
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak
import uproot

if __name__ == '__main__':

    # get command line arguments...
    print(sys.argv)
    fname = sys.argv[1]

    tree = uproot.open(fname)["tree"]

    pp = tree.arrays(["mux", "sigx", "theta", "phi", "diff", "rst"])

    fig, axs = plt.subplots(3,3, figsize=(15,8))
    chans = np.arange(16, dtype=int)+1
    for ii in range(len(pp['mux'])):
        row = int(ii/3)
        col = ii % 3
        axs[row,col].plot(chans, pp['rst'][ii], 'ko')
        titlestr = "mux, sigx, theta, phi, diff = "
        titlestr += f"{pp['mux'][ii]}, {pp['sigx'][ii]}, {pp['theta'][ii]}, {pp['phi'][ii]}, {pp['diff'][ii]}"
        axs[row,col].set_title(titlestr)
    plt.show()
    
    # Data is stored as...
    #dd = {'mux':ak.Array(pUsed[:,0]),
    #      'sigx':ak.Array(pUsed[:,1]),
    #      'theta':ak.Array(pUsed[:,2]),
    #      'phi':ak.Array(pUsed[:,3]),
    #      'diff':ak.Array(pUsed[:,4]),
    #      'rst':ak.Array(rstData)
    #      }
    
    
