import numpy as np
import awkward as ak
import uproot
import matplotlib.pyplot as plt

if __name__ == '__main__':

    tree = uproot.open("data/pressure_scan/20230519/pressure_scan_data.root")["tree"]
    pressures = tree['pressure'].array()
    #rsts, rstRaw = tree.arrays(['rst', 'rstRaw'])
    rst = tree['rst'].array()
    rstRaw = tree['rstRaw'].array()
    
    print(pressures)
    print(rst)
    print(rst[0])

    npres = len(pressures)

    for ii in range(npres):
        #plt.plot(rst[ii]/np.max(rst[ii]), 'k.-')
        plt.plot(rstRaw[ii], 'k.-')

    plt.show()
