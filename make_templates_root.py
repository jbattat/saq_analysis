import sys
import argparse
from datetime import datetime
import numpy as np
import awkward as ak
import uproot
import offset_simulation as off

parser = argparse.ArgumentParser()
parser.add_argument('-m', '--mux', nargs='+')    # offset of UV illumination from center of anode
parser.add_argument('-t', '--theta', nargs='+')  # fiber angle perp to cathode
parser.add_argument('-p', '--phi', nargs='+')    # fiber angle in plane of cathode
parser.add_argument('-s', '--sigx', nargs='+')   # width of UV illumination
parser.add_argument('-d', '--diff', nargs='+')   # diffusion due to drift (mm)
parser.add_argument('-f', '--fout', nargs=1, help='Name of output file')   # diffusion due to drift (mm)
parser.add_argument('-n', '--dry-run', action='store_true', help='Show what would have been done')   # compute how many samples...

# --mux 0 5 0.1
# --theta 60 80 2
# --phi -87 87 3
# --sigx 1.0 2.0 0.5
# --diff 0 4.0 0.2

def makeParamRange(pp):
    # Generate the parameters to scan over
    #mux_min, mux_max, mux_step = float(args.mux[0]), float(args.mux[1]), float(args.mux[2])

    if len(pp) == 1:
        return np.array(pp, dtype=float)

    pmin, pmax, pstep = float(pp[0]), float(pp[1]), float(pp[2])
    nn = int( (pmax-pmin)/pstep ) + 1
    #print(f'nn = {nn}')
    pvals, pstep2 = np.linspace(pmin, pmax, num=nn, endpoint=True, retstep=True)
    if pstep2 != pstep:
        print("Warning: requested and implemented parameter step sizes don't agree")
        print(f"   step requested, step used = {pstep}, {pstep2}")
    return pvals

def makeOutFileName(fout):
    if fout is None:
        timestr = datetime.now().strftime("%Y%m%d%H%M%S")
        fout = f'resets_'+timestr+'.root'
    else:
        fout = fout[0]
    return fout

if __name__ == '__main__':

    # get command line arguments...
    args = parser.parse_args()
    print(args)

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)    
    # FIXME: Check for valid inputs...
    # either one element or three. And valid floats...

    # Generate values to step over
    muxs  = makeParamRange(args.mux)
    print(f"muxs = {muxs}")
    sigxs  = makeParamRange(args.sigx)
    print(f"sigxs = {sigxs}")
    thetas = makeParamRange(args.theta)
    print(f"thetas = {thetas}")
    phis = makeParamRange(args.phi)
    print(f"phis = {phis}")
    diffs = makeParamRange(args.diff)
    print(f"diffs = {diffs}")
    fout = makeOutFileName(args.fout)
    print(f"fout = {fout}")

    nmux = len(muxs)
    nsigx = len(sigxs)
    ntheta = len(thetas)
    nphi = len(phis)
    ndiff = len(diffs)
    ntot = nmux*nsigx*ntheta*nphi*ndiff
    print(f"nmux, nsigx, ntheta, nphi, ndiff = {nmux}, {nsigx}, {ntheta}, {nphi}, {ndiff}")
    print(f"total number of templates: {ntot}")
    if args.dry_run:
        sys.exit()

    nn = 16  # Number of SAQ readout channels
    
    # The "active area"
    xmin = -50.0 # mm
    xmax = -xmin
    ymin = xmin
    ymax = xmax
    nx = 2000
    ny = nx
    
    xx, yy = np.meshgrid( np.linspace(xmin, xmax, nx),
                          np.linspace(ymin, ymax, ny))
    masks = off.compute_masks2(xx, yy)

    rstData = np.zeros((ntot, 16), dtype=float)

    pUsed = np.zeros( (ntot, 5), dtype=float)
    
    #print("rstData = ")
    #print(rstData)


    areas = off.compute_areas()
    #sys.exit()
    convolve = False
    row = 0
    interval = 1000
    print("Starting iteration:")
    for mm in range(ndiff):
        if diffs[mm] != 0:
            convolve = True
            kernel = off.make_diffusion_kernel(xx, yy, diffs[mm], norm=1.0)
        for jj in range(nsigx):
            for kk in range(ntheta):
                for ll in range(nphi):
                    for ii in range(nmux):
                        #print(ii, jj, kk, ll, mm)
                        if row % interval == 0:
                            print(f"           {row}/{ntot}")
                        pUsed[row] = [muxs[ii], sigxs[jj], thetas[kk], phis[ll], diffs[mm]]
                        gg = off.make_gaussian(xx, yy, sigxs[jj], muxs[ii], thetas[kk], phis[ll])
                        if convolve:
                            gg = off.apply_drift_diffusion(gg, kernel)
                        for chan in range(nn):
                            rstData[row,chan] = np.sum(gg*masks[chan])
                        rstData[row] /= areas # normalize based on area of ring
                        rstData[row] /= np.max(rstData[row]) # normalize...
                        row += 1
    #print("rstData = ")
    #print(rstData)
                        
    dd = {'mux':ak.Array(pUsed[:,0]),
          'sigx':ak.Array(pUsed[:,1]),
          'theta':ak.Array(pUsed[:,2]),
          'phi':ak.Array(pUsed[:,3]),
          'diff':ak.Array(pUsed[:,4]),
          'rst':ak.Array(rstData)
          }
    
    ff = uproot.recreate(fout)
    ff['tree'] = dd
    
