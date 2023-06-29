# load in the simulated templates in (multiple) ROOT files
# output a single .pkl file containing a Pandas dataframe with all templates
# Reset data is stored in the pandas dataframe as numpy arrays (not lists)

import sys
import os.path
import numpy as np
import pandas as pd
import ROOT

# FIXME: move this to offset_simulation.py ?
def read_template_data(fnames):
    # fnames can be a string (single file) or a list of strings (multiples files)
    rdf = ROOT.RDataFrame("tree", fnames)

    # entries that are numpy arrays need special treatment
    # so first, get all the float/double entries
    cols = ['mux', 'sigx', 'theta', 'phi', 'diff']
    dd  = rdf.AsNumpy(cols)
    
    # Then handle the RVec entries separately
    dd2 = rdf.AsNumpy(['rst'])
    dd['rst'] = [list(x) for x in dd2['rst']]
    #dd['rst'] = [x for x in dd2['rst']] # should be ok but fails/hangs sometimes. Why???
    # Oh well, can convert from list() to np.array() later...
    
    return pd.DataFrame(dd)
    
if __name__ == '__main__':

    #sim_fname = sys.argv[1:] # FIXME: add check for valid name

    # FIXME: allow users to specify the root files
    sim_fname = ['resets_20230622110548.root', 'resets_20230622110555.root',
                 'resets_20230622110601.root', 'resets_20230622110607.root',
                 'resets_20230622110613.root', 'resets_20230622110619.root',
                 'resets_20230622110623.root', 'resets_20230622110627.root',
                 'resets_20230622110631.root', 'resets_20230622110637.root', 
                 'resets_20230622110643.root', 'resets_20230622110656.root',
                 'resets_20230622110700.root', 'resets_20230622110707.root',
                 'resets_20230622110713.root', 'resets_20230622110716.root']
    #sim_fname = ['resets_20230622110548.root']

    subdirName = './templates'
    sim_fname = [os.path.join(subdirName, ff) for ff in sim_fname]
    print(f'sim_fname = {sim_fname}')

    # read in the template data
    print("Reading template data")
    sim = read_template_data(sim_fname)

    # convert the rst data to numpy arrays in the dataframe
    tmp = [np.array(x) for x in sim['rst']]
    sim['rst'] = [x for x in tmp]

    # Pickle the DataFrame to disk
    foutTemplates = 'resets_BIG_templates.pkl'
    print(f"Saving merged template data as Pandas dataFrame in {foutTemplates}")
    sim.to_pickle(foutTemplates)
