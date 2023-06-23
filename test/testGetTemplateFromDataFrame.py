# load in the simulated templates in (multiple) ROOT files
# output a single .pkl file containing a Pandas dataframe with all templates

import sys
import os.path
import pandas as pd
import numpy as np

if __name__ == '__main__':

    pklName = '../templates/resets_BIG_templates.pkl'
    
    df = pd.read_pickle(pklName)
    print(df)

    mux = 4.7
    sigx = 2.0
    theta = 80.0
    phi = 87.0
    diff = 4.0

    rst = np.array(df.loc[(df['mux'] == mux) &
                          (df['sigx'] == sigx) &
                          (df['theta'] == theta) &
                          (df['phi'] == phi) &
                          (df['diff'] == diff) ]['rst'].values[0])

    print(f"mux:   {df['mux'].unique()}")
    print(f"sigx:  {df['sigx'].unique()}")
    print(f"theta: {df['theta'].unique()}")
    print(f"phi:   {df['phi'].unique()}")
    print(f"diff:  {df['diff'].unique()}")
    print()
    print(f'rst = {rst}')
    
    #print(type(a['rst']))
    #print()
    #print(a['rst'])
    #print()
    #print(np.array(a['rst']))
    #print()
    #print(a['rst'].to_numpy())
    #print()
    #print(np.array(a['rst'].values[0]))

