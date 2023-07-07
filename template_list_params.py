# + Load simulation template data and print out parameter values
import pandas as pd

## Read in the template data
pklName = 'templates/resets_BIG_templates.pkl'
sim = pd.read_pickle(pklName)
print(f"mux:   {sim['mux'].unique()}")
print(f"sigx:  {sim['sigx'].unique()}")
print(f"theta: {sim['theta'].unique()}")
print(f"phi:   {sim['phi'].unique()}")
print(f"diff:  {sim['diff'].unique()}")
