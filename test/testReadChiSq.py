import pandas as pd


dfTempl = pd.read_pickle("../resets_templates.pkl")
dfChi   = pd.read_pickle("../resets_chisq.pkl")

print(dfTempl)
print()
print(dfChi)

