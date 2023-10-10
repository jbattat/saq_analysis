import ROOT

files = ["data/test_cal/10_10_2023_11_19_26.root",
         "data/test_cal/10_10_2023_11_21_29.root"]
         
chain = ROOT.TChain("tt")
for ff in files:
    chain.AddFile(ff)

chain.Scan()
