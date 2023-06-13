import numpy as np
import ROOT
import ROOT.RDF
import awkward as ak
import uproot

nn = 16 # number of SAQ channels

def makeTFileND():
    xx = np.array([1, 2], dtype=np.int32)
    yy = np.array([[10, 11, 12], [20, 21, 22]], dtype=np.float64)

    dd = {'xx':ak.Array(xx),
          'yy':ak.Array(yy)}

    ff = uproot.recreate("junkND.root")
    ff['tree'] = dd

def readTFileND(fname):
    tt = uproot.open("junkND.root:tree")
    print(tt)
    print(tt.keys())
    xx = tt['xx'].array(library='np')
    yy = tt['yy'].array(library='np')
    print(f'xx = {xx}')
    print(f'yy = {yy}')

    
#def readTFile(fname):
#    df = ROOT.RDataFrame("tree", fname)
#    #cols = df.Filter("xx >= 2").AsNumpy(["xx", "yy"]) # retrieve columns "x" and "y" as NumPy arrays
#    print(cols['xx'])
#    print(cols['yy'])
#    print(type(cols))
     
    
def makeTFile1D():
    # This only works for data that is "1D" (one entry per row)
    
    # Let's create some data in numpy arrays
    xx = np.array([1, 2, 3], dtype=np.int32)
    yy = np.array([4, 5, 6], dtype=np.float64)
    
    # Read the data with RDataFrame
    # The column names in the RDataFrame are defined by the keys of the dictionary.
    # Please note that only fundamental types (int, float, ...) are supported and
    # the arrays must have the same length.
    df = ROOT.RDF.MakeNumpyDataFrame({'xx': xx, 'yy': yy})

    # Save the data as a ROOT file
    df.Snapshot('tree', 'junk1D.root')

def readTFile1D(fname):
    df = ROOT.RDataFrame("tree", fname)
    cols = df.Filter("xx >= 2").AsNumpy(["xx", "yy"]) # retrieve columns "x" and "y" as NumPy arrays
    print(cols['xx'])
    print(cols['yy'])
    print(type(cols))
    
if __name__ == '__main__':
    print(f"ROOT.__version__ = {ROOT.__version__}")
    #makeTFile1D()
    #makeTFileND()
    #readTFile1D("junk1D.root")
    readTFileND("junkND.root")
