import pandas as pd
import numpy as np
import saq
import matplotlib.pyplot as plt
import sys
import pickle

midpoint = saq.midpoint
area     = saq.area

resetdata = sys.argv[1]
integral  = pd.DataFrame(pickle.load(open('../template_50.pkl', 'rb')))
chan, resetperarea, reset = np.loadtxt(resetdata, unpack=True, skiprows=1, delimiter=',')

diffusion  = float(input("Enter a diffusion between 0 and 5:"))
offset     = float(input("Enter an offset between 0 and 7: "))

sigma = integral.diffusion.unique()
mu    = integral.offset.unique()
print(sigma)
print(mu)

integral_ind = integral.set_index(["diffusion", "offset"])

def find_parameter(param, var):
    difference = np.zeros(len(param))
    for ii,x in enumerate(param):
        difference[ii] = (abs(var - x))
    return difference.argmin()    

sigma_1 = sigma[find_parameter(sigma, diffusion)]
mu_1    = mu[find_parameter(mu, offset)]
#print(sigma_1, mu_1)

model = integral_ind.loc[sigma_1, mu_1]
print(model[0])
print(resetperarea)
A     = sum(resetperarea[1:10])/sum(model[0][1:10])

plt.plot(midpoint[1:], resetperarea[1:], 'ko', label = 'Data')
plt.plot(midpoint[1:], model[0][1:]*A    , 'bo', label = 'Model')
plt.legend()
plt.show()
