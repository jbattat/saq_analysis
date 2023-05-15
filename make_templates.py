import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle 
from offset_simulation import compute_masks


midpoint = ([ 0.255 ,  1.09  ,  1.855 ,  2.645 ,  3.605 ,  4.605 ,  5.605 ,
        6.855 ,  8.335 , 10.335 , 14.0475, 18.9625, 23.96  , 28.96  ,
       36.44  , 46.45  ])

area =  ([1.28679635e+00, 5.57132005e+00, 8.41638562e+00, 1.48188925e+01,
       2.28456618e+01, 2.87737686e+01, 3.57614560e+01, 6.37454282e+01,
       7.80786305e+01, 1.62677290e+02, 4.33812869e+02, 5.90985650e+02,
       7.52725600e+02, 9.16727496e+02, 2.27871834e+03, 2.86466762e+03])

def make_gaussian(sigx, mux):
    sigy = sigx
    muy = 0  # mm
    
    xmin = -10  # mm
    xmax = 10
    ymin = xmin
    ymax = xmax
    
    nx = 2000
    ny = 2000
    xx, yy = np.meshgrid( np.linspace(xmin, xmax, nx),
                          np.linspace(ymin, ymax, ny))
    masks, area, midpoint  = compute_masks(xx, yy)

    norm = 1
    #norm = 1/(2.0*np.pi *sigma**2)
    gg = np.exp( -(xx-mux)**2/(2*sigx**2) - (yy-muy)**2/(2*sigy**2) )

    return(masks, gg, area)

    
nn = 16
sigma_x = np.linspace(1, 5, 20)
print(sigma_x)
mu_x    = np.linspace(0, 7, 20)
print(mu_x)
gauss_templates = []
int_templates = []
for ii, diff in enumerate(sigma_x):
    for jj, offset in enumerate(mu_x):
        ll = []
        masks, gg, area = make_gaussian(diff, offset)
        gauss_templates.append(gg)

        #determine which subplot to graph on 
        #            irow = jj % len(mux)
        #            icol = ii % len(sigx)
        
        #count the number of resets on each mask
        for kk in range(nn):
            ll.append(np.sum(gg*masks[kk]))
            
            #normalize based on area of ring
        ll = ll/area
        ll = ll/sum(ll)
        int_templates.append(ll)
        int_name = str(diff) + "_" + str(offset) + "_int"
pickle.dump(int_templates, open('test.pkl', 'wb'))
#pickle.dump(gauss_templates, open('test.pkl', 'wb'))
print("pickled") 





