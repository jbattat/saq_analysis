import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import sys

from offset_simulation import compute_masks
"""
Steps as outlined by James
1. import the data and graph it
2. normalize it
  a. I think he said the easy way would be to normalize it based on the number of resets
3. generate a wide variety of offsets and diffusion widths
4. do a course matching (I think using chi**2)
5. do a finer matching (I think using chi**2)
"""
integral  = pickle.load(open('./test.pkl', 'rb'))
#gauss = pickle.load(open('./test.pkl', 'rb'))
print(sum(integral[0]))
print("integrals: ", integral[2])
#print("gaussians: ", gauss[0])
midpoint = ([ 0.255 ,  1.09  ,  1.855 ,  2.645 ,  3.605 ,  4.605 ,  5.605 ,
        6.855 ,  8.335 , 10.335 , 14.0475, 18.9625, 23.96  , 28.96  ,
       36.44  , 46.45  ])

area =  ([1.28679635e+00, 5.57132005e+00, 8.41638562e+00, 1.48188925e+01,
       2.28456618e+01, 2.87737686e+01, 3.57614560e+01, 6.37454282e+01,
       7.80786305e+01, 1.62677290e+02, 4.33812869e+02, 5.90985650e+02,
       7.52725600e+02, 9.16727496e+02, 2.27871834e+03, 2.86466762e+03])
file = open('long_scan/list.txt')
names = file.readlines()

v400_names = []
#### Load data into the file
for x,filename in enumerate(names):
    field = 400
    pressure = int(200 + (x*50))
    var_name ='p' + str(pressure) + '_v' + str(field)
    temp = pd.read_csv(f'long_scan/{filename[:19]}_data.csv')
    locals()[var_name] = temp
    v400_names.append(locals()[var_name])
number = len(names)
data = np.zeros([16, number])
sigma = np.zeros([16, number])

for x, name in enumerate(v400_names):
    data[:,x] = name.iloc[:,1]
    sigma[:,x] = name.iloc[:,2]
    scale = str(1 -  x/15)
    plt.plot(midpoint[1:], data[1:,x], label = 'Pressure: ' + str(200 + x*50) + ', Drift Field: 400 V', color = scale)
plt.xlim(0,20)
plt.legend()
plt.savefig("Drift_Field_400.pdf")
plt.clf()
sigma = np.sqrt(sigma)

def cost(data, model, sigma):
    difference = np.sum((data - model)/sigma)
    return difference


def chi_squared(data, model, sigma):

    chi_squared = np.sum(((data - model)/sigma)**2)
    
    return chi_squared


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


def course_match(raw, sigma_x, mu_x, sigma):
    area =  ([1.28679635e+00, 5.57132005e+00, 8.41638562e+00, 1.48188925e+01,
              2.28456618e+01, 2.87737686e+01, 3.57614560e+01, 6.37454282e+01,
              7.80786305e+01, 1.62677290e+02, 4.33812869e+02, 5.90985650e+02,
              7.52725600e+02, 9.16727496e+02, 2.27871834e+03, 2.86466762e+03])
    nn = 16
    chi = []
    sigma = sigma/area
    a_list = []
    
    for ii, diff in enumerate(sigma_x):
        for jj, offset in enumerate(mu_x):
#            irow = jj % len(mux)
#            icol = ii % len(sigx)
            num = jj + ii*7
            A = sum(integral[num][1:]*raw[1:])/sum(raw[1:]**2)
            raw = raw*A
            a_list.append(A)
            chi.append(cost(raw[1:14], integral[num][1:14], sigma[1:14]))
            #print(chi_squared(raw[1:14], integral[num][1:14], sigma[1:14]))
    fit = np.array(chi).argmin()
    scaling  = a_list[fit]
    sigma_value = int(fit / len(mu_x))
    mu_value    = int(fit % len(mu_x))
    
    return(sigma_x[sigma_value], mu_x[mu_value], fit, sigma, scaling)

for n in range(len(names)):
 
    sigma_x = np.linspace(1, 5, 20)

    mu_x    = np.linspace(0, 7, 20)

    parameters = course_match(data[:,n], sigma_x, mu_x, sigma[:,n])
 
    #print(parameters)
    sigma_2 = np.linspace(parameters[0] - 0.5, parameters[0] + 1, 4)
    mu_2 = np.linspace(parameters[1] - 1, parameters[1] + 1, 4)
    final_fit = course_match(data[:,n],sigma_2, mu_2, sigma[:,n])
#    print(final_fit)

    ll = np.zeros(16)
    masks, gg, area  = make_gaussian(parameters[0], parameters[1])
    for kk in range(16):
        ll[kk] = (np.sum(gg*masks[kk]))

        #normalize based on area of ring
    ll = ll/area
    ll = ll/sum(ll)
  #  print(ll)
     

    plt.plot(midpoint[1:], ll[1:]/max(ll[1:]), 'bo', label = 'model')
    #plt.plot(gg)
    y_error = sigma[1:,n]
    plt.plot(midpoint[1:], data_2[1:,n]/max(data_2[1:,n]), 'ko', label = 'data')
    #plt.errorbar(midpoint[2:], data[2:,n],
#             yerr = y_error,
#             fmt ='o')
    plt.xlim(0, 20)
    plt.legend()
    plt.show()


