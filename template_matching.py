# 16 --> saq.N_SAQ_CHANNELS

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import sys
import saq

from offset_simulation import compute_masks


integral  = pd.DataFrame(pickle.load(open('../template_20.pkl', 'rb')))
sigma_x   = integral.diffusion.unique()
mu_x      = integral.offset.unique()
integral_ind  = integral.set_index(['diffusion', 'offset'])
print(sigma_x)
print(mu_x)



midpoint = saq.midpoint[:]
area = saq.area[:]

#file = open('long_scan/list.txt')
#names = file.readlines()
names  = saq.files_from_list('list.txt')
print(names)

#str.replace(".root", "_myplot.pdf")

#v400_names = []
##### Load data into the file
#for x,filename in enumerate(names):
#    field = 400
#    pressure = int(200 + (x*50))
#    var_name ='p' + str(pressure) + '_v' + str(field)
#    resetdata = filename.replace(".root", "_data.csv")
#    temp = pd.read_csv(resetdata)
#    #temp = pd.read_csv(f'long_scan/{filename[:19]}_data.csv')
#    locals()[var_name] = temp
#    v400_names.append(locals()[var_name])
#number = len(names)
#data = np.zeros([16, number])
#sigma = np.zeros([16, number])

#### Load data into the file
number = len(names)
data = np.zeros([saq.N_SAQ_CHANNELS, number])
sigma = np.zeros([saq.N_SAQ_CHANNELS, number])
for x,filename in enumerate(names):
    field = 400
    pressure = int(200 + (x*50))
    #var_name ='p' + str(pressure) + '_v' + str(field)
    resetdata = filename.replace(".root", "_data.csv")
    #temp = pd.read_csv(resetdata)
    chan, resetperarea, reset = np.loadtxt(resetdata, unpack=True, skiprows=1, delimiter=',')
    #temp = pd.read_csv(f'long_scan/{filename[:19]}_data.csv')
    #locals()[var_name] = temp
    #v400_names.append(locals()[var_name])
    data[:,x] = resetperarea[:]
    sigma[:,x] = np.sqrt(reset[:])/area
    scale = str(1 -  x/15)
    plt.plot(midpoint[1:], data[1:,x], label = 'Pressure: ' + str(200 + x*50) + ', Drift Field: 400 V', color = scale)

plt.xlim(0,20)
plt.legend()
plt.savefig("Drift_Field_400.pdf")
plt.clf()

#for x, name in enumerate(v400_names):
#    data[:,x] = name.iloc[:,1]
#    sigma[:,x] = name.iloc[:,2]
#    scale = str(1 -  x/15)
#    plt.plot(midpoint[1:], data[1:,x], label = 'Pressure: ' + str(200 + x*50) + ', Drift Field: 400 V', color = scale)
#plt.xlim(0,20)
#plt.legend()
#plt.savefig("Drift_Field_400.pdf")
#plt.clf()
#sigma = np.sqrt(sigma)

def cost(data, model, sigma):
    difference = np.sum(np.abs(data - model)/sigma)
    return difference


def chi_squared(data, model, sigma):

    chi_squared = np.sum(((data - model)/sigma)**2)
    
    return chi_squared


#def make_gaussian(sigx, mux):
#    sigy = sigx
#    muy = 0  # mm
#    
#    xmin = -10  # mm
#    xmax = 10
#    ymin = xmin
#    ymax = xmax
#    
#    nx = 2000
#    ny = 2000
#    xx, yy = np.meshgrid( np.linspace(xmin, xmax, nx),
#                          np.linspace(ymin, ymax, ny))
#    masks, area, midpoint  = compute_masks(xx, yy)
#
#    norm = 1
#    #norm = 1/(2.0*np.pi *sigma**2)
#    gg = np.exp( -(xx-mux)**2/(2*sigx**2) - (yy-muy)**2/(2*sigy**2) )
#
#    return(masks, gg, area)


def course_match(raw, sigma_x, mu_x, sigma):
    area = saq.area[:]
    nn = saq.N_SAQ_CHANNELS
    chi = []
    a_list = []
    for ii, diff in enumerate(sigma_x):
        for jj, offset in enumerate(mu_x):
#            irow = jj % len(mux)
#            icol = ii % len(sigx)
            template = integral_ind.loc[diff, offset]
            #A = sum(template[0][1:])/sum(raw[1:])
            A = sum((template[0][1:14]*raw[1:14])/(sigma[1:14]**2))/sum((raw[1:14]**2)/(sigma[1:14]**2))
            raw_norm  = raw[:]*A
            a_list.append(A)
            chi.append(chi_squared(raw_norm[1:9], template[0][1:9], sigma[1:9]))
            #chi.append(cost(raw[1:14], integral[num][1:14], sigma[1:14]))
            #print(chi_squared(raw[1:14], integral[num][1:14], sigma[1:14]))
    fit = np.array(chi).argmin()
#    plt.imshow(chi, cmap = 'gray')
#    plt.show()
#    print(chi)
    scaling  = a_list[fit]
    return(fit, integral.diffusion[fit], integral.offset[fit], scaling)
"""
example = integral.template[10]*50
sigma_ex = np.sqrt(integral.template[10])/area
print(sigma_ex)
#plt.plot(midpoint, integral.template[10])
noise = np.random.normal(0, .1, example.shape)
example = example + noise
parameter_ex = course_match(example, sigma_x, mu_x, sigma_ex)

print(parameter_ex)
plt.plot(midpoint, integral_ind.loc[parameter_ex[1], parameter_ex[2]][0], 'bo')
plt.plot(midpoint, example*parameter_ex[3], 'ko')
plt.show()
"""
for n in range(len(names)):

    parameters = course_match(data[:,n], sigma_x, mu_x, sigma[:,n])

    print(parameters)
    #sigma_2 = np.linspace(parameters[0] - 0.5, parameters[0] + 1, 4)
    #mu_2 = np.linspace(parameters[1] - 1, parameters[1] + 1, 4)
    #final_fit = course_match(data[:,n],sigma_2, mu_2, sigma[:,n])
    #print(final_fit)
    ll = integral_ind.loc[parameters[1], parameters[2]]
 

    #ll = np.zeros(16)
    #masks, gg, area  = make_gaussian(parameters[0], parameters[1])
    #for kk in range(16):
    #    ll[kk] = (np.sum(gg*masks[kk]))
    #
    #    #normalize based on area of ring
    #ll = ll/area
    #ll = ll/sum(ll)
    #print(ll)
    scale = 1/parameters[3]
    plt.plot(midpoint[1:], ll[0][1:]*scale, 'bo', label = 'model')
    #plt.plot(gg)
    y_error = sigma[1:,n]
    sigma[1:5,] = sigma[1:5,]*10
    #plt.plot(midpoint[1:], data[1:,n]/max(data[1:,n]), 'ko', label = 'data')
    plt.errorbar(midpoint[1:], data[1:,n], yerr = sigma[1:,n], fmt = 'ko', label = 'data')
    #plt.errorbar(midpoint[2:], data[2:,n],
#             yerr = y_error,
#             fmt ='o')
    plt.xlim(0, 20)
    plt.ylabel("Resets per mm^2")
    plt.xlabel("Radial distance (mm)")
    plt.legend()
    plt.show()


