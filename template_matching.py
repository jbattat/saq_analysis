# 16 --> saq.N_SAQ_CHANNELS

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import sys
import saq

from offset_simulation import compute_masks


integral  = pd.DataFrame(pickle.load(open('./template_50.pkl', 'rb')))
sigma_x   = integral.diffusion.unique()
mu_x      = integral.offset.unique()
integral_ind  = integral.set_index(['diffusion', 'offset'])

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
    scale = str(1 -  (x+1)/15)
    #plt.plot(midpoint[1:], data[1:,x], label = ('Run ' + str(x+1)), color = scale)
    plt.plot(midpoint[1:], data[1:,x], label = ('Pressure ' + str(200 +50*x)), color = scale)
plt.xlim(0,20)
plt.xlabel("Radial Distance (mm)")
plt.ylabel("Number of Resets averaged by area (mm^2)")
plt.title("Diffusion Scan over Pressure Range of 200 torr to 750 torr")
plt.legend()
plt.savefig("Drift_Field_400.pdf")
plt.clf()

def chi_squared(data, model, sigma):

    chi_squared = np.sum(((data - model)/sigma)**2)
    #chi_squared = np.sum(((data - model)**2))
    
    return chi_squared

def fixed_offset(raw, mu, sigma):
    nn = saq.N_SAQ_CHANNELS
    chi = []
    a_list = []
    for ii, diff in enumerate(sigma_x):
        template = integral_ind.loc[diff, mu]
        start = 4
        A = sum((template[0][start:10]*raw[start:10])/(sigma[start:10]**2))/sum((raw[start:10]**2)/(sigma[start:10]**2))
        template = template/A
        a_list.append(A)
        chi.append((chi_squared(raw[start:9], template[0][start:9], sigma[start:9])))
    fit = np.array(chi).argmin()
    scaling = a_list[fit]
    return(fit, sigma_x[fit], scaling, chi[fit])
        
def coarse_match(raw, sigma_x, mu_x, sigma):
    area = saq.area[:]
    nn = saq.N_SAQ_CHANNELS
    chi = []
    a_list = []
    for ii, diff in enumerate(sigma_x):
        for jj, offset in enumerate(mu_x):
#            irow = jj % len(mux)
#            icol = ii % len(sigx)
            template = integral_ind.loc[diff, offset]
            #A = sum(template[0][:])/sum(raw[:])
            #start = next((i for i, x in enumerate(sigma) if x), None)
            start = 1
            A = sum((template[0][start:10]*raw[start:10])/(sigma[start:10]**2))/sum((raw[start:10]**2)/(sigma[start:10]**2))
            template = template/A            

            a_list.append(A)
            chi.append((chi_squared(raw[start:9], template[0][start:9], sigma[start:9])))
        
    fit = np.array(chi).argmin()
#    plt.imshow(chi, cmap = 'gray')
#    plt.show()
#    print(chi)
    scaling  = a_list[fit]
    return(fit, integral.diffusion[fit], integral.offset[fit], scaling)

"""
example  = []
sigma_ex = []

for n in range(10):
    example.append(integral.template[n*3]*50)
#    plt.plot(midpoint, integral.template[10])
    example_n = example[n][:]
    noise = np.random.normal(0, 5, example_n.shape)
    example_n = example_n + noise
    parameter_ex = course_match(example_n, sigma_x, mu_x, sigma_ex)

    print(parameter_ex)
    plt.plot(midpoint, integral_ind.loc[parameter_ex[1], parameter_ex[2]][0], label = "model")
    plt.plot(midpoint, example_n*parameter_ex[3], 'ko', label = "data")
    plt.plot(midpoint, example[n]*parameter_ex[3], 'bo', label = "clean")
    plt.legend()
    plt.show()

"""

axis = int(np.ceil(np.sqrt(number)))
fig, axs = plt.subplots(3, 4, layout = 'constrained')
chis = []

#calculate the best offset for this set of data
for ii, offset in enumerate(mu_x):
    chi_total = 0
    for n in range(len(names)):
        chi_min = 0
        chi_min = fixed_offset(data[:,n], offset, sigma[:,n])
        chi_total = chi_total + chi_min[3]
    chis.append(chi_total)
mu = np.array(chis).argmin()
print("best offset:", mu_x[mu])

#match the width for each set of data
for n in range(len(names)):
    parameters = fixed_offset(data[:,n], mu_x[mu], sigma[:,n])
    ll = integral_ind.loc[parameters[1], mu_x[mu]]
 
    icol = int(n % axis)
    irow = int(n/  axis)
    scale = 1/parameters[2]
    y_error = sigma[1:,n]

    axs[irow, icol].plot(midpoint[1:], ll[0][1:]*scale, 'b', label = 'model')
    axs[irow, icol].errorbar(midpoint[1:], data[1:,n], yerr = sigma[1:,n], fmt = 'ko', label = 'data')
    axs[irow, icol].set_title(f'{n + 1}, Width: {round(parameters[1], 2)}, Offset: {round(mu_x[mu], 2)}')
    #axs[irow, icol].set_title(f'Run {n+1}')

    axs[irow, icol].set_xlim(0, 20)

#    plt.errorbar(midpoint[1:], data[1:,n], yerr = sigma[1:,n], fmt = 'ko', label = 'data')
    #plt.errorbar(midpoint[2:], data[2:,n],
#             yerr = y_error,
#             fmt ='o')
    plt.xlim(0, 20)
    plt.ylabel("Resets per mm^2")
    plt.xlabel("Radial distance (mm)")
#    plt.legend()
#    plt.show()

plt.show()


"""
parameters = coarse_match(data[:,5], sigma_x, mu_x, sigma[:,5])
scale      = 1/parameters[3]
ll = []
ll.append(integral.iloc[parameters[0]- 200][2])
ll.append(integral.iloc[parameters[0] - 100][2])
ll.append(integral.iloc[parameters[0]][2])
ll.append(integral.iloc[parameters[0] + 100][2])
ll.append(integral.iloc[parameters[0] + 200][2])
                                                  
plt.plot(midpoint[1:], data[1:, 5], label = "data")
for n in ll:
    print(n)
    plt.plot(midpoint[1:], n[1:], label = f"fit + {-200 + 100*n}")
plt.legend()
plt.show()

"""
"""
#Stuff to plot the templates
plot_sigma = [0, 10, 20, 30, 40]
plot_mu = [0, 12, 24, 35]

fig, axs = plt.subplots(len(plot_sigma), len(plot_mu), layout = 'constrained')

for ii, sig in enumerate(plot_sigma):
    for jj, mu in enumerate(plot_mu):
        ss = np.array(integral_ind.loc[sigma_x[sig], mu_x[mu]])
        print(len(ss[0]))
        irow = ii % len(plot_sigma)
        icol = jj % len(plot_mu)

        axs[irow, icol].plot(midpoint, ss[0], 'bo')
        axs[irow, icol].set_title(f'Width: {round(sigma_x[sig], 2)}, Offset: {round(mu_x[mu], 2)}')
        #axs[irow, icol].set_title(f'Run {n+1}')

        axs[irow, icol].set_xlim(0, 20)

#    plt.errorbar(midpoint[1:], data[1:,n], yerr = sigma[1:,n], fmt = 'ko', label = 'data')
    #plt.errorbar(midpoint[2:], data[2:,n],
#             yerr = y_error,
#             fmt ='o')
    plt.xlim(0, 20)
    plt.ylabel("Resets per mm^2")
    plt.xlabel("Radial distance (mm)")
#    plt.legend()
#    plt.show()

plt.show()
"""
