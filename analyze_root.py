# robust fileroot determination
# run from any directory, save outputs to dir that has .root file

import ROOT
import sys
import csv

import math
import numpy as np
import matplotlib.pyplot as plt
import struct
from array import array
from ROOT import TCanvas

import saq
from scipy.optimize import curve_fit

###########################################################################
#################### calculate the widths of the anodes ###################
###########################################################################
area = saq.area[:]
midpoint = saq.midpoint[:]

save_data = np.zeros((saq.N_SAQ_CHANNELS, 3))
save_data[:,0] = np.arange(saq.N_SAQ_CHANNELS)+1  # channel numbers 1-16

###########################################################################
################# import the data from the ROOT file ######################
###########################################################################

#import the ROOT file that we want to analyze
input_file = sys.argv[1]

#import the tree from the ROOT file and break it into arrays of the times and the channels
rdf = ROOT.RDataFrame("tt", input_file)
data = rdf.AsNumpy()
#print("Number of Resets: ", len(data["Timestamp"]))

#print(data)
#print(type(data))

ts = data["Timestamp"]
masks = data["ChMask"]

# make a quick way to ensure the channel we want is in the mask
m = lambda ch, mask: 1 << ch & mask

###########################################################################
############# Conver the time of the resets into seconds ##################
###########################################################################

clock_rate = 30.3e6 #the frequency of the zybo board
n = 0 #number of loops through the clock

time_sec = np.zeros(len(ts))
#Convert the individual entries to show the time in seconds
for i in range(len(ts)):
    cutoff = (2**32)-1 #the value at which the zybo's clock resets back to zero
    loop_value = (cutoff)/clock_rate #the number of seconds it takes to go through one loop

    ts[i] = ts[i]/clock_rate + (n*loop_value) #the actual conversion
    #add a count to the number of loops if the time has reset and fix that value
    if (int(ts[i]) - int(ts[i-1])) < 0 and i !=0:
         n+=1
         ts[i] = ts[i]/clock_rate + n*loop_value #the actual conversion

total_run = ts[-1]
#print("Total runtime is: ", total_run, "seconds")
#print("Seconds/reset: ", total_run/len(data["Timestamp"])) 
###########################################################################
########## Calculate the distribution of resets and graph it ##############
###########################################################################

# create a list of the channels and all of their resets
chResets = [[t for t, mask in zip(ts, masks) if m(ch, mask)] for ch in range(saq.N_SAQ_CHANNELS)]
channel_dist = np.zeros(saq.N_SAQ_CHANNELS)

#average the number of resets by the area of the ring
for ch, resets in enumerate(chResets):
    save_data[ch, 2] = len(resets)
    channel_dist[ch] = (len(resets)/area[ch])
    #print("Resets per channel", str(ch + 1), ": ", len(resets))
save_data[:,1] = channel_dist
#print(save_data)
#channel_dist = channel_dist/max(channel_dist)
dark = [1.,0.23290509, 0.15205463, 0.17366971, 0.05632563, 0.04472116,
        0.03598277, 0.02018649, 0.03296155, 0.01581746, 0.00741562, 0.00979818,
        0.01196661, 0.01263316, 0.,         0.        ]
#channel_dist = channel_dist - dark 
#print("# of resets per channel - averaged", repr(channel_dist))

#print(channel_dist)
plt.plot(midpoint, channel_dist, 'o')
plt.xlim(0, 20)
plt.xlabel("Distance from center (mm)")
plt.ylabel("Number of resets per channel averaged over area")
plt.title("Distribution of resets")
outfile = input_file.replace('.root', '_dist.pdf')
plt.savefig(outfile)
plt.clf()
#plt.show()

###########################################################################
################### Fit the Distribution of Resets  #######################
###########################################################################
"""
def func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

#mean = sum(midpoint[4:]*channel_dist[4:])/sum(channel_dist[4:])
mean = sum(midpoint*channel_dist)/sum(channel_dist)                  
#sigma = sum(channel_dist[4:]*(midpoint[4:]-mean**2)/sum(channel_dist[4:]))
#print(sigma)
sigma = sum(channel_dist*(midpoint-mean)**2)/sum(channel_dist)
popt, pcov = curve_fit(func, midpoint, channel_dist, p0=[max(channel_dist),mean,sigma])
#popt, pcov = curve_fit(func, midpoint[4:], channel_dist[4:], p0=[max(channel_dist),mean,1.5])
print(popt)
plt.plot(midpoint, channel_dist, 'o', label = 'data')
xm = np.linspace(0, 15, 100)
ym = func(xm, popt[0], popt[1], popt[2])
plt.plot(xm, ym, c='r', label = 'model')
plt.xlabel("Distance from center (mm)")
plt.ylabel("Number of resets averaged over annulus area")
plt.xlim(0, 15)
plt.legend()
outfile = input_file.replace('.root', '_model.pdf')
plt.savefig(outfile)
#plt.show()
"""
###########################################################################
################### Calculate the average current  ########################
###########################################################################

#calculate the time between resets for each channel
channel1_resets = chResets[0]
rtd = [ [] for _ in range(16)]
current  = [ [] for _ in range(16)]
for ch in range(16):
    rtd[ch] = np.zeros(len(chResets[ch]))
    current[ch] = np.zeros(len(chResets[ch]))

elements = len(max(chResets, key=len))
#tdr = np.zeros([16, elements])
#for ch in range(16):
#    tdr[:ch]  = np.linspace(0, total_run, elements)

voltage = np.array([258, 238, 262, 266, 237, 249, 257, 236, 253, 240, 255, 262, 245, 258, 241, 247])*4
cap = 1e-11
#print(chResets[ch])
hold = 0
for ch in range(13):
    for r in range(len(chResets[ch])-1):
        prev = chResets[ch][r+1]
        curr = chResets[ch][r]
        temp  = prev - curr
        rtd[ch][r] = temp
        if rtd[ch][r] != 0:
            current[ch][r] = voltage[ch]*cap/rtd[ch][r]/1e-9
        else:
            current[ch][r] = 0
    plt.step(chResets[ch][:-1], current[ch][:-1], where ='post', label = str(ch+1))

#plt.step(chResets[0][:-1], current[0][:-1], where ='post', label = "Channel 1")
plt.title("Reconstructed Current")
plt.xlabel("Time (s)")
plt.ylabel("Current (nA)")
plt.legend()
outfile = input_file.replace('.root', '_current.pdf')
plt.savefig(outfile)
#plt.show()

###########################################################################
############################ Save Data  ###################################
###########################################################################

header = ["SAQ_Channel_Numbers", "Resets_Averaged_by_Area", "Raw_Number_of_Resets"]

outfile = input_file.replace('.root', '_data.csv')
with open(outfile, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # write multiple rows
    writer.writerows(save_data)








