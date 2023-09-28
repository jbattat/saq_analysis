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
    time_sec[i] = ts[i]/clock_rate + (n*loop_value) #the actual conversion 
    #add a count to the number of loops if the time has reset and fix that value
    if (int(ts[i]) - int(ts[i-1])) < 0 and i !=0:
         n+=1
         time_sec[i] = float(ts[i]/clock_rate) + float(n*loop_value) #the actual conversion
total_run = time_sec[-1]

print("Total runtime is: ", total_run, "seconds")
#print("Seconds/reset: ", total_run/len(data["Timestamp"])) 
###########################################################################
########## Calculate the distribution of resets and graph it ##############
###########################################################################

# create a list of the channels and all of their resets
chResets = [[t for t, mask in zip(time_sec, masks) if m(ch, mask)] for ch in range(saq.N_SAQ_CHANNELS)]
#chResets = [[t for t, mask in zip(time_sec, masks)] for ch in range(saq.N_SAQ_CHANNELS)]
channel_dist = np.zeros(saq.N_SAQ_CHANNELS)

#average the number of resets by the area of the ring
for ch, resets in enumerate(chResets):
    save_data[ch, 2] = len(resets)
    channel_dist[ch] = (len(resets)/area[ch])
    print("Resets per channel", str(ch + 1), ": ", len(resets))
save_data[:,1] = channel_dist


#print("# of resets per channel - averaged", repr(channel_dist))

#print(channel_dist)
#plt.plot(midpoint, channel_dist, 'o')
plt.xlim(0, 20)
plt.xlabel("Distance from center (mm)")
plt.ylabel("Number of resets per channel averaged over area")
plt.title("Distribution of resets")
outfile = input_file.replace('.root', '_dist.pdf')
plt.savefig(outfile)
#plt.show()
plt.clf()




###########################################################################
################### Calculate the average current  ########################
###########################################################################

#calculate the time between resets for each channel
channel1_resets = chResets[0]
rtd = [ [] for _ in range(16)]
current  = [ [] for _ in range(16)]

mask = [ [] for _ in range(16)]    
mean = np.zeros(saq.N_SAQ_CHANNELS)
std  = np.zeros(saq.N_SAQ_CHANNELS)

for ch in range(16):
    rtd[ch] = np.zeros(len(chResets[ch]))
    current[ch] = np.zeros(len(chResets[ch]))
    mask[ch] = np.zeros(len(chResets[ch]))

elements = len(max(chResets, key=len))
#tdr = np.zeros([16, elements])
#for ch in range(16):
#    tdr[:ch]  = np.linspace(0, total_run, elements)

voltage = np.array([.258, .238, .262, .266, .237, .249, .257, .236, .253, .240, .255, .262, .245, .258, .241, .247])*4
cap = 1e-11
#print(chResets[ch])
hold = 0
reset_fixed = np.zeros(16)
fig, axs = plt.subplots(2,5, figsize=(12,8))
    
for ch in range(10):
    remove_r = []
    for r in range(len(chResets[ch])-1):
        prev = chResets[ch][r+1]
        curr = chResets[ch][r]
        temp  = prev - curr
        rtd[ch][r] = temp
        if rtd[ch][r] < 1e-4:
            remove_r.append(r)
        else:
            if rtd[ch][r] != 0:
                current[ch][r] = voltage[ch]*cap/rtd[ch][r]/1e-9
            else:
                current[ch][r] = 0
    current[ch] = np.delete(current[ch], remove_r)
    chResets[ch] = np.delete(chResets[ch], remove_r)
    reset_fixed[ch] = len(chResets[ch])/area[ch]
    #plt.step(chResets[ch][1:-1], current[ch][1:-1], where ='post', label = str(ch+1))
    #plt.plot(chResets[ch][1:-1], current[ch][1:-1], 'o')
    row = int(ch/5)
    col = ch%5
    lower_bound = 1e-4
    upper_bound = 1e-2
    mask = (rtd[ch] >= lower_bound) & (rtd[ch] <= upper_bound)
    mean[ch] = np.mean(rtd[ch][mask])

    std[ch] = np.std(rtd[ch][mask])


    #mean[ch] = np.mean(rtd[ch])

    #std[ch] = np.std(rtd[ch])
    print("ch", ch+1, " ,mean: ",  mean[ch], " ,std: ", std[ch])
    axs[row,col].hist(rtd[ch][mask], bins = 100)
    #axs[row,col].hist(rtd[ch], bins = 100)
    axs[row,col].axvline(mean[ch], color = 'red')
    axs[row,col].axvline(mean[ch]-std[ch] , color = 'green')
    axs[row,col].axvline(mean[ch]+std[ch] , color = 'green')
    

#    print(rtd[ch])

#print(max(rtd[0]))
#plt.hist(rtd[0], 100)
#plt.xlim(0, 0.5)
#print("average rtd: ", np.median(rtd[0]))
plt.show()
sys.exit()
#plt.step(chResets[1][:-1], current[1][:-1], where ='post', label = "Channel 1")
plt.title("Reconstructed Current")
plt.xlabel("Time (s)")
plt.ylabel("Current (nA)")
plt.legend()

outfile = input_file.replace('.root', '_current.pdf')
plt.savefig(outfile)
plt.show()
plt.clf()
plt.plot(midpoint, reset_fixed, 'o')
plt.show()
plt.hist(rtd[7], bins = 50)
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








