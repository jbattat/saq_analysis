# robust fileroot determination
# run from any directory, save outputs to dir that has .root file

import sys
import numpy as np
import matplotlib.pyplot as plt
import ROOT
sys.path.append("../")
import saq

###########################################################################
#################### calculate the widths of the anodes ###################
###########################################################################
area = saq.area[:]
midpoint = saq.midpoint[:]

###########################################################################
################# import the data from the ROOT file ######################
###########################################################################

#import the ROOT file that we want to analyze
input_file = sys.argv[1]

#import the tree from the ROOT file and break it into arrays of the times and the channels
rdf = ROOT.RDataFrame("tt", input_file)
data = rdf.AsNumpy()

ts = data["Timestamp"]
masks = data["ChMask"]

# make a quick way to ensure the channel we want is in the mask
m = lambda ch, mask: 1 << ch & mask

###########################################################################
############# Convert the time of the resets into seconds #################
###########################################################################
# FIXME: move this to saq.py

clock_rate = saq.CLOCK_FREQ #the frequency of the zybo board
nn = 0 # number of loops through the clock counter (did it wrap?)

cutoff = (2**32)-1 # Counter is 32 bits (highest counter value)
loop_time = cutoff/clock_rate # (seconds) Time for counter to wrap

time_sec = np.zeros(len(ts))
#Convert the individual entries to show the time in seconds
for i in range(len(ts)):
    time_sec[i] = ts[i]/clock_rate + (nn*loop_time) #the actual conversion 
    #add a count to the number of loops if the time has reset and fix that value
    if (int(ts[i]) - int(ts[i-1])) < 0 and i !=0:
         nn+=1
         time_sec[i] = float(ts[i]/clock_rate) + float(nn*loop_time) #the actual conversion
total_run = time_sec[-1]

print("Total runtime is: ", total_run, "seconds")

###########################################################################
########## Calculate the distribution of RTDs and graph it   ##############
###########################################################################

# create a list of the channels and all of their resets
chResets = [[t for t, mask in zip(time_sec, masks) if m(ch, mask)] for ch in range(saq.N_SAQ_CHANNELS)]

#calculate the time between resets for each channel
rtd = [ [] for _ in range(saq.N_SAQ_CHANNELS)]
mask = [ [] for _ in range(saq.N_SAQ_CHANNELS)]  # RTDs to ignore
mean = np.zeros(saq.N_SAQ_CHANNELS) 
std  = np.zeros(saq.N_SAQ_CHANNELS)

for ch in range(saq.N_SAQ_CHANNELS):
    rtd[ch] = np.zeros(len(chResets[ch]))
    mask[ch] = np.zeros(len(chResets[ch]))

fig, axs = plt.subplots(2,5, figsize=(12,8))
    
for ch in range(10): # only calibrated first 10 channels...
    for r in range(len(chResets[ch])-1):
        prev = chResets[ch][r+1]
        curr = chResets[ch][r]
        rtd[ch][r] = prev - curr
    lower_bound = 1e-4
    upper_bound = 1e-2
    mask = (rtd[ch] >= lower_bound) & (rtd[ch] <= upper_bound)
    mean[ch] = np.mean(rtd[ch][mask])

    std[ch] = np.std(rtd[ch][mask])
    #std[ch] = np.std(rtd[ch])

    print("ch", ch+1, " ,mean: ",  mean[ch], " ,std: ", std[ch])

    row = int(ch/5)
    col = ch%5
    axs[row,col].hist(rtd[ch][mask], bins = 100)
    #axs[row,col].hist(rtd[ch], bins = 100)
    axs[row,col].axvline(mean[ch], color = 'red')
    axs[row,col].axvline(mean[ch]-std[ch] , color = 'green')
    axs[row,col].axvline(mean[ch]+std[ch] , color = 'green')
    
outfile = input_file.replace('.root', '_hists.pdf')
plt.savefig(outfile)
