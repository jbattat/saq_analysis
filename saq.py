import numpy as np
import ROOT
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker

N_SAQ_CHANNELS = 16
CLOCK_FREQ = 30.3e6 # Hz (frequency of the zybo board clock)
# fixme: CLOCK_FREQ comes from data file directly!
N_CLOCK_BITS = 32 # used to determine clock wraps


# Measured using GerbView by James Battat on 2023 June 8
# for each channel, rmin is the inner radius and rmax is the outer radius
# inner and outer radii (to midpoints of gaps between channels) in millimeters
rmin = [ 0.000,  0.670,  1.386,  2.106,  2.981,  4.005,  4.994,  6.015, 
         7.481,  9.994, 12.497, 15.023, 19.996, 24.962, 30.026, 39.977 ]

rmax = [ 0.670,  1.386,  2.106,  2.981,  4.005,  4.994,  6.015,  7.481,
         9.994, 12.497, 15.023, 19.996, 24.962, 30.026, 39.977, 50.065]
rmin = np.array(rmin)
rmax = np.array(rmax)

# for backwards compatibility
annuli = [ [rmin[ii], rmax[ii]] for ii in range(len(rmin)) ]

# These values are old -- do not use...
#rmin = [ 0.000, 0.640,  1.478,  2.205,  3.095,  4.105,  5.100,  6.115,
#         7.595, 9.089, 11.590, 16.505, 21.460, 26.460, 31.495, 41.440]
#rmax = [ 0.640, 1.4775,  2.205,  3.095,  4.105,  5.100,  6.115,  7.595,
#         9.085, 11.590, 16.505, 21.460, 26.460, 31.495, 41.440, 51.275]

area = np.pi*(rmax**2 - rmin**2)

#area_OLD =  ([1.28679635e+00, 5.57132005e+00, 8.41638562e+00, 1.48188925e+01,
#              2.28456618e+01, 2.87737686e+01, 3.57614560e+01, 6.37454282e+01,
#              7.80786305e+01, 1.62677290e+02, 4.33812869e+02, 5.90985650e+02,
#              7.52725600e+02, 9.16727496e+02, 2.27871834e+03, 2.86466762e+03])

radius_of_channel = np.sqrt(0.5*(rmax**2 + rmin**2))

# DEFUNCT: I don't know what "midpoint" is meant to be.
midpoint = ([ 0.255,  1.090,  1.855,  2.645,   3.605,   4.605,  5.605,
              6.855,  8.335, 10.335, 14.0475, 18.9625, 23.960, 28.960,
              36.44, 46.45  ])

# File utility functions
def get_metadata(filename, verbose=False):
    """ returns a dictionary with the meta data from a SAQ ROOT file """
    metadata = ROOT.RDataFrame("mt", filename).AsNumpy()
    if verbose:
        print(metadata)
    return {'VERSION':metadata['Version'],
            'ZYBO_FRQ':metadata['Zybo_FRQ'],
            'FILE_TIMESTAMP':metadata['Date'],
            'SAQ_DIV':metadata['SAQ_DIV']
            }

def get_raw_data(filename):
    # open up ttree into an rdataframe --> dictionary
    # dict keys are the branch names of the TTree
    return ROOT.RDataFrame("tt", filename).AsNumpy()

def resets_by_channel(data):
    """ returns the timestamp of each reset as a list of lists: resets[chan][reset_i] """

    # create a list of the channels and all of their resets
    return [[t for t, mask in zip(data["Timestamp"], data["ChMask"]) if m(ch, mask)]
            for ch in range(N_SAQ_CHANNELS)]

def unwrap_resets(resets_wrapped, counter_bits):
    """ unwrap reset timestamps
    
    inputs:
    .   resets: counter values corresponding to recorded resets
    .   counter_bits: number of bits for the counter (e.g. 32)
    
    output:
    .   the unwrapped (monotonically increasing) counter value
    """
    return np.unwrap(resets_wrapped, period=0.5*2**counter_bits)


def rtds_of_resets(resets, SAQ_DIV, ZYBO_FRQ):
    """ given a list of reset timestamps, compute the reset time differences (RTDs)
    
    Assumes that the reset timestamps have been unwrapped, but are in clock counts
    (not time in seconds)

    inputs:
    .    resets: list of resets timestampe (clock counts) ***unwrapped***
    .    SAQ_DIV: divisor factor applied to the raw FPGA clock
    .    ZYBO_FRQ: frequency of FPGA clock

    e.g. if the zybo clock is 20MHz (50ns period) and the
    divisor is 1,000 then a timestamp of 500 corresponds to
    500*1000 ticks of the 20MHz clock, or 500*100*50ns = 25ms
    """
    rsts = np.array(resets)
    rtds = rsts[1:]-rsts[0:-1]
    rtds *= SAQ_DIV / ZYBO_FRQ
    return rtds

def channel_made_reset(chan, mask):
    """ check if a specific channel is included in the channel mask """
    # Returns True if channel "chan" is present in the mask "mask"
    return ((1 << chan) & mask) > 0

def channels_in_reset_mask(mask):
    chlist =  [n for n in range(16) if ((mask >> n) & 1)]
    return chlist

#####################################
# basic analyis routines

def study_resets(data):
    # in progress... want to make sure there are no
    # outlier channels (ones that reset too often, or are high when other chans are high)
    # don't expect many resets in which multiple channels reset at the same time...

    # reset masks are saved in data['ChMask']
    n_resets = len(data['ChMask'])

    print(data['ChMask'])
    
    zeros   = []
    singles = []
    chs = []
    multiplicity = np.zeros(N_SAQ_CHANNELS+1, dtype=np.int32)  # tally # of hits from 0 to 16
    
    hots  = np.zeros(16, dtype=np.int32) # "hot" pixels...? (# of times pixel N fired when 2+ pixels fired)
    hots2 = np.zeros(16, dtype=np.int32) # "hot" pixels...? (# of times pixel N fired when exactly 2 pixels fired)
    hots3 = np.zeros(16, dtype=np.int32) # "hot" pixels...? (# of times pixel N fired when exactly 3 pixels fired)
    hots4p = np.zeros(16, dtype=np.int32) # "hot" pixels...? (# of times pixel N fired when 4 or more pixels fired)

    freq = np.zeros(16, dtype=np.int32) # how many times did each pixel fire?
    
    for ii, mask in enumerate(data['ChMask']):
        chs.append(channels_in_reset_mask(mask))
        nhits = len(chs[-1])
        multiplicity[nhits] += 1

        for ch in chs[-1]:
            freq[ch] += 1
        
        # which channels are present in multi-hit events?
        if nhits > 1:
            for ch in chs[-1]:
                hots[ch] += 1
                
            if nhits == 2:
                for ch in chs[-1]:
                    hots2[ch] += 1
            elif nhits == 3:
                for ch in chs[-1]:
                    hots3[ch] += 1
            elif nhits > 3:
                for ch in chs[-1]:
                    hots4p[ch] += 1
                
    print('multiplicity:')
    print(multiplicity)
    print()
    print('hots:')
    print(hots)
    print('hots2:')
    print(hots2)
    print('hots3:')
    print(hots3)
    print('freq:')
    print(freq)
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(15,8))
    axs[0,0].bar(np.arange(17, dtype=np.int32), multiplicity)
    #axs[0,0].xaxis.set_major_locator(MaxNLocator(integer=True))
    loc = plticker.MultipleLocator(base=1.0) 
    axs[0,0].xaxis.set_major_locator(loc)
    axs[0,0].set_xlim(-0.5,16.5)
    axs[0,0].set_xlabel("Number of channels that fired")
    axs[0,0].set_ylabel("Frequency")
    axs[0,0].set_title(f"Channel Multiplicity (Total # of resets: {n_resets})")

    axs[0,1].bar(np.arange(16, dtype=np.int32), freq)
    axs[0,1].xaxis.set_major_locator(loc)
    axs[0,1].set_xlim(-0.5,15.5)
    axs[0,1].set_xlabel("Channel Number")
    axs[0,1].set_ylabel("Number of times each channel fired")
    #axs[0,1].set_title(f"Channel Multiplicity (Total # of resets: {n_resets})")

    axs[1,1].bar(np.arange(16, dtype=np.int32), hots4p)
    axs[1,1].xaxis.set_major_locator(loc)
    axs[1,1].set_xlim(-0.5,15.5)
    axs[1,1].set_xlabel("Channel Number")
    axs[1,1].set_ylabel("Number of times each channel fired")
    axs[1,1].set_title(f"Hot channels(?): (requires 4+ channels to have fired)")

    plt.savefig('junk.pdf')

    
def files_from_list(flist):
    with open(flist) as fh:
        files = fh.readlines()
    files = [x.strip() for x in files]
    #print(files)
    return files

def convert_time(time):
    clock_rate = 30.3e6 #the frequency of the zybo board
    n = 0 #number of loops through the clock

    time_sec = np.zeros(len(time))
    #Convert the individual entries to show the time in seconds
    for i in range(len(time)):
        cutoff = (2**32)-1 #the value at which the zybo's clock resets back to zero
        loop_value = (cutoff)/clock_rate #the number of seconds it takes to go through one loop
        time_sec[i] = time[i]/clock_rate + (n*loop_value) #the actual conversion 
        #add a count to the number of loops if the time has reset and fix that value
        if (int(time[i]) - int(time[i-1])) < 0 and i !=0:
            n+=1
            time_sec[i] = float(time[i]/clock_rate) + float(n*loop_value) #the actual conversio
    return time_sec

def calculate_rtd(resets):
    #calculate the time between resets for each channel
    rtd = np.zeros(len(resets)) 
    hold = 0

    for r in range(len(resets)-1):
        prev   = resets[r+1]
        curr   = resets[r]
        temp   = prev - curr
        rtd[r] = temp
    return rtd
##############################################################
# Compute the area of each annular ring on the anode
# The code below was used to generate the hard-coded "area" numbers above

##width is the dimension of the electrode and spacing is the distance between them 
#width = np.array([0.51, 0.64, 0.64, 0.82, 0.94, 0.94, 0.96, 1.42, 1.42, 2.42, 4.825, 4.825, 4.83, 4.83, 9.65, 9.65])
#spacing = np.array([0, 0.26, 0.125, 0.06, 0.08, 0.06, 0.05, 0.06, 0.06, 0.08, 0.09, 0.09, 0.17, 0.17, 0.24, 0.36, 0])
#inner_radii = np.zeros(N_SAQ_CHANNELS)
#outer_radii = np.zeros(N_SAQ_CHANNELS)
#area = np.zeros(N_SAQ_CHANNELS)
#middle_gap = np.zeros(N_SAQ_CHANNELS)
#sum_r = 0
#inner = 0
#
##turn the widths and spacings into an inner and outer radii
#for x in range(len(width)):
#    inner_radii[x] = sum_r + spacing[x]
#    outer_radii[x] = inner_radii[x] + width[x]
#    middle_gap[x] = outer_radii[x] + spacing[x+1]/2
#    sum_r = sum_r + width[x] + spacing[x]
#    
#
##calculate the mid_point of each ring
#mid_point =  (inner_radii + outer_radii)/2
#print(repr(mid_point))
##calculate the area of each ring
#for x in range(N_SAQ_CHANNELS):
#    if x != 0:
#        area[x] = (middle_gap[x]**2 - middle_gap[x-1]**2)*np.pi
#    else:
#        area[x] = np.pi*middle_gap[x]**2
##area = (outer_radii**2 - inner_radii**2)*np.pi
##print(area)

