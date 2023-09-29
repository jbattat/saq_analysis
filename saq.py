import numpy as np

N_SAQ_CHANNELS = 16
CLOCK_FREQ = 30.3e6 # Hz (frequency of the zybo board clock)

area =  ([1.28679635e+00, 5.57132005e+00, 8.41638562e+00, 1.48188925e+01,
       2.28456618e+01, 2.87737686e+01, 3.57614560e+01, 6.37454282e+01,
       7.80786305e+01, 1.62677290e+02, 4.33812869e+02, 5.90985650e+02,
       7.52725600e+02, 9.16727496e+02, 2.27871834e+03, 2.86466762e+03])

midpoint = ([ 0.255 ,  1.09  ,  1.855 ,  2.645 ,  3.605 ,  4.605 ,  5.605 ,
        6.855 ,  8.335 , 10.335 , 14.0475, 18.9625, 23.96  , 28.96  ,
       36.44  , 46.45  ])



# utility functions

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

