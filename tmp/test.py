'''
---------------------------
    pisciarelli.py
---------------------------
Created on 08.02.2016
Last modified on 09.02.2016
Author: Marc Wieland
Description: Single-well pump test to derive transmissivity from waterlevel measurements.
             0. Resample and/or slice waterlevel and airpressure timeseries (optional).
             1. Correct waterlevel with airpressure data
             2. Auto-detect drawdown events in waterlevel timeseries.
             3. Compute Transmissivity (T) as follows: 
                 T = 2.3 * Q / 4 * Pi * dsw
                 Q = pumping rate
                 dsw: change in drawdown per log-cycle      
Reference: Cooper and Jacob (1946): A generalized method for evaluating formation constants
           and summarizing well-field history. Transactions of the American Geophysical
           Union, 27 (5), 526-534.

TODO: x 1. test plot raw wL data (and also corrected wL data) for April 10,2014 -> see if there are bumps
      2. interpolate timeseries data (e.g., spline) to 5min (1min)
      3. extract all wL values into csv for A) pump events (between min->max) and 
                                            B) drawdown events (between max->min)
----
'''

import time
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import timeseries.helper as helper
import matplotlib.dates as mdates

# Parameters to set ########################################################################################################
dat_wl = 'PITC (DCX22).csv' # water level data [mbar]
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_pitc/marc/'
###
t_slice = True    # temporal slicing
slice_t1 = '2014-07-01 07:00'   # slicing start time
slice_t2 = '2014-07-01 15:00'   # slicing end time
###
t_int = True # temporal interpolation 
ir = "1min" # temporal interpolation rate (e.g., 5min, H, D, M)
if_d = 'cubic' # temporal interpolation function for drawdown events 
if_f = 'cubic'  # temporal interpolation function for fill events                
                # {'linear', 'time', 'index', 'values', 'nearest', 'zero', 'slinear', , 'quadratic',
                #  'cubic', 'barycentric', 'krogh', 'polynomial', 'spline', 'piecewise_polynomial', 'pchip'}
###
peak_lookahead = 1 #200     # distance to look ahead from a peak candidate to determine if it is the actual peak
############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

####################################
### Data input and preprocessing ###
####################################
# Read data from csv to pandas dataframe
df1 = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 

# Convert waterlevel (wL) from [mbar] to [m]
df1['wL#'] = df1['wL#'] * 0.0101972

# Create waterlevel timeseries
dates1 = pd.DatetimeIndex(df1['longDATE'].values)
wL = pd.Series(df1['wL#'].values, index=dates1)

if t_slice is True:
    # Slice timeseries
    wL = wL[slice_t1:slice_t2]

##########################################################
### Compute transmissivity (T) (Cooper and Jacob 1946) ###
##########################################################
# TODO: this part is not yet correct i think
fig = plt.figure()
# Plot raw waterlevel timeseries
plt.plot(wL.index, wL.values, 'ro', mew=2, ms=8, label='raw (5min)')

#######################
### Event detection ###
#######################
# Detect drawdown events in waterlevel timeseries and return h1, h2 and dt
events = helper.drawdown_event_detection(wL, peak_lookahead)
print events

#####################################################################
### Separate waterlevels for each event between drawdown and fill ###
#####################################################################
fig = plt.figure()
# Plot raw data wL values
plt.plot(wL.index, wL.values, 'wo', mew=2, ms=8, label='wL_corr (10min)')

wL_d = None
wL_f = None

for e in range(events.shape[0]):
    if e == 0:
        # Get waterlevels within first drawdown event
        wL_d = wL[events['h1_t'][e]:events['h2_t'][e]]
        if t_int is True:
            wL_d = wL_d.resample(ir, how='mean')    
            wL_d = wL_d.interpolate(method=if_d)
        wL_d.plot() 
        wL_d = wL_d.to_frame()
        wL_d['event'] = e 
        if e != events.shape[0]-1:
            # Get waterlevels within first fill event
            wL_f = wL[events['h2_t'][e]:events['h1_t'][e+1]]
            if t_int is True:
                wL_f = wL_f.resample(ir, how='mean')
                wL_f = wL_f.interpolate(method=if_f)
            wL_f.plot()
            wL_f = wL_f.to_frame()
            wL_f['event'] = e
    else:
        # Get waterlevels within all other drawdown events
        wL_d_ = wL[events['h1_t'][e]:events['h2_t'][e]]
        if t_int is True:
            wL_d_ = wL_d_.resample(ir, how='mean')    
            wL_d_ = wL_d_.interpolate(method=if_d)       
        wL_d_.plot()
        wL_d_ = wL_d_.to_frame()
        wL_d_['event'] = e 
        wL_d = wL_d.append(wL_d_)
        if e != events.shape[0]-1:
            # Get waterlevels within all other fill events
            wL_f_ = wL[events['h2_t'][e]:events['h1_t'][e+1]]
            if t_int is True:
                wL_f_ = wL_f_.resample(ir, how='mean')
                wL_f_ = wL_f_.interpolate(method=if_f)
            wL_f_.plot()
            wL_f_ = wL_f_.to_frame()
            wL_f_['event'] = e 
            wL_f = wL_f.append(wL_f_)

# Plot event markers
plt.plot(events['h1_t'], events['h1'], 'r+', mew=2, ms=8)
plt.plot(events['h2_t'], events['h2'], 'g+', mew=2, ms=8)
plt.xlabel('Time')
plt.ylabel('wL [m]')
plt.grid()
plt.tight_layout()
plt.savefig(wdir + 'int_cubic_separated.png', dpi=300)
plt.show()
plt.close()

# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'