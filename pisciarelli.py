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
import pandas as pd
import matplotlib.pyplot as plt
import timeseries.helper as helper

# Parameters to set ########################################################################################################
dat_wl = 'PITC (DCX22).csv'
dat_bar = 'AGN3baro.csv'
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_pitc/marc/'
###
t_res = True # temporal resampling 
sr = "15min" # temporal resampling rate (e.g., 5min, H, D, M)
stat = 'mean' # temporal resampling statistics 
###
t_slice = True    # temporal slicing
slice_t1 = '2014-09-20 00:01'   # slicing start time
slice_t2 = '2014-09-20 23:59'   # slicing end time
###
peak_lookahead = 1 #200     # distance to look ahead from a peak candidate to determine if it is the actual peak
###
Q = 18  # Pumping rate [cubic meters per hour]
############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

####################################
### Data input and preprocessing ###
####################################
# Read data from csv to pandas dataframe
df1 = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 
df2 = pd.read_table(wdir + dat_bar, header=0, delimiter=';', parse_dates=True) 

# Create waterlevel timeseries
dates1 = pd.DatetimeIndex(df1['longDATE'].values)
df1['wL#'] = df1['wL#'] / 100   # convert waterlevel (wL) from [cm] to [m]
wL = pd.Series(df1['wL#'].values, index=dates1)

# Create airpressure timeseries
dates2 = pd.DatetimeIndex(df2['longDATE'].values)
#df2['aP'] = df2['aP'] * 10.19716    # convert airpressure (aP) from [bar] to [mH2O]
df2['aP'] = df2['aP'] / 100     # convert airpressure (aP) from [cmH2O] to [mH2O]
aP = pd.Series(df2['aP'].values, index=dates2)

if t_res is True:
    # Resample timeseries
    wL = wL.resample(sr, how=stat)
    aP = aP.resample(sr, how=stat)
    
if t_slice is True:
    # Slice timeseries
    wL = wL[slice_t1:slice_t2]
    aP = aP[slice_t1:slice_t2]

# Correct waterlevel with airpressure data
wL = wL - aP

#######################
### Event detection ###
#######################
# Detect drawdown events in waterlevel timeseries and return h1, h2 and dt
events = helper.drawdown_event_detection(wL, peak_lookahead)
print events

##########################################################
### Compute transmissivity (T) (Cooper and Jacob 1946) ###
##########################################################
# TODO: this part is not yet correct i think
for e in range(events.shape[0]):
    # Get all waterlevels within a drawdown event
    wL_e = wL[events['h1_t'][e]:events['h2_t'][e]]
    print wL_e
    
    # Get delta t (dt) and delta waterlevel (dwL) (=single drawdown values)
    dwL = []
    dt = []
    for i in range(len(wL_e)):
        if i != 0:
            dwL.append(wL_e.values[i-1] - wL_e.values[i])
            dt.append(wL_e.index[i] - wL_e.index[i-1])
    dwL_e = pd.DataFrame({'dwL' : dwL,
                          'dt' : np.cumsum(dt)})
    dwL_e['dt'] = dwL_e['dt'] / np.timedelta64(1, 'm')    # convert timedelta64 to float
    print dwL_e
    
    # Fit straight line to the drawdown values on semilog scale
    x = dwL_e['dt']
    x_ln = np.log(x)
    y = dwL_e['dwL']
    a, b = np.polyfit(x_ln, y, 1)
    
    # Compute transmissivity (T) over one logarithmic cycle
    y_max = a * 1 + b
    y_min = a * 10 + b
    dsw = y_max - y_min   # [m]
    T = 2.303 * Q / (4 * np.pi * dsw)
    print T
    '''
    # Plot drawdown vs time on semilog scale
    plt.figure()
    #plt.plot(x_ln, a * x_ln + b, '--')
    x_plt = np.array(range(1,11,1)).reshape((10, 1))
    plt.plot(x_plt, a * x_plt + b, 'b-')
    plt.plot(x_ln, y, 'go', label='water level')
    #plt.plot(1, y_max, 'bo')
    #plt.plot(10, y_min, 'bo')
    plt.title('Pisciarelli Tennis Club (KELLER DCX22)')
    plt.xlabel('Time')
    plt.ylabel('Drawdown [m]')
    plt.legend(loc='best')
    #plt.xscale('log')
    plt.grid()
    plt.tight_layout()
    plt.savefig(wdir + 'pisc_wL_semilog_' + str(e) + '.png' , dpi=300)
    plt.show()
    plt.close()
    '''
####################
### Plot results ###
####################
fig = plt.figure()

ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()

# Plot waterlevel timeseries VS air pressure timeseries
if t_res is True:
    # Waterlevel
    wL_plt = wL.plot(ax=ax1, kind='line', color='blue', label = 'wL', 
            title='Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ')')
    plt_save = wdir + 'pisc_wL_aP_' + sr + stat + '.png'
else:
    # Waterlevel
    wL_plt = wL.plot(ax=ax1, kind='line', color='blue', label = 'wL', 
            title='Pisciarelli Tennis Club (KELLER DCX22)')   
    plt_save = wdir + 'pisc_wL_aP.png'  
wL_plt.set_xlabel('Time')
wL_plt.set_ylabel('wL [m]')
# Air pressure
aP_plt = aP.plot(ax=ax2, kind='line', color='red', label = 'aP',)
aP_plt.set_ylabel('aP [mH2O]')

ax1.legend(loc='lower left')
ax2.legend(loc='lower right')
plt.grid()
plt.tight_layout()
plt.savefig(plt_save, dpi=300)
plt.show()
plt.close()

# Plot waterlevel timeseries with event markers
if t_res is True:
    plt.plot(wL.index, wL.values, 'b-', label='water level')
    plt.title('Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ')')
    plt_save = wdir + 'pisc_wL_events_' + sr + stat + '.png'
else:
    plt.plot(wL.index, wL.values, 'b-', label='water level')
    plt.title('Pisciarelli Tennis Club (KELLER DCX22)')
    plt_save = wdir + 'pisc_wL_events.png'  
# Add event markers
plt.plot(events['h1_t'], events['h1'], 'r+', mew=2, ms=8, label='h1')
plt.plot(events['h2_t'], events['h2'], 'g+', mew=2, ms=8, label='h2')
#plt.xticks(rotation=15)
plt.xlabel('Time')
plt.ylabel('wL [m]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig(plt_save, dpi=300)
plt.show()
plt.close()

'''
plt.plot(dwL_e['dt'], dwL_e['dwL'], 'bo', label='water level')
# TODO: plot the fitted straight line
plt.plot(dwL_e['dt'], a * dwL_e['dt'] + b, 'r-')
plt.title('Pisciarelli Tennis Club (KELLER DCX22)')
plt.xlabel('Time [min]')
plt.ylabel('Drawdown [m]')
plt.legend(loc='best')
plt.xscale('log')
plt.grid()
plt.tight_layout()
plt.savefig(wdir + 'pisc_wL_semilogTEST.png' , dpi=300)
plt.show()
plt.close()
'''
# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'