'''
---------------------------
    pisciarelli.py
---------------------------
Created on 08.02.2016
Last modified on 22.02.2016
Author: Marc Wieland
Description: Single-well pump test to derive transmissivity and storage coefficient from waterlevel measurements.
             0. Convert units of input data to [m]
             1. Resample and/or slice timeseries (optional)
             2. Correct waterlevel with airpressure data (optional)
             3. Convert water level to [m below ground]
             4. Auto-detect events in waterlevel timeseries.
             5. Separate waterlevels for each event between drawdown and fill
             6. Interpolate drawdown and fill parts separately (optional)
             7. Compute Transmissivity (T) in [m^2/s]:
                 T = Q / (4 * np.pi * a)
                and Storage coefficient (S) in [m] according to Zheng et al.:
                 S = -4 * T * b / a * (r * r)
                 Q = pumping rate
                 a = slope of straight line fit
                 b = intercept of straight line fit
                 r = radius
             8. Write events, waterlevels drawdown and waterlevels fill into csv files and plot figures  
Reference: Zheng, Guo and Lei (2005): An improved straight-line fitting method for analyzing pumping 
            test recovery data. Ground water, 43 (6), 939-942.
           Cooper and Jacob (1946): A generalized method for evaluating formation constants
            and summarizing well-field history. Transactions of the American Geophysical
            Union, 27 (5), 526-534.
----
'''

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeseries.helper as helper

# Parameters to set ########################################################################################################
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_pitc/marc/' # Work directory
dat_wl = 'PITC (DCX22).csv' # Water level data [mbar]
###
c_ap = True # Correct water level with air pressure data
dat_bar = 'AGN3baro.csv'    # Air pressure data [cmH2O]
###
t_slice = True    # Temporal slicing
slice_t1 = '2014-07-01 07:00'   # Slicing start time
slice_t2 = '2014-07-01 15:00'   # Slicing end time
###
t_res = True    # Temporal resampling 
sr = "10min"    # Temporal resampling rate (e.g., 5min, H, D, M)
stat = 'mean'   # Temporal resampling statistics 
###
t_int = False   # Temporal interpolation 
ir = "1min"     # Temporal interpolation rate (e.g., 5min, H, D, M)
if_d = 'cubic'  # Temporal interpolation function for drawdown events 
if_f = 'cubic'  # Temporal interpolation function for fill events                
                # {'linear', 'time', 'index', 'values', 'nearest', 'zero', 'slinear', , 'quadratic',
                #  'cubic', 'barycentric', 'krogh', 'polynomial', 'spline', 'piecewise_polynomial', 'pchip'}
###
peak_lookahead = 1  # Event detection: distance to look ahead from a peak candidate to determine the actual peak
###
Q = 0.006944    # Pumping rate [cubic meters per second]
r = 0.5         # Radius [m]
g = 33.7        # Ground level [m] (used to compute meters below ground)
############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

####################################
### Data input and preprocessing ###
####################################
# Create waterlevel timeseries
df1 = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 
dates1 = pd.DatetimeIndex(df1['longDATE'].values)
df1['wL#'] = df1['wL#'] * 0.0101972   # convert waterlevel (wL) from [mbar] to [m]
wL = pd.Series(df1['wL#'].values, index=dates1)
if t_slice is True:
    # Slice waterlevel timeseries
    wL = wL[slice_t1:slice_t2]
if t_res is True:
    # Resample waterlevel timeseries
    wL = wL.resample(sr, how=stat)

if c_ap is True:
    # Create airpressure timeseries
    df2 = pd.read_table(wdir + dat_bar, header=0, delimiter=';', parse_dates=True) 
    dates2 = pd.DatetimeIndex(df2['longDATE'].values)
    df2['aP'] = df2['aP'] / 100     # convert airpressure (aP) from [cmH2O] to [mH2O]
    aP = pd.Series(df2['aP'].values, index=dates2)
    if t_slice is True:
        # Slice airpressure timeseries
        aP = aP[slice_t1:slice_t2]
    if t_res is True:
        # Resample airpressure timeseries
        aP = aP.resample(sr, how=stat)
    # Correct waterlevel with airpressure data
    wL = wL - aP

# Convert waterlevel to meters below ground
wL = wL - g

#######################
### Event detection ###
#######################
# Detect events in waterlevel timeseries and return h1, h2 and dt
events = helper.drawdown_event_detection(wL, peak_lookahead)
print events

###############################################################
### Get waterlevels for drawdown and fill around each event ###
###############################################################
wL_d = None
wL_f = None

for e in range(events.shape[0]):
    if e == 0:
        # Get waterlevels within first drawdown event
        wL_d = wL[events['h1_t'][e]:events['h2_t'][e]]
        if t_int is True:
            # Interpolate timeseries
            wL_d = wL_d.resample(ir, how='mean')    
            wL_d = wL_d.interpolate(method=if_d)
        #wL_d.plot() 
        wL_d = wL_d.to_frame()
        wL_d['event'] = e 
        if e != events.shape[0]-1:
            # Get waterlevels within first fill event
            wL_f = wL[events['h2_t'][e]:events['h1_t'][e+1]]
            if t_int is True:
                # Interpolate timeseries
                wL_f = wL_f.resample(ir, how='mean')
                wL_f = wL_f.interpolate(method=if_f)
            #wL_f.plot()
            wL_f = wL_f.to_frame()
            wL_f['event'] = e
    else:
        # Get waterlevels within all other drawdown events
        wL_d_ = wL[events['h1_t'][e]:events['h2_t'][e]]
        if t_int is True:
            # Interpolate timeseries
            wL_d_ = wL_d_.resample(ir, how='mean')    
            wL_d_ = wL_d_.interpolate(method=if_d)       
        #wL_d_.plot()
        wL_d_ = wL_d_.to_frame()
        wL_d_['event'] = e 
        wL_d = wL_d.append(wL_d_)
        if e != events.shape[0]-1:
            # Get waterlevels within all other fill events
            wL_f_ = wL[events['h2_t'][e]:events['h1_t'][e+1]]
            if t_int is True:
                # Interpolate timeseries
                wL_f_ = wL_f_.resample(ir, how='mean')
                wL_f_ = wL_f_.interpolate(method=if_f)
            #wL_f_.plot()
            wL_f_ = wL_f_.to_frame()
            wL_f_['event'] = e 
            wL_f = wL_f.append(wL_f_)
if wL_d is not None:
    wL_d.columns = ['wL', 'event']
if wL_f is not None:
    wL_f.columns = ['wL', 'event']

#################################################################################
### Compute Transmissivity (T) and Storage coefficient (S) (Zheng et al 2005) ###
#################################################################################
# Define columns for T and S
events['T'] = 0
events['S'] = 0

# Get fill event ids
eid = wL_f['event'].unique()

for e in range(len(eid)):
    # Get tp for fill event (converted to float and unit in [s])
    tp = events[events.index == eid[e]]['dt'] / np.timedelta64(1, 's')
    tp = tp.values[0]
    
    # Get waterlevels within fill event
    wL_e = wL_f[wL_f['event'] == eid[e]]
    wL_e.columns = ['wL', 'event']
    
    # Get delta waterlevel (dwL) and delta t (dt)
    dwL = []
    dt = []
    for i in range(len(wL_e)):
        if i != 0:
            dwL.append(wL_e['wL'].values[i] - wL_e['wL'].values[i-1])
            dt.append(wL_e.index[i] - wL_e.index[i-1])
    
    # Compute cummulative dwL and dt (note: this would be "s'" and "t'" in Zheng et al.)
    dwL_e = pd.DataFrame({'dwL' : np.cumsum(dwL),
                          'dt' : np.cumsum(dt)})
    # Convert timedelta64 to float and unit in [s]
    dwL_e['dt'] = dwL_e['dt'] / np.timedelta64(1, 's')
    dwL_e['tp'] = tp
    #print dwL_e
    
    # Fit straight line
    x = dwL_e['dt'] * np.log(dwL_e['dt'] + dwL_e['tp'] / dwL_e['dt'])
    y = dwL_e['dt'] * dwL_e['dwL']
    
    a, b = np.polyfit(x, y, 1)
    
    # Compute transmissivity (T) in [m^2/s]
    T = Q / (4 * np.pi * a)
    events['T'][events.index == eid[e]] = T
    
    # Compute storage coefficient (S) in [m]
    S = -4 * T * b / a * (r * r) 
    events['S'][events.index == eid[e]] = S
    
    # Plot event water level
    plt.plot(wL_d[wL_d['event'] == eid[e]].index, wL_d[wL_d['event'] == eid[e]]['wL'].values, 'b-', label='water level')
    plt.plot(wL_d[wL_d['event'] == eid[e]].index, wL_d[wL_d['event'] == eid[e]]['wL'].values, 'bo')
    plt.plot(wL_f[wL_f['event'] == eid[e]].index, wL_f[wL_f['event'] == eid[e]]['wL'].values, 'g-', label='water level')
    plt.plot(wL_f[wL_f['event'] == eid[e]].index, wL_f[wL_f['event'] == eid[e]]['wL'].values, 'go')
    plt.title('Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ') - Event ' + str(eid[e]))
    plt.xlabel('Time')
    plt.ylabel('wL [m]')
    plt.grid()
    plt.tight_layout()
    plt.savefig(wdir + 'pisc_wL_event_' + str(eid[e]) + '.png', dpi=300)
    #plt.show()
    plt.close()
    
    # Plot event straight line fit
    plt.plot(x, a * x + b, '--')
    plt.plot(x, y, 'ro')
    plt.title('Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ') - Event ' + str(eid[e]))
    plt.annotate('T = ' + str(T) + ' m^2/s', xy=(1000, 1000), xytext=(1000, 1000))
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.tight_layout()
    plt.savefig(wdir + 'pisc_linefit_event_' + str(eid[e]) + '.png', dpi=300)
    #plt.show()
    plt.close()

####################
### Plot results ###
####################
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()

# Plot waterlevel timeseries with event markers and airpressure data
if t_res is True:
    title='Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ')'
    plt_save = wdir + 'pisc_wL_' + sr + stat + '.png'    
else:
    title='Pisciarelli Tennis Club (KELLER DCX22)'
    plt_save = wdir + 'pisc_wL.png'
ax1.plot(wL.index, wL.values, 'b-', label='wL')
# Add event markers
ax1.plot(events['h1_t'], events['h1'], 'r+', mew=2, ms=8, label='h1')
ax1.plot(events['h2_t'], events['h2'], 'g+', mew=2, ms=8, label='h2')
ax1.set_ylabel('wL [m]')
ax1.legend(loc='lower left')
if c_ap is True:
    # Add air pressure
    ax2.plot(aP.index, aP.values, 'r-', label='aP')
    ax2.set_ylabel('aP [mH2O]')
    ax2.legend(loc='lower right')
#plt.xticks(rotation=15)
plt.title(title)
plt.xlabel('Time')
plt.grid()
plt.tight_layout()
plt.savefig(plt_save, dpi=300)
plt.show()
plt.close()

# Reduce events to only those where T and S could be computed (so where we have drawdown and fill parts)
events = events[events['T'] != 0.0]

# Plot Transmissivity (T) over all events
plt.plot(events['h2_t'], events['T'], 'r-')
plt.plot(events['h2_t'], events['T'], 'ro')
plt.title('Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ')')
plt.xlabel('Time')
plt.ylabel('Transmissivity [m^2/s]')
plt.grid()
plt.tight_layout()
plt.savefig(wdir + 'pisc_transmissivity.png', dpi=300)
#plt.show()
plt.close()

# Plot Storage coefficient (S) over all events
plt.plot(events['h2_t'], events['S'], 'r-')
plt.plot(events['h2_t'], events['S'], 'ro')
plt.title('Pisciarelli Tennis Club (KELLER DCX22 - ' + sr + ' ' + stat + ')')
plt.xlabel('Time')
plt.ylabel('Storage coefficient [-]')
plt.grid()
plt.tight_layout()
plt.savefig(wdir + 'pisc_storagecoefficient.png', dpi=300)
#plt.show()
plt.close()

######################
### Results output ###
######################
print events

# Write events dataframe to csv
events.to_csv(wdir + 'pisc_events.csv', sep=';')

# Write water levels for drawdown and fill events to csv
if wL_d is not None:
    wL_d.columns = ['wL', 'event']
    wL_d.to_csv(wdir + 'pisc_events_wL_drawdown.csv', sep=';')
if wL_f is not None:
    wL_f.columns = ['wL', 'event']
    wL_f.to_csv(wdir + 'pisc_events_wL_fill.csv', sep=';')

# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'

'''
##########################################################
### Compute transmissivity (T) (Cooper and Jacob 1946) ###
##########################################################
for e in range(events.shape[0]):
    # TODO: this part is not yet correct i think
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
