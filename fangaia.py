'''
---------------------------
        fangaia.py
---------------------------
Created on 27.01.2016
Last modified on 29.01.2016
Author: Marc Wieland
Description: Hydraulic refill test to derive conductivity from waterlevel measurements.
             0. Resample and/or slice waterlevel timeseries (optional).
             1. Auto-detect drawdown events in waterlevel timeseries.
             2. Compute Hydraulic conductivity (K) as follows: 
                 K= r / 4 dt * 2.303 * log(h1h2) 
                 r = radius (constant from field  experiment)
                 h1: waterlevel at start of an event (after refill)
                 h2: waterlevel at end of an event (after drawdown)
                 dt: time difference between h1 and h2
             3. Correct K with external evapotranspiration data.
Reference: 
----
'''

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import timeseries.helper as helper

# Parameters to set ###############################################################################################
dat_wl = 'FANG (DIVER) HH compensated.csv'  # waterlevel data
dat_ev = 'FANG EP.csv'  # evapotranspiration data
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_fang/marc/' # work directory
r = 5   # radius [m] (constant from field  experiment)
###
t_res = False # temporal resampling 
sr = "D" # temporal resampling rate (e.g., 5min, H, D, M)
stat = 'mean' # temporal resampling statistics 
###
t_slice = False    # temporal slicing
slice_t1 = '10/01/2013 00:00'   # slicing start time
slice_t2 = '11/30/2013 00:00'   # slicing end time
###
peak_lookahead = 200 #200     # distance to look ahead from a peak candidate to determine if it is the actual peak
###################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

####################################
### Data input and preprocessing ###
####################################
# Read waterlevel data from csv to pandas dataframe
df_wl = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 

# Read evaporation data from csv to pandas dataframe
df_ev = pd.read_table(wdir + dat_ev, header=0, delimiter=';') 

# Create waterlevel timeseries
wL = pd.Series(df_wl['wL'].values, index=pd.DatetimeIndex(df_wl['longDATE'].values))

if t_res is True:
    # Resample waterlevel timeseries
    wL = wL.resample(sr, how=stat)
    
if t_slice is True:
    # Slice waterlevel timeseries
    wL = wL[slice_t1:slice_t2]

#######################
### Event detection ###
#######################
# Detect drawdown events in waterlevel timeseries and return h1, h2 and dt
events = helper.drawdown_event_detection(wL, peak_lookahead)

# Convert dt from [ns] to [s]
dt = np.array((events['dt'].values / 1000000000), dtype=np.float64)

##############################
### Hydraulic conductivity ###
##############################
# Compute hydraulic conductivity (K) and convert from [m/s] to [Darcy]
events['K'] = r / (4 * np.pi * dt) * 2.303 * np.log(events['h1'] / events['h2'])
events['K'] = events['K'] * 100000

########################################
### Corrected hydraulic conductivity ###
########################################
# Lookup monthly evapotranspiration values (ETpot) for the event start (h1_t) and end months (h2_t)
events['et'] = events['h2_t'].dt.month.replace(range(1,13), df_ev['ETpot'].tolist())

# Correct h2 waterlevel with ET values for each event 
# TODO: check formula with Heiko (this is one to one what he did in excel)
dtd = dt / 86400 # convert dt from [s] to [days] and c
events['h2corr'] = events['h2'] + events['et'] / 1000 * dtd / 30.25

# Compute hydraulic conductivity (K) and convert from [m/s] to [Darcy]
events['Kcorr'] = r / (4 * np.pi * dt) * 2.303 * np.log(events['h1'] / events['h2corr'])
events['Kcorr'] = events['Kcorr'] * 100000
print events

# TODO: redo the correction more precisely
# Get start and end month and read the according ET values from lookup table
#events['et_h1_t'] = events['h1_t'].dt.month.replace(range(1,13), df_ev['ETpot'].tolist())
#events['et_h2_t'] = events['h2_t'].dt.month.replace(range(1,13), df_ev['ETpot'].tolist())

# Compute ET for the duration of the event (dt) - taking into account the number of days per month
#events['et'] = 

# Correct h2 waterlevel with ET values for each event 
#events['h2corr'] = events['h2'] + events['et']

# Compute hydraulic conductivity (K) and convert from [m/s] to [Darcy]
#events['Kcorr'] = r / (4 * np.pi * dt) * 2.303 * np.log(events['h1'] / events['h2corr'])
#events['Kcorr'] = events['Kcorr'] * 100000

####################
### Plot results ###
####################
# Plot waterlevel timeseries
fig = plt.figure()
if t_res is True:
    plt.plot(wL.index, wL.values, 'b-', label='water level')
    plt.title('Solfatara Fangaia (ceraDIVER - ' + sr + ' ' + stat + ')')
    plt_save = wdir + 'fang_wL_T' + sr + stat + '.png'
else:
    plt.plot(wL.index, wL.values, 'b-', label='water level')
    plt.title('Solfatara Fangaia (ceraDIVER)')
    plt_save = wdir + 'fang_wL_T.png'  
# Add event markers
plt.plot(events['h1_t'], events['h1'], 'r+', mew=2, ms=8, label='h1')
plt.plot(events['h2_t'], events['h2'], 'g+', mew=2, ms=8, label='h2')
plt.xticks(rotation=15)
plt.xlabel('Time')
plt.ylabel('wL_comp [cm]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig(plt_save, dpi=300)
plt.show()
plt.close()

# Plot hydraulic conductivity for detected drawdown events
if t_res is True:
    plt.plot(events['h1_t'], events['K'], 'ro', mew=2, ms=8, label='hyd. cond.')
    plt.plot(events['h1_t'], events['Kcorr'], 'go', mew=2, ms=8, label='hyd. cond. corr.')
    plt.title('Solfatara Fangaia (ceraDIVER - ' + sr + ' ' + stat + ')')
    plt_save = wdir + 'fang_K_T' + sr + stat + '.png'
else:
    plt.plot(events['h1_t'], events['K'], 'ro', mew=2, ms=8, label='hyd. conductivity')
    plt.plot(events['h1_t'], events['Kcorr'], 'go', mew=2, ms=8, label='hyd. cond. corr.')
    plt.title('Solfatara Fangaia (ceraDIVER)')
    plt_save = wdir + 'fang_K_T.png'  
plt.xticks(rotation=15)
plt.xlabel('Time')
plt.ylabel('K [Darcy]')
plt.legend(loc='best')
plt.grid()
plt.tight_layout()
plt.savefig(plt_save, dpi=300)
plt.show()
plt.close()

######################
### Results output ###
######################
# write events dataframe to csv
if t_res is True:
    events.to_csv(wdir + 'fang_events_' + sr + stat + '.csv', sep=';')
else:
    events.to_csv(wdir + 'fang_events.csv', sep=';')

# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'