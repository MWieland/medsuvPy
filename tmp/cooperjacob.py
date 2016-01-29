'''
---------------------------
    cooperjacob.py
---------------------------
Created on 26.01.2016
Last modified on 26.01.2016
Author: Marc Wieland
Description: Transmissivity (T) is estimated from the pumping rate (Q) and the change in drawdown 
             per log-cycle (Delta s) from the following equation:   T = 2.3 * Q / 4 * Pi * Delta s
             where, Delta s is the change in drawdown per log-cycle (L).
             
             
             r = distance from the discharging well
             t = time elapsed since start of discharge
             T = transmissibility of the aquifier
             S = coefficient of storage
             Q = discharge of the well
             
             
             
Reference: 
----
'''

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters to set ########################################################################################################
dat = 'FANG (DIVER) HH compensated.csv'
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_fang/marc/'
L = "H" # temporal resampling rate (hour, day, month, etc.)
stat = 'mean' # temporal resampling statistics 
Q = 600   # Pumping rate
############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

# read data from csv to pandas dataframe
df = pd.read_table(wdir + dat, header=0, delimiter=';', parse_dates=True) 

# create waterlevel timeseries
dates = pd.DatetimeIndex(df['longDATE'].values)
wL = pd.Series(df['wL'].values, index=dates)

# slice timeseries
#series = series['10/01/2013 00:00':'11/30/2013 00:00']

#######################################################
# compute transmissivity (T) per time interval (Cooper and Jacob 1946)
# TODO: where do i get Q from? how is ds computed? how do i pick the events?
ds = np.abs(wL.resample(L, how='min') - wL.resample(L, how='max'))
T = 2.303 * Q / 4 * np.pi * ds
#print T

# compute hydraulic conductivity (K) per time interval (Cooper and Jacob 1946)
#dt = 30
#r = 5
#K = r / (4 * np.pi * dt) * 2.303 * np.log(wL.resample(L, how='max') / wL.resample(L, how='min'))
#print K

#######################################################
# plot waterlevel timeseries
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
wL_plt = wL.resample(L, how=stat).plot(ax=ax1, kind='line', color='blue',
                                       label='water level', title='Solfatara Fangaia (ceraDIVER - ' + L + ' ' + stat + ')')
#wL_plt = wL.plot(ax=ax1, kind='line', color='black',
#                 label='water level', title='Solfatara Fangaia (ceraDIVER)') #logx=True)

#print wL
#print np.log(wL)


# TEST: get events identified by Heiko
events = pd.Series(df['event'].values, index=dates).dropna()
print events
# add event starts to plots
wL_plt.axvline(events.index[0], color='r', label='event start')
wL_plt.axvline(events.index[1], color='r')
wL_plt.axvline(events.index[2], color='r')
wL_plt.axvline(events.index[3], color='r')
wL_plt.axvline(events.index[4], color='r')
wL_plt.axvline(events.index[5], color='r')
# TODO: add event ends to plots
#wL_plt.axvline(events.index[0] + pd.DateOffset(days=12), color='g', label='event end')
#wL_plt.axvline(events.index[1] + pd.DateOffset(days=12), color='g')
#wL_plt.axvline(events.index[2] + pd.DateOffset(days=12), color='g')
#wL_plt.axvline(events.index[3] + pd.DateOffset(days=12), color='g')
#wL_plt.axvline(events.index[4] + pd.DateOffset(days=12), color='g')
#wL_plt.axvline(events.index[5] + pd.DateOffset(days=12), color='g')

wL_plt.set_xlabel('Time')
wL_plt.set_ylabel('wL_comp [cm]')
wL_plt.legend(loc='best')

# plot transmissivity timeseries
#ax2 = fig.add_subplot(3,1,2)
#T_plt = T.plot(ax=ax2, kind='line', title='Transmissivity (' + L + ')')
#T_plt.set_xlabel('Time')
#T_plt.set_ylabel('T [cm/' + L + ']')

# plot hydraulic conductivity timeseries
#ax3 = fig.add_subplot(3,1,3)
#K_plt = K.plot(ax=ax3, kind='line', title='Hydraulic conductivity (' + L + ')')
#K_plt.set_xlabel('Time')
#K_plt.set_ylabel('K [Darcy]')

plt.tight_layout()
plt.savefig(wdir + 'fang_wL_T_' + L + stat + '.png', dpi=300)
plt.show()
#######################################################

# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'