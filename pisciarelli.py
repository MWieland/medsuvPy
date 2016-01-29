'''
---------------------------
    pisciarelli.py
---------------------------
Created on 26.01.2016
Last modified on 26.01.2016
Author: Marc Wieland
Description:      
Reference: Transmissivity (T) is estimated from the pumping rate (Q) and the change in drawdown 
             per log-cycle (Delta s) from the following equation:   T = 2.3 * Q / 4 * Pi * Delta s
             where, Delta s is the change in drawdown per log-cycle (L).
             
             
             r = distance from the discharging well
             t = time elapsed since start of discharge
             T = transmissibility of the aquifier
             S = coefficient of storage
             Q = discharge of the well
             
TODO: so far only a few test plots and an example formula. to be finished        
----
'''

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters to set ########################################################################################################
dat_wl = 'PITC (DCX22).csv'
dat_bar = 'AGN3baro.csv'
wdir = '/media/datadrive_/projects/medsuv_heiko/exp_pitc/marc/'
L = "D" # temporal resampling rate (hour, day, month, etc.)
stat = 'mean' # temporal resampling statistics 
############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

# read data from csv to pandas dataframe
df1 = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 
df2 = pd.read_table(wdir + dat_bar, header=0, delimiter=';', parse_dates=True) 

# create waterlevel timeseries
dates1 = pd.DatetimeIndex(df1['longDATE'].values)
wL = pd.Series(df1['wL#'].values, index=dates1)

# create airpressure timeseries
dates2 = pd.DatetimeIndex(df2['longDATE'].values)
aP = pd.Series(df2['aP'].values, index=dates2)

# slice timeseries
#wL = wL['08/01/2014 00:00':'08/31/2014 23:59']
#aP = aP['08/01/2014 00:00':'08/31/2014 23:59']

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
wL_plt = wL.resample(L, how=stat).plot(ax=ax1, kind='line', color='red', label = 'wL',
                                       title='Pisciarelli Tennis Club (KELLER DCX22 - ' + L + ' ' + stat + ')')
#wL_plt = wL.plot(ax=ax1, kind='line', color='black',
#                 label='water level', title='Pisciarelli Tennis Club (KELLER DCX22)') #logx=True)
wL_plt.set_xlabel('Time')
wL_plt.set_ylabel('wL# [cm]')

# plot air pressure timeseries
ax2 = ax1.twinx()
aP_plt = aP.resample(L, how=stat).plot(ax=ax2, kind='line', color='blue', label = 'aP',)
aP_plt.set_ylabel('aP [bar]')

ax1.legend(loc='lower left')
ax2.legend(loc='lower right')

plt.tight_layout()
plt.savefig(wdir + 'pisc_wL_aP_T_' + L + stat + '.png', dpi=300)
plt.show()
#######################################################

# Get run time
endtime = time.time()
time_total = endtime-starttime
print str(time_total) + ' sec'