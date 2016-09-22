'''
Created on 13 Jun 2016
@author: mwieland

Results reported in Zheng et al.:
T=30.73 ft2/min (=2.85 m2/min)
S=0.0696

"The residual drawdown data were recorded in an
observation well located at r=100 ft from the pumping
well. The well was pumped td=800 min with a discharge
rate of Q=162.9 ft3/min, and the recovery period was also re-
corded for 800 min after the pump was turned off."

BUG FIXES:
NOTE1: use recovery data (= data after the pump was switched off)!
NOTE2: remove first row of data (zero delta value)
NOTE3: adjust formula for x -> add brackets!
NOTE4: keep units consistent throughout the script (ft / min OR m / min)
'''
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters to set ########################################################################################################
site_name = 'Observation Well (Recovery)'
instrument = 'USDI Reference Data'
###
wdir = '/home/mwieland/eclipse_workspace/python/medsuv/testdata/' # Work directory
dat_wl = 'usdi_observationwell_recovery.csv' # File that holds the water level data 
col_wl_dt = 'dt_min'    # Column that holds the observed delta t [min]
col_wl_dl = 'dwl_m' # Column that holds the observed delta water level [m]

Q = 4.61           # Pumping rate [m3/min]
r = 30.48          # Radius [m]
#Q = 162.9         # Pumping rate [ft3/min]
#r = 100           # Radius [ft]

############################################################################################################################

starttime=time.time()
print 'Starttime : ' + str(time.strftime("%H:%M:%S"))

df1 = pd.read_table(wdir + dat_wl, header=0, delimiter=';', parse_dates=True) 

# Get cummulative dwL and dt (note: this would be "s'" and "t'" in Zheng et al.)
dwL_e = pd.DataFrame({'dwL' : df1[col_wl_dl],
                      'dt' : df1[col_wl_dt],
                      'tp' : 800})

# Plot data
plt.plot(dwL_e['dt'], dwL_e['dwL'], '--')
plt.title(site_name + ' (' + instrument + ')')
plt.xlabel('dt')
plt.ylabel('dwL')
plt.grid()
plt.tight_layout()
plt.show()
plt.close()

# Remove first entry from dataframe (zero delta value)
dwL_e = dwL_e.ix[1:]
print dwL_e

# Fit straight line
x = dwL_e['dt'] * np.log((dwL_e['dt'] + dwL_e['tp']) / dwL_e['dt'])
y = dwL_e['dt'] * dwL_e['dwL']

a, b = np.polyfit(x, y, 1)
print 'a : ' + str(a)
print 'b : ' + str(b)

# Compute transmissivity (T) in [m^2/min]
T = Q / (4 * np.pi * a)

# Compute storage coefficient (S) in [m]
S = -(4 * T * b) / (a * (r * r))

print 'T : ' + str(T)
print 'S : ' + str(S)

# Plot event straight line fit
plt.plot(x, a * x + b, '--')
plt.plot(x, y, 'ro')
plt.title(site_name + ' (' + instrument + ')')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid()
plt.tight_layout()
plt.show()
plt.close()
