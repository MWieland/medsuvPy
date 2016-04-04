'''
Created on Apr 1, 2016
# TODO: I think that the input that Heiko sent me is not Waterlevels but delta waterlevels (=residual drawdowns=s)
#       This would change the whole game and not allow my script to process it properly -> here is a test where I
#       input directly the values as delta waterlevel (dwL) and delta time (dt)
#       -> the result for T is pretty close to the one reported in the paper, but S is still off...
@author: marc
'''
#TODO: test method against original groundwater manual data 1977 as in method paper of zheng et al.


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

Q = 0.05        # Pumping rate [cubic meters per second]
r = 50          # Radius [m]
tp = 4200       # Pumping time [s]

# Read waterlevel (wL) [m] and delta t (dt) [min]
wL = [0.593,0.506,0.415,0.354,0.309,0.276,0.249,0.208,0.178,0.138,0.111,0.085,0.068,0.045,0.032,0.023]
dt_ = [0,5,10,15,20,25,30,40,50,70,90,120,150,210,270,330]

# convert dt in [s]
dt_ = np.multiply(dt_,60)

# Get delta waterlevel (dwL)
dwL = []
dt = []
x=0
for i in range(len(wL)):
    if i != 0:
        x += 1
        dwL.append(wL[i] - wL[i-1])
        dt.append(dt_[x])

# Compute cummulative dwL and dt (note: this would be "s'" and "t'" in Zheng et al.)
# TODO: not sure about this part here if i really need to use the cummulative here -> the axis values of the
#       resulting plot seems extremely high
#       this may not affect the a value (and thus T) but it influences the b value (and thus S)!
dwL_e = pd.DataFrame({'dwL' : np.cumsum(np.abs(dwL)),
                      'dt' : dt})
                      
dwL_e['tp'] = tp
print dwL_e

# Fit straight line
x = dwL_e['dt'] * np.log(dwL_e['dt'] + dwL_e['tp'] / dwL_e['dt'])
y = dwL_e['dt'] * dwL_e['dwL']

#x = x.fillna(0)
print x
print y

a, b = np.polyfit(x, y, 1)

print a
print b
  
# Compute transmissivity (T) in [m^2/s]
T = Q / (4 * np.pi * a)

# Compute storage coefficient (S) in [m]
# TODO: needed to adjust this equation - old one was wrong (brackets!)
S = -(4 * T * b) / (a * (r * r))

print T
print S


# Plot event straight line fit
plt.plot(x, a * x + b, '--')
plt.plot(x, y, 'ro')
plt.annotate('T = ' + str(T) + ' m^2/s', xy=(1000, 1000), xytext=(1000, 1000))
plt.xlabel('X')
plt.ylabel('Y')
plt.grid()
plt.tight_layout()
#plt.savefig(wdir + 'linefit_event_' + str(eid[e]) + '.png', dpi=300)
plt.show()
plt.close()

