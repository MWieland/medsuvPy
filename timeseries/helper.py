'''
---------------------------
        helper.py
---------------------------
Created on 27.01.2016
Last modified on 28.01.2016
Author: Marc Wieland
Description: Contains helper functions related to timeseries analysis.
----
'''

#import sys
import numpy as np
import pandas as pd

i = 10000
x = np.linspace(0, 3.5 * np.pi, i)
y = (0.3*np.sin(x) + np.sin(1.3 * x) + 0.9 * np.sin(4.2 * x) + 0.06 *
    np.random.randn(i))


def _datacheck_peakdetect(x_axis, y_axis):
    if x_axis is None:
        x_axis = range(len(y_axis))
    
    if len(y_axis) != len(x_axis):
        raise (ValueError, 
                'Input vectors y_axis and x_axis must have same length')
    
    #needs to be a numpy array
    y_axis = np.array(y_axis)
    x_axis = np.array(x_axis)
    return x_axis, y_axis


def peakdetect(y_axis, x_axis = None, lookahead = 300, delta=0):
    """
    Detects local maxima and minima in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    arguments:
    y_axis: A list containg the signal over which to find peaks
    x_axis (optional): A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead (optional): distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta (optional): this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    returns: two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
    
    Taken from Python code at: 
    https://gist.github.com/sixtenbe/1178136#file-peakdetect-py
    
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    """
    max_peaks = []
    min_peaks = []
    dump = []   #Used to pop the first hit which almost always is false
       
    # check input data
    x_axis, y_axis = _datacheck_peakdetect(x_axis, y_axis)
    # store data length for later use
    length = len(y_axis)
    
    
    #perform some checks
    if lookahead < 1:
        raise ValueError, "Lookahead must be '1' or above in value"
    if not (np.isscalar(delta) and delta >= 0):
        raise ValueError, "delta must be a positive number"
    
    #maxima and minima candidates are temporarily stored in
    #mx and mn respectively
    mn, mx = np.Inf, -np.Inf
    
    #Only detect peak if there is 'lookahead' amount of points after it
    for index, (x, y) in enumerate(zip(x_axis[:-lookahead], 
                                        y_axis[:-lookahead])):
        if y > mx:
            mx = y
            mxpos = x
        if y < mn:
            mn = y
            mnpos = x
        
        ####look for max####
        if y < mx-delta and mx != np.Inf:
            #Maxima peak candidate found
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].max() < mx:
                max_peaks.append([mxpos, mx])
                dump.append(True)
                #set algorithm to only find minima now
                mx = np.Inf
                mn = np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
                continue
            #else:  #slows shit down this does
            #    mx = ahead
            #    mxpos = x_axis[np.where(y_axis[index:index+lookahead]==mx)]
        
        ####look for min####
        if y > mn+delta and mn != -np.Inf:
            #Minima peak candidate found 
            #look ahead in signal to ensure that this is a peak and not jitter
            if y_axis[index:index+lookahead].min() > mn:
                min_peaks.append([mnpos, mn])
                dump.append(False)
                #set algorithm to only find maxima now
                mn = -np.Inf
                mx = -np.Inf
                if index+lookahead >= length:
                    #end is within lookahead no more peaks can be found
                    break
            #else:  #slows shit down this does
            #    mn = ahead
            #    mnpos = x_axis[np.where(y_axis[index:index+lookahead]==mn)]
    
    #Remove the false hit on the first value of the y_axis
    try:
        if dump[0]:
            max_peaks.pop(0)
        else:
            min_peaks.pop(0)
        del dump
    except IndexError:
        #no peaks were found, should the function return empty lists?
        pass
        
    return [max_peaks, min_peaks]


def rolling_window(a, size):
    shape = a.shape[:-1] + (a.shape[-1] - size + 1, size)
    strides = a.strides + (a. strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)


def drawdown_event_detection(timeseries, peak_lookahead):
    """
    Performs drawdown event detection in waterlevel timeseries using
    a standard peak detection function.
    
    arguments:
    timeseries: Pandas timeseries object that holds the waterlevel as values.
    peak_lookahead: Distance to look ahead from a peak candidate to determine if it is the actual peak.
                    '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    
    returns: Pandas dataframe with h1: waterlevel at event start time.
                                   h1_t: event start time.
                                   h2: waterlevel at event end time.
                                   h2_t: event end time.
                                   dt: timedifference between start and end.                                
    """
    # detect max and min peaks in timeseries
    peak_max, peak_min = peakdetect(timeseries.values, timeseries.index, lookahead=peak_lookahead, delta=0.30)
    if len(peak_max) == 0 or len(peak_min) == 0:
        print 'Sorry no peaks found in timeseries. Try with a smaller peak_lookahead value.'
        exit()
    else:
        # convert resulting lists to timeseries
        peak_max = np.array(peak_max)
        peak_min = np.array(peak_min)
        peak_max_dates = pd.DatetimeIndex(peak_max[:,0])
        peak_min_dates = pd.DatetimeIndex(peak_min[:,0])
        peak_max = pd.Series(peak_max[:,1], index=peak_max_dates)
        peak_min = pd.Series(peak_min[:,1], index=peak_min_dates)
        
        # create helper timeseries to identify drawdown event starts and ends from peaks
        peak_max_ = np.zeros(shape=(peak_max.shape[0],1), dtype=int)
        peak_min_ = np.zeros(shape=(peak_min.shape[0],1), dtype=int)
        peak_max_.fill(1)   # identifier for event starts
        peak_min_.fill(2)   # identifier for event ends
        peak_max_ = pd.Series(peak_max_[:,0], index=peak_max_dates)
        peak_min_ = pd.Series(peak_min_[:,0], index=peak_min_dates)
        peaks_ = peak_max_.combine_first(peak_min_)
        peaks_ = np.all(rolling_window(peaks_.values, 2) == [1, 2], axis=1)
        
        # extract event-relevant peaks from peaks timeseries using helper object
        peaks = peak_max.combine_first(peak_min)
        h1 = []
        h1_t = []
        h2 = []
        h2_t = []
        for i in range(len(peaks_)):
            if peaks_[i] == 1:
                # if helper object is True it is an event
                h1.append(peaks.values[i])
                h1_t.append(peaks.index[i])
                h2.append(peaks.values[i+1])
                h2_t.append(peaks.index[i+1])
        
        # combine peaks into events object
        events = {'h1' : h1,
                  'h1_t' : h1_t,
                  'h2' : h2,
                  'h2_t' : h2_t}
        events = pd.DataFrame(events)
        
        # compute time difference between events start and end    
        events['dt'] = events['h2_t'] - events['h1_t']
        
        return events