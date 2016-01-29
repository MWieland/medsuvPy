if __name__ == "__main__":
    import numpy as np
    from peakdetect import peakdetect
    from math import pi
    import pylab
    
    i = 10000
    x = np.linspace(0,3.7*pi,i)
    y = (0.3*np.sin(x) + np.sin(1.3 * x) + 0.9 * np.sin(4.2 * x) + 0.06 * 
    np.random.randn(i))
    y *= -1
    
    _max, _min = peakdetect(y, x, 750, 0.30)
    xm = [p[0] for p in _max]
    ym = [p[1] for p in _max]
    xn = [p[0] for p in _min]
    yn = [p[1] for p in _min]
    
    plot = pylab.plot(x, y)
    pylab.hold(True)
    pylab.plot(xm, ym, 'r+')
    pylab.plot(xn, yn, 'g+')
    
    
    pylab.show()