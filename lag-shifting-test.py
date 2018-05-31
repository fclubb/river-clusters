# Script for testing the shifting of profiles based on the lag.
# From an ipython notebook on correlation of time series at:
# https://currents.soest.hawaii.edu/ocn_data_analysis/_static/SEM_EDOF.html

# Our standard imports:
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# Access to many standard distributions:
import scipy.stats as ss

# for shifting
from scipy.ndimage.interpolation import shift

# test fft
from numpy.fft import fft, ifft, fft2, ifft2, fftshift

# cross correlation.
# make 2 time series of 500 points
def cross_correlation_using_numpy(y1, y2):

    # make an array with the lags
    lags = np.arange(-npts + 1, npts)
    ccov = np.correlate(y1 - y1.mean(), y2 - y2.mean(), mode='full')
    ccor = ccov / (npts * y1.std() * y2.std())

    maxlag = lags[np.argmax(ccor)]
    print ("max correlation at lag{}".format(maxlag))

    # shift y2 so to reduce the lag
    shifted_y2 = shift(y2, maxlag, cval=np.nan)

    plt.plot(x,y1,c='b')
    plt.plot(x,shifted_y2,c='r')
    plt.show()

def cross_correlation_using_fft(y1, y2):
    f1 = fft(y1)
    f2 = fft(np.flipud(y2))
    cc = np.real(ifft(f1 * f2))
    return fftshift(cc)

# shift &lt; 0 means that y starts 'shift' time steps before x # shift &gt; 0 means that y starts 'shift' time steps after x
def compute_shift(y1, y2):
    assert len(y1) == len(y2)
    c = cross_correlation_using_fft(y1, y2)
    assert len(c) == len(y1)
    zero_index = int(len(y1) / 2) - 1
    max_lag = int(zero_index - np.argmax(c))
    print max_lag

    # shift y2 so to reduce the lag
    shifted_y2 = shift(y2, max_lag, cval=np.nan)

    plt.plot(x,y1,c='b')
    plt.plot(x,shifted_y2,c='r')
    plt.show()
    return max_lag

npts = 500
x = np.linspace(0, 50, npts)
y1 = 5 * np.sin(x/2) + np.random.randn(npts)
y2 = 5 * np.cos(x/2) + np.random.randn(npts)
compute_shift(y1, y2)
