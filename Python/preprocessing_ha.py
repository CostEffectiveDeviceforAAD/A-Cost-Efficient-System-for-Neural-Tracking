# -------------------------------------- Preprocessing --------------------------------------#
import librosa
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
#from OpenBCI_lsl import *
import numpy as np
from pylab import *

def FIR_filter(data, lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut/nyq
    high = highcut/nyq
    #b, a = butter(order, [low, high], btype='band')
    b = signal.firwin(order, [low, high], window='hamming', pass_zero='bandpass')
    y = lfilter(b, 1.0, data)
    return y

def Preproccessing(win,srate, low, high, order):

    #### Re-reference ####
    for t in range(0, len(win.T)):
        re = np.mean(win[:, t])
        win[:, t] = win[:, t] - re

    #### Filtering ####

    win = FIR_filter(win, low, high, srate, order)  # order 조정필요.

    #### Resampling ####

    win = librosa.resample(win, srate, 64)

    ## Z-scoring

    win = (win - win.mean()) / win.std()

    return win

#Plot frequency and phase response
def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * log10 (abs(h))
    subplot(211)
    plot(w/max(w),h_dB)
    ylim(-150, 5)
    ylabel('Magnitude (db)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Frequency response')
    subplot(212)
    h_Phase = unwrap(arctan2(imag(h),real(h)))
    plot(w/max(w),h_Phase)
    ylabel('Phase (radians)')
    xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    title(r'Phase response')
    subplots_adjust(hspace=0.5)