# -------------------------------------- Preprocessing --------------------------------------#
import librosa
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
from scipy.stats import stats
#from OpenBCI_lsl import *
import numpy as np
from pylab import *

def butter_bandpass_filter(data, lowcut, highcut, fs, order):   # IIR
    nyq = 0.5 * fs
    low = lowcut/nyq
    high = highcut/nyq
    b, a = signal.butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

def FIR_filter(data, lowcut, highcut, fs, order):   # FIR
    nyq = 0.5 * fs
    b = signal.firwin(order, [lowcut, highcut], window='hamming', pass_zero='bandpass', nyq = nyq)   # filter design
    y = filtfilt(b, 1.0, data)
    return y

def Preproccessing(win,srate, low, high, order):

    #### Re-reference ####
    for t in range(0, len(win.T)):  # win : channel by samples
        re = np.mean(win[:, t])
        win[:, t] = win[:, t] - re

    #### Filtering ####
    win = FIR_filter(win, low, high, srate, order)

    #### Resampling ####
    win = librosa.resample(win, srate, 64)

    ## Z-scoring
    win = (win - win.mean()) / win.std()

    return win


