# -------------------------------------- Preprocessing --------------------------------------#
import librosa
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
from scipy.stats import stats
#from OpenBCI_lsl import *
import numpy as np
from pylab import *

def FIR_filter(data, lowcut, highcut, fs, order):   # FIR
    b = signal.firwin(order, [lowcut, highcut], window='hamming', pass_zero='bandpass', fs = fs)   # filter design
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


