# -------------------------------------- Preprocessing --------------------------------------#
import librosa
from scipy import signal
from scipy.signal import butter, lfilter, resample, filtfilt
#from OpenBCI_lsl import *
import numpy as np

def butter_bandpass_filter(data, lowcut, highcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut/nyq
    high = highcut/nyq
    b, a = butter(order, [low, high], btype='band')
    y = filtfilt(b, a, data)
    return y

def Preproccessing(win,srate, low, high, order):

    #### Re-reference ####
    for t in range(0, srate * 15):
        re = np.mean(win[:, t])
        win[:, t] = win[:, t] - re

    #### Filtering ####

    win = butter_bandpass_filter(win, low, high, srate, order)  # order 조정필요.

    #### Resampling ####

    win = librosa.resample(win, srate, 64)

    ## Z-scoring

    win = (win - win.mean()) / win.std()

    return win